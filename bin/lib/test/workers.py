"""test:workers — per-tag worker test runner with on-demand fixture downloads."""

import argparse
import os
import shutil
import subprocess
import sys
import tarfile
import urllib.request

from lib._docker import bake_target, require_drdb
from lib._runtime import REPO_ROOT, Globals, run, stderr
from lib.test._shared import coverage_command

# S3 bucket holding worker test fixtures (tarballs, fastq/sra/cel/PCL files).
WORKER_TEST_DATA_REPO = "https://s3.amazonaws.com/data-refinery-test-assets"

# Tag → bake target. agilent reuses affymetrix's image; qn and janitor reuse
# smasher's. All other tags map to their own image.
WORKER_TAG_TO_IMAGE = {
    "salmon": "salmon",
    "transcriptome": "transcriptome",
    "no_op": "no_op",
    "downloaders": "downloaders",
    "smasher": "smasher",
    "illumina": "illumina",
    "agilent": "affymetrix",
    "affymetrix": "affymetrix",
    "qn": "smasher",
    "janitor": "smasher",
    "compendia": "compendia",
}

# Default execution order, mirroring the legacy script's `worker_images` loop.
WORKER_TAGS = [
    "salmon",
    "transcriptome",
    "no_op",
    "downloaders",
    "smasher",
    "illumina",
    "agilent",
    "affymetrix",
    "qn",
    "janitor",
    "compendia",
]


def _fetch_url(url, dest):
    """Download `url` to `dest`, creating parent dirs as needed."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    if Globals.verbose:
        stderr(f"+ download {url} → {dest}")
    if Globals.dry_run:
        return
    urllib.request.urlretrieve(url, dest)


def _fetch_if_missing(name, dest):
    """Download `<repo>/<name>` to `dest` only if `dest` doesn't already exist."""
    if dest.exists():
        return
    stderr(f"Downloading {name}")
    _fetch_url(f"{WORKER_TEST_DATA_REPO}/{name}", dest)


def _fetch_tarball(name, extract_to):
    """Download `<repo>/<name>` to a temp file, extract into `extract_to`, then delete the tar."""
    extract_to.mkdir(parents=True, exist_ok=True)
    archive = extract_to / name
    stderr(f"Downloading + extracting {name}")
    _fetch_url(f"{WORKER_TEST_DATA_REPO}/{name}", archive)
    if Globals.dry_run:
        return
    with tarfile.open(archive) as tf:
        tf.extractall(extract_to)
    archive.unlink()


def _fetch_salmon_data(test_volume):
    salmon_dir = test_volume / "salmon_tests"
    # `newer` marker means the tarball was updated upstream — re-download.
    if not salmon_dir.exists() or (salmon_dir / "newer").exists():
        if salmon_dir.exists():
            shutil.rmtree(salmon_dir)
        _fetch_tarball("salmon_tests_newer.tar.gz", test_volume)

    # SalmonTools — wipe + re-extract every run (matches legacy).
    tools_dir = test_volume / "salmontools"
    if tools_dir.exists():
        shutil.rmtree(tools_dir)
    _fetch_tarball("salmontools_test_data.tar.gz", test_volume)

    rna_dir = test_volume / "raw/TEST/SALMON"
    for fname in ("ERR1562482_1.fastq.gz", "ERR1562482_2.fastq.gz", "ERR1562482.sra"):
        _fetch_if_missing(fname, rna_dir / fname)


def _fetch_affymetrix_data(test_volume):
    cel_dir = test_volume / "raw/TEST/CEL"
    for fname in (
        "GSM1426071_CD_colon_active_1.CEL",
        "GSM45588.CEL",
        "GSM1364667_U_110208_7-02-10_S2.CEL",
    ):
        _fetch_if_missing(fname, cel_dir / fname)
    pcl_dir = test_volume / "TEST/PCL"
    for fname in (
        "GSM1426071_CD_colon_active_1.PCL",
        "GSM45588.PCL",
        "GSM1364667_U_110208_7-02-10_S2.PCL",
    ):
        _fetch_if_missing(fname, pcl_dir / fname)


def _fetch_transcriptome_data(test_volume):
    tx_dir = test_volume / "raw/TEST/TRANSCRIPTOME_INDEX/AEGILOPS_TAUSCHII"
    for fname in ("aegilops_tauschii_short.fa.gz", "aegilops_tauschii_short.gtf.gz"):
        _fetch_if_missing(fname, tx_dir / fname)
    _fetch_if_missing(
        "Homo_sapiens_testdata.gtf",
        test_volume / "raw/TEST/TRANSCRIPTOME_INDEX/Homo_sapiens_testdata.gtf",
    )


def _fetch_illumina_data(test_volume):
    ilu_dir = test_volume / "raw/TEST/ILLUMINA"
    for fname in (
        "GSE22427_non-normalized.txt",
        "GSE54661_non_normalized.txt",
        "GSE106321_non-normalized.txt",
        "GSE33814_trimmed_non-normalized.txt",
        "GSE112517_non-normalized.txt",
        "GSE48023_trimmed_non-normalized.txt",
        "GSE41355_non-normalized.txt",
        "GSE100301_non-normalized.txt",
    ):
        _fetch_if_missing(fname, ilu_dir / fname)
    _fetch_if_missing("Ad-Cre-2.AVG_Signal.tsv", ilu_dir / "reference/Ad-Cre-2.AVG_Signal.tsv")


def _fetch_agilent_data(test_volume):
    _fetch_if_missing(
        "GSM466597_95899_agilent.txt",
        test_volume / "raw/TEST/AGILENT_TWOCOLOR/GSM466597_95899_agilent.txt",
    )


def _fetch_no_op_data(test_volume):
    raw_dir = test_volume / "raw/TEST/NO_OP"
    for fname in (
        "GSM557500-tbl-1.txt",
        "GSM1234847_sample_table.txt",
        "GSM1089291-tbl-1.txt",
        "GSM1089291-tbl-1-modified.txt",
    ):
        _fetch_if_missing(fname, raw_dir / fname)

    # GSM1234847_sample_table_headerless.txt is derived locally: drop the header.
    headerless = raw_dir / "GSM1234847_sample_table_headerless.txt"
    if not headerless.exists():
        stderr("Generating GSM1234847_sample_table_headerless.txt")
        if not Globals.dry_run:
            with (raw_dir / "GSM1234847_sample_table.txt").open() as src, headerless.open(
                "w"
            ) as dst:
                for i, line in enumerate(src):
                    if i > 0:
                        dst.write(line)

    exp_dir = test_volume / "TEST/NO_OP/EXPECTED"
    for fname in (
        "gene_converted_GSM557500-tbl-1.txt",
        "GSM269747.PCL",
        "gene_converted_GSM1234847-tbl-1.txt",
        "gene_converted_GSM1089291-tbl-1.txt",
    ):
        _fetch_if_missing(fname, exp_dir / fname)


def _fetch_smasher_data(test_volume):
    """Smasher fixtures — also pulled in by the compendia tag."""
    pcl_dir = test_volume / "PCL"
    for fname in (
        "GSM1237810_T09-1084.PCL",
        "GSM1237812_S97-PURE.PCL",
        "GSM1238108-tbl-1.txt",
        "GSM1487313_liver.PCL",
        "SRP149598_gene_lengthScaledTPM.tsv",
        "GSM1084806-tbl-1.txt",
        "GSM1084807-tbl-1.txt",
        "SRR1731761_output_gene_lengthScaledTPM.tsv",
        "SRR1731762_output_gene_lengthScaledTPM.tsv",
        "danio_target.tsv",
    ):
        _fetch_if_missing(fname, pcl_dir / fname)
    bad_dir = test_volume / "BADSMASH"
    for fname in ("big.PCL", "small.PCL", "bad.PCL"):
        _fetch_if_missing(fname, bad_dir / fname)
    quant_dir = test_volume / "QUANT"
    for fname in ("smasher-test-quant.sf", "smasher-test-truncated-quant.sf"):
        _fetch_if_missing(fname, quant_dir / fname)


def _fetch_qn_data(test_volume):
    qn_dir = test_volume / "QN"
    for fname in ("1.tsv", "2.tsv", "3.tsv", "4.tsv", "5.tsv", "6.tsv", "7.tsv"):
        _fetch_if_missing(fname, qn_dir / fname)


def _fetch_compendia_data(test_volume):
    # Compendia tests reuse the smasher fixtures, plus their own list-file driven downloads.
    _fetch_smasher_data(test_volume)
    for kind, list_file in (("MICROARRAY", "microarray.txt"), ("RNASEQ", "rnaseq.txt")):
        list_dir = test_volume / f"raw/TEST/{kind}"
        if (list_dir / list_file).exists():
            continue
        list_dir.mkdir(parents=True, exist_ok=True)
        src_list = REPO_ROOT / "workers/tests/data" / list_file
        shutil.copy(src_list, list_dir / list_file)
        stderr(f"Downloading {kind} files listed in {list_file}")
        if Globals.dry_run:
            continue
        with (list_dir / list_file).open() as f:
            for url in f:
                url = url.strip()
                if not url:
                    continue
                dest = list_dir / url.rsplit("/", 1)[-1]
                _fetch_url(url, dest)
    _fetch_if_missing("danio_target.tsv", test_volume / "QN/danio_target.tsv")


WORKER_DATA_FETCHERS = {
    "salmon": _fetch_salmon_data,
    "affymetrix": _fetch_affymetrix_data,
    "transcriptome": _fetch_transcriptome_data,
    "illumina": _fetch_illumina_data,
    "agilent": _fetch_agilent_data,
    "no_op": _fetch_no_op_data,
    "smasher": _fetch_smasher_data,
    "qn": _fetch_qn_data,
    "compendia": _fetch_compendia_data,
    # downloaders/janitor have no dedicated test data block in the legacy script.
}


def _ensure_worker_test_data(tag, test_volume):
    """Run every fetcher whose tag matches `tag` (or all when tag is None)."""
    test_volume.mkdir(parents=True, exist_ok=True)
    if tag is None:
        for fetcher in WORKER_DATA_FETCHERS.values():
            fetcher(test_volume)
    elif tag in WORKER_DATA_FETCHERS:
        WORKER_DATA_FETCHERS[tag](test_volume)
    # Compendia explicitly inherits smasher data even when run alone — fetcher above handles that.


def _run_worker_image(image_name, tag, forwarded_args, test_volume):
    """Run the coverage-wrapped tests inside the worker image via docker run."""
    extra = ["--tag", tag] + list(forwarded_args)
    cmd = [
        "docker",
        "run",
        "--rm",
        "--network",
        "refinebio_default",
        "--env",
        "AWS_ACCESS_KEY_ID",
        "--env",
        "AWS_SECRET_ACCESS_KEY",
        "--env-file",
        "workers/environments/test",
        "--memory=5G",
        "--platform",
        "linux/amd64",
        "--volume",
        f"{test_volume}:/home/user/data_store",
    ]
    cmd.append("--tty")
    if sys.stdout.isatty():
        cmd.append("--interactive")
    image_tag = os.environ.get("SYSTEM_VERSION", "local")
    dockerhub = os.environ.get("DOCKERHUB_REPO", "ccdlstaging")
    cmd.append(f"{dockerhub}/dr_{image_name}:{image_tag}")
    cmd.extend(["bash", "-c", coverage_command(extra)])
    return run(cmd)


def cmd_test_workers(argv):
    p = argparse.ArgumentParser(
        prog="rbio test:workers",
        description=(
            "Run the workers test suite. With -t TAG, builds + tests just one "
            "worker; without, iterates over every worker tag in order."
        ),
        epilog=(
            "tags: " + ", ".join(WORKER_TAGS) + "\n\n"
            "wraps:\n"
            "  (per tag) docker buildx bake <image>\n"
            "  (per tag) docker run --network refinebio_default --rm <image> bash -c '<coverage runner>'"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "-t",
        "--tag",
        choices=WORKER_TAGS,
        help="run tests for a single worker tag (default: every tag in order).",
    )
    args, forwarded = p.parse_known_args(argv)

    if (rc := require_drdb("test:workers")) != 0:
        return rc

    # Workers have their own test_volume next to the worker code (matches the
    # legacy `workers/run_tests.sh` mount). Committed fixtures + VCR cassettes
    # live under here; api/common/foreman use the top-level test_volume/.
    test_volume = REPO_ROOT / "workers" / "test_volume"
    _ensure_worker_test_data(args.tag, test_volume)

    # Open up bind-mount perms so the worker container (which runs as the
    # image's baked user, not the host user) can write into test_volume —
    # both the bind-mount root AND every host-created subdir (PCL/, QN/,
    # raw/TEST/..., etc.) need to be writable to whatever UID the container
    # runs as. Some test fixtures live in committed subdirs we can't chmod
    # (the chmod returns nonzero on those), so check=False and ignore
    # stderr; what matters is that fresh dirs we created above get opened up.
    # Mirrors the legacy script's `chmod -R a+rwX test_volume`.
    if not Globals.dry_run:
        subprocess.run(
            ["chmod", "-R", "a+rwX", str(test_volume)],
            check=False,
            stderr=subprocess.DEVNULL,
        )

    # Mirror the legacy script's mock-AWS-creds-for-smasher behavior: smasher
    # and compendia tests use VCR cassettes recorded against bogus credentials.
    if args.tag in (None, "smasher", "compendia"):
        os.environ["AWS_ACCESS_KEY_ID"] = "XXX"
        os.environ["AWS_SECRET_ACCESS_KEY"] = "XXX"

    tags = [args.tag] if args.tag else list(WORKER_TAGS)
    for tag in tags:
        image = WORKER_TAG_TO_IMAGE[tag]
        if (rc := bake_target(image)) != 0:
            return rc
        if (rc := _run_worker_image(image, tag, forwarded, test_volume)) != 0:
            return rc
    return 0


COMMANDS = [
    ("test:workers", cmd_test_workers, "run the workers test suite (per-tag with -t)"),
]
