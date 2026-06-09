"""dev:* — higher-level local dev workflows that orchestrate jobs."""

import argparse
import json
import subprocess

from lib._requirements import epilog_from, requires
from lib._runtime import Globals, run, stderr

WORKER_IMAGES = [
    "affymetrix",
    "agilent",
    "compendia",
    "downloaders",
    "illumina",
    "no_op",
    "salmon",
    "smasher",
    "transcriptome",
]

JOB_SUBCOMMANDS = ["run_downloader_job", "run_processor_job"]


# Maps processor job names to the worker image that handles them. Keep in sync
# with workers/data_refinery_workers/processors/management/commands/run_processor_job.py.
_PROCESSOR_JOB_IMAGE = {
    "AFFY_TO_PCL": "affymetrix",
    "AGILENT_TWOCOLOR_TO_PCL": "affymetrix",
    "SALMON": "salmon",
    "ILLUMINA_TO_PCL": "illumina",
    "TRANSCRIPTOME_INDEX_LONG": "transcriptome",
    "TRANSCRIPTOME_INDEX_SHORT": "transcriptome",
    "NO_OP": "no_op",
}


@requires(tools=["docker"])
def cmd_dev_pipeline(argv):
    p = argparse.ArgumentParser(
        prog="rbio dev:pipeline",
        description=(
            "Survey an accession locally, then drain the resulting downloader "
            "and processor jobs one at a time against the local stack."
        ),
        epilog=epilog_from(
            cmd_dev_pipeline,
            "wraps:\n"
            "  docker compose run --rm foreman python3 manage.py survey_all --accession ACC\n"
            "  (loop) docker compose run --rm foreman python3 manage.py get_job_to_be_run\n"
            "  (loop) docker compose run --rm <worker> python3 manage.py run_{processor,downloader}_job ...",
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("accession_code", help="experiment accession to survey")
    args = p.parse_args(argv)

    if (
        rc := run(
            [
                "docker",
                "compose",
                "run",
                "--rm",
                "foreman",
                "python3",
                "manage.py",
                "survey_all",
                "--accession",
                args.accession_code,
            ]
        )
    ) != 0:
        return rc

    while True:
        job = _get_job_to_run()
        if not job:
            return 0
        if (rc := _run_job(job)) != 0:
            return rc


def _get_job_to_run():
    if Globals.dry_run:
        return None
    # The management command prints the JSON as its last non-empty stdout line.
    result = subprocess.run(
        [
            "docker",
            "compose",
            "run",
            "--rm",
            "foreman",
            "python3",
            "manage.py",
            "get_job_to_be_run",
        ],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        return None
    lines = [ln for ln in result.stdout.splitlines() if ln.strip()]
    if not lines:
        return None
    try:
        return json.loads(lines[-1])
    except json.JSONDecodeError:
        return None


def _run_job(job):
    job_name = job["job_name"]
    job_id = job["job_id"]
    print(f"Running {job_name} Job with id {job_id}!")

    if job["job_type"] == "DownloaderJob":
        subcommand = "run_downloader_job"
        image = "downloaders"
    else:
        subcommand = "run_processor_job"
        image = _PROCESSOR_JOB_IMAGE.get(job_name, "downloaders")

    return run(
        [
            "docker",
            "compose",
            "run",
            "--rm",
            image,
            "python3",
            "manage.py",
            subcommand,
            f"--job-name={job_name}",
            f"--job-id={job_id}",
        ]
    )


@requires(tools=["docker"])
def cmd_dev_job(argv):
    p = argparse.ArgumentParser(
        prog="rbio dev:job",
        description="Run a single downloader/processor job in a worker container.",
        epilog=epilog_from(
            cmd_dev_job,
            "examples:\n"
            "  rbio dev:job downloaders run_downloader_job --job-name=SRA --job-id=12345\n"
            "  rbio dev:job affymetrix run_processor_job --job-name=AFFY_TO_PCL --job-id=54321\n"
            "\n"
            "agilent uses the affymetrix image (same binary).\n"
            "\n"
            "wraps: docker compose run --rm <image> python3 manage.py <subcommand> <args>",
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "image",
        choices=WORKER_IMAGES,
        help="worker image to run the job in",
    )
    p.add_argument(
        "subcommand",
        choices=JOB_SUBCOMMANDS,
        help="management command to run",
    )
    p.add_argument(
        "manage_args",
        nargs=argparse.REMAINDER,
        help="forwarded to manage.py (e.g. --job-name=SRA --job-id=12345)",
    )
    args = p.parse_args(argv)

    # agilent and affymetrix share the affymetrix image.
    image = "affymetrix" if args.image == "agilent" else args.image

    return run(
        [
            "docker",
            "compose",
            "run",
            "--rm",
            image,
            "python3",
            "manage.py",
            args.subcommand,
            *args.manage_args,
        ]
    )


@requires(tools=["docker"])
def cmd_dev_janitor(argv):
    p = argparse.ArgumentParser(
        prog="rbio dev:janitor",
        description="Run the janitor (cleanup) management command in a smasher container.",
        epilog=epilog_from(
            cmd_dev_janitor,
            "wraps: docker compose run --rm smasher python3 manage.py run_janitor",
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.parse_args(argv)
    return run(
        [
            "docker",
            "compose",
            "run",
            "--rm",
            "smasher",
            "python3",
            "manage.py",
            "run_janitor",
        ]
    )


COMMANDS = [
    ("dev:pipeline", cmd_dev_pipeline, "survey an accession + drain its job queue locally"),
    ("dev:job", cmd_dev_job, "run a single downloader/processor job"),
    ("dev:janitor", cmd_dev_janitor, "run the janitor cleanup command in a smasher container"),
]
