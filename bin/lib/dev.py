"""dev:* — local dev stack lifecycle (compose up/down + ad-hoc image runs)."""

import argparse
import hashlib
import json
import os
import subprocess

from lib._docker import bake_target, require_drdb
from lib._runtime import REPO_ROOT, Globals, run, stderr

DEV_UP_DEFAULT_SERVICES = ["api", "postgres", "elasticsearch"]
DEV_UP_OPTIONAL_SERVICES = ["foreman"]
DEV_UP_ALL_SERVICES = sorted(set(DEV_UP_DEFAULT_SERVICES) | set(DEV_UP_OPTIONAL_SERVICES))

# Maps the job_name field returned by foreman's get_job_to_be_run to the
# worker image that knows how to run it. Downloader jobs always use the
# downloaders image, so they're handled separately and not in this map.
PROCESSOR_IMAGE_FOR_JOB_NAME = {
    "AFFY_TO_PCL": "affymetrix",
    "SALMON": "salmon",
    "ILLUMINA_TO_PCL": "illumina",
    "TRANSCRIPTOME_INDEX_LONG": "transcriptome",
    "TRANSCRIPTOME_INDEX_SHORT": "transcriptome",
    "NO_OP": "no_op",
}


def cmd_dev_up(argv):
    p = argparse.ArgumentParser(
        prog="rbio dev:up",
        description="Start the local dev stack.",
        epilog=(
            "services:\n"
            f"  default:   {', '.join(DEV_UP_DEFAULT_SERVICES)}\n"
            f"  optional:  {', '.join(DEV_UP_OPTIONAL_SERVICES)}\n"
            "\n"
            "  workers are not part of dev:up — they're ephemeral job\n"
            "  containers, started on demand by other commands.\n"
            "\n"
            "wraps: docker compose up -d <services>"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "-s",
        "--service",
        dest="services",
        action="append",
        default=[],
        metavar="SERVICE",
        choices=DEV_UP_ALL_SERVICES,
        help="service to start (repeatable). default: the full default stack.",
    )
    p.add_argument(
        "--wait", action="store_true", help="block until services pass their healthchecks"
    )
    args = p.parse_args(argv)

    services = list(args.services) if args.services else list(DEV_UP_DEFAULT_SERVICES)

    if not Globals.dry_run:
        missing = check_missing_local_images(services)
        if missing:
            stderr("rbio dev:up: missing local image(s):")
            for img in sorted(missing):
                stderr(f"  {img}")
            stderr("")
            stderr("build them first:  rbio build")
            return 1

    cmd = ["docker", "compose", "up", "-d"]
    if args.wait:
        cmd.append("--wait")
    cmd.extend(services)
    return run(cmd)


def check_missing_local_images(services):
    """Return refinebio-built images that compose needs but aren't in the local daemon."""
    cfg = subprocess.run(
        ["docker", "compose", "config", "--format", "json"],
        capture_output=True,
        text=True,
        cwd=str(REPO_ROOT),
    )
    if cfg.returncode != 0:
        return []
    try:
        config = json.loads(cfg.stdout)
    except json.JSONDecodeError:
        return []
    needed = [config.get("services", {}).get(s, {}).get("image", "") for s in services]
    needed = [i for i in needed if "/dr_" in i]
    if not needed:
        return []
    ls = subprocess.run(
        ["docker", "image", "ls", "--format", "{{.Repository}}:{{.Tag}}"],
        capture_output=True,
        text=True,
    )
    local = set(ls.stdout.splitlines()) if ls.returncode == 0 else set()
    return [i for i in needed if i not in local]


def cmd_dev_down(argv):
    p = argparse.ArgumentParser(
        prog="rbio dev:down",
        description=(
            "Stop the local dev stack. Volumes are preserved by default — "
            "your postgres data is safe."
        ),
        epilog="wraps: docker compose down",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "-v",
        "--volumes",
        action="store_true",
        help="also remove volumes (DESTRUCTIVE — wipes the db)",
    )
    p.add_argument("--rmi", action="store_true", help="remove images built by compose")
    args = p.parse_args(argv)

    cmd = ["docker", "compose", "down"]
    if args.volumes:
        cmd.append("--volumes")
    if args.rmi:
        cmd.extend(["--rmi", "local"])
    return run(cmd)


def cmd_dev_manage(argv):
    p = argparse.ArgumentParser(
        prog="rbio dev:manage",
        description=(
            "Build a sub-project image and run a Django management command in it. "
            "The image is built directly (not via compose) so it can mount the "
            "api/volume dir and pick up host AWS credentials."
        ),
        epilog=(
            "examples:\n"
            "  rbio dev:manage -i api_local -s api search_index --rebuild -f\n"
            "  rbio dev:manage -i foreman survey_all --accession E-GEOD-1\n"
            "\n"
            "wraps: docker build -f <service>/dockerfiles/Dockerfile.<image>; docker run python3 manage.py ..."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("-i", "--image", required=True, help="image variant (e.g. api_local, foreman)")
    p.add_argument(
        "-s",
        "--service",
        default="foreman",
        help="sub-project directory holding the Dockerfile (default: foreman)",
    )
    p.add_argument(
        "manage_args",
        nargs=argparse.REMAINDER,
        help="arguments forwarded to manage.py",
    )
    args = p.parse_args(argv)

    dockerfile = REPO_ROOT / args.service / "dockerfiles" / f"Dockerfile.{args.image}"
    if not Globals.dry_run and not dockerfile.is_file():
        stderr(f"rbio dev:manage: missing dockerfile: {dockerfile}")
        return 1
    if (rc := require_drdb("dev:manage")) != 0:
        return rc

    volume_dir = _ensure_volume_dir()
    if (rc := _build_dr_shell_image(dockerfile)) != 0:
        return rc
    return run(_dr_shell_run_argv(args, volume_dir))


def _ensure_volume_dir():
    """The api/volume dir is mounted into dr_shell; create + chmod if missing."""
    volume_dir = REPO_ROOT / "api" / "volume"
    volume_dir.mkdir(parents=True, exist_ok=True)
    subprocess.call(["chmod", "-R", "a+rwX", str(volume_dir)])
    return volume_dir


def _build_dr_shell_image(dockerfile):
    """Build the dockerfile into a single throwaway-tagged image (dr_shell)."""
    system_version = os.environ.get("SYSTEM_VERSION") or _branch_hash()
    dockerhub_repo = os.environ.get("DOCKERHUB_REPO", "ccdlstaging")
    return run(
        [
            "docker",
            "build",
            "--build-arg",
            f"DOCKERHUB_REPO={dockerhub_repo}",
            "--build-arg",
            f"SYSTEM_VERSION={system_version}",
            "--file",
            str(dockerfile),
            "--platform",
            "linux/amd64",
            "--tag",
            "dr_shell",
            ".",
        ]
    )


def _dr_shell_run_argv(args, volume_dir):
    """Compose the `docker run` argv for the dr_shell container."""
    env_file = REPO_ROOT / args.service / "environments" / "local"

    argv = ["docker", "run"]
    # Join the compose network so `database` (alias for drdb) and
    # `elasticsearch` (dres) resolve via docker's embedded DNS. The network
    # is created lazily by compose; require_drdb above ensures it exists.
    argv += ["--network", "refinebio_default"]
    # `--env KEY` (no value) inherits from the host env without baking the
    # secret into the docker process listing.
    argv += ["--env", "AWS_ACCESS_KEY_ID"]
    argv += ["--env", "AWS_SECRET_ACCESS_KEY"]
    argv += ["--env-file", str(env_file)]
    argv += ["--platform", "linux/amd64"]
    argv += ["--tty"]
    argv += ["--volume", f"{volume_dir}:/home/user/data_store"]
    argv += ["--volume", "/tmp:/tmp"]
    if os.isatty(1):
        argv += ["--interactive"]
    argv += ["dr_shell", "python3", "manage.py", *args.manage_args]
    return argv


def _branch_hash():
    """SHA-1 of the current branch name; fallback for SYSTEM_VERSION. Matches scripts/common.sh::get_branch_hash."""
    result = subprocess.run(
        ["git", "rev-parse", "--abbrev-ref", "HEAD"],
        capture_output=True,
        text=True,
        cwd=str(REPO_ROOT),
    )
    branch = (result.stdout.strip() or "unknown") + "\n"
    return hashlib.sha1(branch.encode()).hexdigest()


def cmd_dev_shell(argv):
    p = argparse.ArgumentParser(
        prog="rbio dev:shell",
        description="Open an interactive Python shell in the foreman container.",
        epilog="wraps: docker buildx bake foreman; docker compose run --rm foreman python3 manage.py shell",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.parse_args(argv)
    if (rc := bake_target("foreman")) != 0:
        return rc
    return run(["docker", "compose", "run", "--rm", "foreman", "python3", "manage.py", "shell"])


def cmd_dev_pipeline(argv):
    p = argparse.ArgumentParser(
        prog="rbio dev:pipeline",
        description=(
            "Run the full ingest pipeline end-to-end for a single accession: "
            "surveyor -> downloader -> processor."
        ),
        epilog=(
            "wraps:\n"
            "  foreman/run_management_command.sh survey_all --accession <code>\n"
            "  (loop) foreman/run_management_command.sh get_job_to_be_run\n"
            "  (loop) workers/run_job.sh -i <image> run_downloader_job|run_processor_job ..."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "accession_code",
        help="experiment accession (or organism for transcriptome indices) to ingest",
    )
    args = p.parse_args(argv)
    if (rc := require_drdb("dev:pipeline")) != 0:
        return rc

    manage = str(REPO_ROOT / "foreman" / "run_management_command.sh")
    run_job = str(REPO_ROOT / "workers" / "run_job.sh")

    if (rc := run([manage, "survey_all", "--accession", args.accession_code])) != 0:
        return rc

    while True:
        job, rc = _next_job(manage)
        if rc is not None:
            return rc
        if job is None:
            return 0
        if (rc := _run_pipeline_job(run_job, job)) != 0:
            return rc


def _next_job(manage_sh):
    """Returns (job_dict, None) when a job is available; (None, 0) when done;
    (None, rc) on error. Foreman prints the job spec as the last non-blank line."""
    result = subprocess.run([manage_sh, "get_job_to_be_run"], capture_output=True, text=True)
    if result.returncode != 0:
        stderr("rbio dev:pipeline: failed to fetch next job from foreman")
        return None, result.returncode
    lines = [line for line in result.stdout.splitlines() if line.strip()]
    if not lines:
        return None, 0
    try:
        return json.loads(lines[-1]), None
    except json.JSONDecodeError:
        return None, 0


def _run_pipeline_job(run_job_sh, job):
    job_name = job["job_name"]
    job_id = job["job_id"]
    print(f"Running {job_name} Job with id {job_id}!")

    if job["job_type"] == "DownloaderJob":
        subcommand = "run_downloader_job"
        image_name = "downloaders"
    else:
        subcommand = "run_processor_job"
        image_name = PROCESSOR_IMAGE_FOR_JOB_NAME.get(job_name)
        if not image_name:
            stderr(f"rbio dev:pipeline: no image mapped for processor job '{job_name}'")
            return 1

    return run(
        [
            run_job_sh,
            "-i",
            image_name,
            subcommand,
            f"--job-name={job_name}",
            f"--job-id={job_id}",
        ]
    )


COMMANDS = [
    ("dev:up", cmd_dev_up, "start the dev stack"),
    ("dev:down", cmd_dev_down, "stop the dev stack"),
    ("dev:manage", cmd_dev_manage, "build an image + run a Django management command in it"),
    ("dev:shell", cmd_dev_shell, "open an interactive python shell in foreman"),
    ("dev:pipeline", cmd_dev_pipeline, "run the survey -> download -> process pipeline locally"),
]
