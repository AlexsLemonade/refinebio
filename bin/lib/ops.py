"""ops:* — production batch operations (terminate jobs, RAM-hour reports)."""

import argparse
import datetime
import os

from lib._runtime import REPO_ROOT, Globals, stderr


def _load_env_file(path):
    """Mutate os.environ from a KEY=VALUE file (skip blank lines + #-comments)."""
    with open(path) as f:
        for raw in f:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            key, _, value = line.partition("=")
            if key:
                os.environ[key] = value


def cmd_ops_kill_jobs(argv):
    p = argparse.ArgumentParser(
        prog="rbio ops:kill-jobs",
        description=(
            "Terminate every running / queued AWS Batch job in the production "
            "worker queues. Reads AWS_REGION + REFINEBIO_JOB_QUEUE_WORKERS_NAMES "
            "from the current environment (or --env-file)."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--env-file",
        default=str(REPO_ROOT / "infrastructure" / "prod_env"),
        help="KEY=VALUE file to load before running (default: infrastructure/prod_env)",
    )
    args = p.parse_args(argv)

    if args.env_file and os.path.exists(args.env_file):
        _load_env_file(args.env_file)

    try:
        region = os.environ["AWS_REGION"]
        queue_names = os.environ["REFINEBIO_JOB_QUEUE_WORKERS_NAMES"].split(",")
    except KeyError as e:
        stderr(f"rbio ops:kill-jobs: missing required env var {e}")
        return 1

    if Globals.dry_run:
        stderr(f"DRY RUN: would terminate jobs in {len(queue_names)} queue(s) in {region}")
        return 0

    import boto3

    batch = boto3.client("batch", region_name=region)
    statuses = ["SUBMITTED", "PENDING", "RUNNABLE", "STARTING", "RUNNING"]
    for queue in queue_names:
        for status in statuses:
            token = None
            while True:
                kwargs = {"jobQueue": queue, "jobStatus": status}
                if token:
                    kwargs["nextToken"] = token
                page = batch.list_jobs(**kwargs)
                for job in page["jobSummaryList"]:
                    print(f"Deleting job: {job['jobId']}")
                    batch.terminate_job(jobId=job["jobId"], reason="kill_all_jobs")
                token = page.get("nextToken")
                if not token:
                    break
    return 0


def _ram_hours_for_job(job):
    if not (job.start_time and job.end_time and job.ram_amount):
        return 0
    hours = (job.end_time - job.start_time).total_seconds() / 3600
    gb = job.ram_amount / 1024
    return hours * gb


def _ram_hours_for_sample(pyrefinebio, accession_code, utc_start):
    sample = pyrefinebio.Sample.get(accession_code)
    if not sample.original_files:
        print(f"Found no original files for sample {accession_code}")
        return [accession_code, 0, 0, 0]

    downloader_hours = 0
    processor_hours = 0
    seen_dl, seen_p = set(), set()
    for original_file_id in sample.original_files:
        of = pyrefinebio.OriginalFile.get(original_file_id)
        for j in of.downloader_jobs:
            if j.id not in seen_dl and j.start_time and j.start_time > utc_start:
                downloader_hours += _ram_hours_for_job(j)
                seen_dl.add(j.id)
        for j in of.processor_jobs:
            if j.id not in seen_p and j.start_time and j.start_time > utc_start:
                processor_hours += _ram_hours_for_job(j)
                seen_p.add(j.id)
    return [accession_code, downloader_hours, processor_hours, downloader_hours + processor_hours]


def cmd_ops_ram_hours(argv):
    p = argparse.ArgumentParser(
        prog="rbio ops:ram-hours",
        description=(
            "Calculate how many RAM-hours an experiment or sample consumed.\n"
            "A RAM-hour is 1GB of memory for 1 hour."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("accession_code", help="experiment accession (or sample with --sample)")
    p.add_argument(
        "--start-date",
        default="2021-06-23",
        type=lambda s: datetime.datetime.strptime(s, "%Y-%m-%d"),
        help="ignore jobs started before this date (YYYY-MM-DD). default: 2021-06-23",
    )
    p.add_argument(
        "--sample",
        action="store_true",
        help="treat accession_code as a sample (default: experiment)",
    )
    args = p.parse_args(argv)

    import pyrefinebio
    import pytz

    utc_start = pytz.utc.localize(args.start_date)

    print(
        ", ".join(
            ["accession_code", "downloader_ram_hours", "processor_ram_hours", "total_ram_hours"]
        )
    )
    if args.sample:
        row = _ram_hours_for_sample(pyrefinebio, args.accession_code, utc_start)
        print(", ".join(map(str, row)))
    else:
        experiment = pyrefinebio.Experiment.get(args.accession_code)
        for sample in experiment.samples:
            row = _ram_hours_for_sample(pyrefinebio, sample.accession_code, utc_start)
            print(", ".join(map(str, row)))
    return 0


COMMANDS = [
    ("ops:kill-jobs", cmd_ops_kill_jobs, "terminate all queued/running AWS Batch jobs"),
    ("ops:ram-hours", cmd_ops_ram_hours, "report RAM-hours consumed by an experiment/sample"),
]
