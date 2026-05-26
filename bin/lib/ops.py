"""ops:* — manual ops tools that act on AWS infrastructure."""

import argparse
import json
import os
import re
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime

from lib._requirements import epilog_from, requires
from lib._runtime import stderr

_STATUSES = ["SUBMITTED", "PENDING", "RUNNABLE", "STARTING", "RUNNING"]


@requires(
    env=["DEPLOY_USER"],
    env_optional={"DEPLOY_ENV": "dev"},
    tools=["aws"],
)
def cmd_ops_kill_jobs(argv):
    p = argparse.ArgumentParser(
        prog="rbio ops:kill-jobs",
        description="Terminate Batch jobs for a user/env scope. Interactive by default.",
        epilog=epilog_from(
            cmd_ops_kill_jobs,
            "Queue scope = {DEPLOY_USER}-{DEPLOY_ENV}\n"
            "Examples of scope:  circleci-prod, circleci-staging, david-dev\n"
            "\n"
            "AWS auth from standard env (AWS_PROFILE or static creds).\n"
            "AWS region from AWS_REGION. Verify with `rbio debug:env`.\n"
            "\n"
            "By default: prompts to pick queue(s), then confirms before terminating.\n"
            "Pass --all + --yes for unattended use.\n"
            "\n"
            "wraps: aws batch describe-job-queues | list-jobs | terminate-job",
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--all",
        action="store_true",
        help="skip the picker; target every queue in scope",
    )
    p.add_argument(
        "--yes",
        "-y",
        action="store_true",
        help="skip the destructive confirmation prompt",
    )
    args = p.parse_args(argv)

    user = os.environ.get("DEPLOY_USER")
    if not user:
        stderr("ops:kill-jobs: DEPLOY_USER required (run `rbio debug:env` to check)")
        return 1
    env = os.environ.get("DEPLOY_ENV", "dev")
    scope = f"{user}-{env}"

    queues = _discover_queues(scope)
    if not queues:
        stderr(f"ops:kill-jobs: no queues found for scope '{scope}'")
        return 1

    check_time = datetime.now().astimezone().strftime("%H:%M:%S %Z")
    jobs_by_queue = _collect_jobs(queues)
    counts = {q: len(ids) for q, ids in jobs_by_queue.items()}

    selected = queues if args.all else _select_queues(queues, counts, check_time)
    if not selected:
        print("No queues selected. Exiting.")
        return 0

    selected_jobs = {q: jobs_by_queue[q] for q in selected}
    if not args.yes and not _confirm_kill(selected_jobs, counts, check_time):
        print("Aborted.")
        return 0

    for queue, ids in selected_jobs.items():
        for job_id in ids:
            print(f"Terminating: {job_id}")
            subprocess.run(
                [
                    "aws",
                    "batch",
                    "terminate-job",
                    "--job-id",
                    job_id,
                    "--reason",
                    "rbio ops:kill-jobs",
                ],
                check=True,
            )
    return 0


def _discover_queues(scope):
    """Return queue names whose terraform-managed name suffix matches the scope."""
    result = subprocess.run(
        ["aws", "batch", "describe-job-queues", "--output", "json"],
        capture_output=True,
        text=True,
        check=True,
    )
    data = json.loads(result.stdout)
    # Queue naming follows infrastructure/batch/job_queue.tf:
    #   data-refinery-batch-<role>-queue-<user>-<stage>[-<index>]
    pattern = re.compile(rf"-queue-{re.escape(scope)}(-\d+)?$")
    return sorted(q["jobQueueName"] for q in data["jobQueues"] if pattern.search(q["jobQueueName"]))


def _collect_jobs(queues):
    """List jobs in parallel across (queue × in-flight status), return {queue: [id, ...]}."""
    print(f"Listing jobs across {len(queues)} queue(s)...")
    by_queue = {q: [] for q in queues}
    tasks = [(q, s) for q in queues for s in _STATUSES]
    with ThreadPoolExecutor(max_workers=min(20, len(tasks))) as ex:
        for queue, ids in ex.map(lambda t: (t[0], _list_job_ids(*t)), tasks):
            by_queue[queue].extend(ids)
    return by_queue


def _list_job_ids(queue, status):
    ids = []
    next_token = None
    while True:
        cmd = [
            "aws",
            "batch",
            "list-jobs",
            "--job-queue",
            queue,
            "--job-status",
            status,
            "--output",
            "json",
        ]
        if next_token:
            cmd.extend(["--starting-token", next_token])
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        data = json.loads(result.stdout)
        ids.extend(j["jobId"] for j in data.get("jobSummaryList", []))
        next_token = data.get("nextToken")
        if not next_token:
            break
    return ids


def _select_queues(queues, counts, check_time):
    if not sys.stdin.isatty():
        stderr("ops:kill-jobs: not a TTY; pass --all to skip picker")
        return []
    max_name = max(len(q) for q in queues)
    max_count_width = max(len(str(c)) for c in counts.values())
    print(f"\nQueues in scope ({len(queues)})  [estimated counts from {check_time}]:")
    for i, q in enumerate(queues, 1):
        c = counts[q]
        print(f"  {i}. {q:<{max_name}}    [{c:>{max_count_width}} jobs]")
    sel = input("\nSelect queue(s) [#, comma-separated, 'all', or empty to cancel]: ").strip()
    if not sel:
        return []
    if sel.lower() in ("all", "a"):
        return queues
    try:
        indices = [int(x.strip()) - 1 for x in sel.split(",")]
        return [queues[i] for i in indices if 0 <= i < len(queues)]
    except (ValueError, IndexError):
        stderr(f"Invalid selection: {sel!r}")
        return []


def _confirm_kill(selected_jobs, counts, check_time):
    if not sys.stdin.isatty():
        stderr("ops:kill-jobs: not a TTY; pass --yes to skip confirmation")
        return False
    selected = list(selected_jobs.keys())
    max_name = max(len(q) for q in selected)
    max_count_width = max(len(str(counts[q])) for q in selected)
    print(
        f"\nAbout to terminate jobs in {len(selected)} queue(s)  "
        f"[estimated counts from {check_time}]:"
    )
    for q in selected:
        c = counts[q]
        print(f"  {q:<{max_name}}    [{c:>{max_count_width}} jobs]")
    response = input("\nContinue? [y/N]: ").strip().lower()
    return response in ("y", "yes")


COMMANDS = [
    ("ops:kill-jobs", cmd_ops_kill_jobs, "terminate Batch jobs for a user/env scope"),
]
