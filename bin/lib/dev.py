"""dev:* — local dev stack lifecycle (compose up/down)."""

import argparse
import json
import subprocess

from lib._runtime import REPO_ROOT, Globals, run, stderr

DEV_UP_DEFAULT_SERVICES = ["api", "postgres", "elasticsearch"]
DEV_UP_OPTIONAL_SERVICES = ["foreman"]
DEV_UP_ALL_SERVICES = sorted(set(DEV_UP_DEFAULT_SERVICES) | set(DEV_UP_OPTIONAL_SERVICES))


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


COMMANDS = [
    ("dev:up", cmd_dev_up, "start the dev stack"),
    ("dev:down", cmd_dev_down, "stop the dev stack"),
]
