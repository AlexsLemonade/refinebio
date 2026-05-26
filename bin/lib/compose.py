"""compose:* — wrappers around docker compose for the local dev stack."""

import argparse
import json
import subprocess

from lib._requirements import epilog_from, requires
from lib._runtime import REPO_ROOT, Globals, run, stderr

COMPOSE_UP_DEFAULT_SERVICES = ["api", "postgres", "elasticsearch"]
COMPOSE_UP_OPTIONAL_SERVICES = ["foreman"]
COMPOSE_UP_ALL_SERVICES = sorted(
    set(COMPOSE_UP_DEFAULT_SERVICES) | set(COMPOSE_UP_OPTIONAL_SERVICES)
)


@requires(tools=["docker"])
def cmd_compose_up(argv):
    p = argparse.ArgumentParser(
        prog="rbio compose:up",
        description="Start the local dev stack.",
        epilog=(
            "services:\n"
            f"  default:   {', '.join(COMPOSE_UP_DEFAULT_SERVICES)}\n"
            f"  optional:  {', '.join(COMPOSE_UP_OPTIONAL_SERVICES)}\n"
            "\n"
            "  workers are not part of compose:up — they're ephemeral job\n"
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
        choices=COMPOSE_UP_ALL_SERVICES,
        help="service to start (repeatable). default: the full default stack.",
    )
    p.add_argument(
        "--wait", action="store_true", help="block until services pass their healthchecks"
    )
    args = p.parse_args(argv)

    services = list(args.services) if args.services else list(COMPOSE_UP_DEFAULT_SERVICES)

    if not Globals.dry_run:
        missing = check_missing_local_images(services)
        if missing:
            stderr("rbio compose:up: missing local image(s):")
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


@requires(tools=["docker"])
def cmd_compose_down(argv):
    p = argparse.ArgumentParser(
        prog="rbio compose:down",
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


@requires(tools=["docker"])
def cmd_compose_manage(argv):
    p = argparse.ArgumentParser(
        prog="rbio compose:manage",
        description="Run a Django management command in a service container.",
        epilog=epilog_from(
            cmd_compose_manage,
            "examples:\n"
            "  rbio compose:manage foreman shell\n"
            "  rbio compose:manage foreman clear_database\n"
            "  rbio compose:manage api search_index --rebuild -f\n"
            "\n"
            "shell is just another management command (python3 manage.py shell),\n"
            "so it works through the same path as everything else.\n"
            "\n"
            "wraps: docker compose run --rm <service> python3 manage.py <args>",
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("service", help="compose service (foreman, api, workers, ...)")
    p.add_argument(
        "args",
        nargs=argparse.REMAINDER,
        help="forwarded to python3 manage.py",
    )
    args = p.parse_args(argv)
    return run(
        [
            "docker",
            "compose",
            "run",
            "--rm",
            args.service,
            "python3",
            "manage.py",
            *args.args,
        ]
    )


COMMANDS = [
    ("compose:up", cmd_compose_up, "start the dev stack"),
    ("compose:down", cmd_compose_down, "stop the dev stack"),
    ("compose:manage", cmd_compose_manage, "run a management command in a service container"),
]
