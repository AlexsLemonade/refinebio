"""Docker-side helpers: bake target invocation + container pre-flights."""

import subprocess

from lib._runtime import Globals, run, stderr


def bake_target(target):
    """Build a single bake target into the local docker daemon."""
    return run(["docker", "buildx", "bake", "-f", "docker-bake.hcl", target, "--load"])


def is_container_running(name):
    """True if a container with this name is up."""
    result = subprocess.run(
        ["docker", "container", "inspect", name, "--format", "{{.State.Running}}"],
        capture_output=True,
        text=True,
    )
    return result.returncode == 0 and result.stdout.strip() == "true"


def require_container(cmd_name, container, friendly):
    """Pre-flight: error with a hint if the named container isn't running. Skipped under --dry-run."""
    if Globals.dry_run or is_container_running(container):
        return 0
    stderr(f"rbio {cmd_name}: {friendly} ({container}) is not running")
    stderr("start it with:  rbio dev:up")
    return 1


def require_drdb(cmd_name):
    return require_container(cmd_name, "drdb", "postgres")


def require_dres(cmd_name):
    return require_container(cmd_name, "dres", "elasticsearch")


def require_services(cmd_name, friendly_to_container):
    """Pre-flight for commands needing multiple services up. Reports all missing in one message."""
    if Globals.dry_run:
        return 0
    missing = [
        name
        for name, container in friendly_to_container.items()
        if not is_container_running(container)
    ]
    if not missing:
        return 0
    required = " + ".join(friendly_to_container.keys())
    missing_str = ", ".join(missing)
    stderr(f"rbio {cmd_name}: requires {required} (not running: {missing_str})")
    stderr("start them with:  rbio dev:up")
    return 1
