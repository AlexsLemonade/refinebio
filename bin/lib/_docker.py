"""Docker-side helpers: bake target invocation + container pre-flights."""

import os
import subprocess

from lib._runtime import Globals, run, stderr
from lib.image import WORKER_R_ENV_NAMES, resolve_base, resolve_worker_r_env


def resolve_bake_env():
    """Resolve every content-hash env var the bake file references.
    Returns (env_patch, ok). UBUNTU_DIGEST + BASE_HASH + each worker's
    <WORKER>_R_ENV_HASH are resolved unconditionally; bake picks which
    reach which target. Callers that already supplied a value via the
    environment take precedence over auto-resolution."""
    env = {}
    ubuntu_digest = os.environ.get("UBUNTU_DIGEST")
    base_hash, digest = resolve_base(ubuntu_digest=ubuntu_digest)
    if base_hash is None:
        return {}, False
    env["UBUNTU_DIGEST"] = digest
    env["BASE_HASH"] = base_hash
    for worker in WORKER_R_ENV_NAMES:
        var = f"{worker.upper()}_R_ENV_HASH"
        env[var] = os.environ.get(var) or resolve_worker_r_env(worker, base_hash=base_hash)
    return env, True


def bake_target(target):
    """Build a single bake target into the local docker daemon."""
    hash_env, ok = resolve_bake_env()
    if not ok:
        return 1
    env = {**os.environ, **hash_env}
    return run(
        ["docker", "buildx", "bake", "-f", "docker-bake.hcl", target, "--load"],
        env=env,
    )


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
