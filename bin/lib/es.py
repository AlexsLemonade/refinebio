"""es:* — elasticsearch index + container ops."""

import argparse
import subprocess
import time

from lib._docker import bake_target, require_dres
from lib._runtime import Globals, run, stderr


def cmd_es_rebuild(argv):
    p = argparse.ArgumentParser(
        prog="rbio es:rebuild",
        description="Rebuild the elasticsearch index from scratch.",
        epilog=(
            "wraps:\n"
            "  docker buildx bake api_local\n"
            "  docker compose run --rm api python3 manage.py search_index --rebuild -f"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.parse_args(argv)
    if (rc := require_dres("es:rebuild")) != 0:
        return rc
    if (rc := bake_target("api_local")) != 0:
        return rc
    return run(
        [
            "docker",
            "compose",
            "run",
            "--rm",
            "api",
            "python3",
            "manage.py",
            "search_index",
            "--rebuild",
            "-f",
        ]
    )


def _wait_for_es(timeout=30):
    """Poll ES cluster health until reachable, up to `timeout` seconds."""
    deadline = time.time() + timeout
    while time.time() < deadline:
        result = subprocess.run(
            ["docker", "exec", "dres", "curl", "-sf", "http://localhost:9200/_cluster/health"],
            capture_output=True,
        )
        if result.returncode == 0:
            return 0
        time.sleep(1)
    stderr(f"rbio es:reset: elasticsearch did not become ready within {timeout}s")
    return 1


def cmd_es_reset(argv):
    p = argparse.ArgumentParser(
        prog="rbio es:reset",
        description=(
            "Tear down + recreate the elasticsearch container (clears any stale "
            "in-container state), then rebuild the search index."
        ),
        epilog=(
            "wraps:\n"
            "  docker compose rm -fsv elasticsearch\n"
            "  docker compose up -d elasticsearch\n"
            "  (poll until ready)\n"
            "  rbio es:rebuild"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.parse_args(argv)
    if (rc := run(["docker", "compose", "rm", "-fsv", "elasticsearch"])) != 0:
        return rc
    if (rc := run(["docker", "compose", "up", "-d", "elasticsearch"])) != 0:
        return rc
    if not Globals.dry_run and (rc := _wait_for_es()) != 0:
        return rc
    return cmd_es_rebuild([])


COMMANDS = [
    ("es:rebuild", cmd_es_rebuild, "rebuild the elasticsearch index"),
    ("es:reset", cmd_es_reset, "wipe + recreate the elasticsearch container then rebuild"),
]
