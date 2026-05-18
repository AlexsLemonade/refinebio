"""common:* — data_refinery_common schema + sdist workflow."""

import argparse
import os

from lib._docker import bake_target, require_drdb
from lib._runtime import REPO_ROOT, Globals, run, stderr


def is_common_stale():
    """True if any common/data_refinery_common/*.py is newer than the most recent common/dist sdist."""
    sdist_dir = REPO_ROOT / "common" / "dist"
    sdists = list(sdist_dir.glob("*.tar.gz"))
    if not sdists:
        return True
    newest_sdist = max(s.stat().st_mtime for s in sdists)
    src = REPO_ROOT / "common" / "data_refinery_common"
    if not src.is_dir():
        return False
    newest_src = max((p.stat().st_mtime for p in src.rglob("*.py")), default=0)
    return newest_src > newest_sdist


def require_common_fresh(cmd_name):
    """Pre-flight: fail with a hint if common/ source is newer than the last sdist build.

    Caller is responsible for honoring a --allow-stale-common flag; this helper
    is unconditional.
    """
    if not is_common_stale():
        return 0
    stderr(f"rbio {cmd_name}: common/ source has changed since the last sdist build")
    stderr("  fix:    rbio common:update")
    stderr(f"  bypass: rbio {cmd_name} --allow-stale-common")
    return 1


def _compose_run_migrations(*manage_args):
    """`docker compose run --rm migrations python3 manage.py <args>`."""
    return run(
        [
            "docker",
            "compose",
            "run",
            "--rm",
            "migrations",
            "python3",
            "manage.py",
            *manage_args,
        ]
    )


def cmd_common_make_migrations(argv):
    p = argparse.ArgumentParser(
        prog="rbio common:make-migrations",
        description="Generate Django migration files for data_refinery_common.",
        epilog=(
            "wraps:\n"
            "  docker buildx bake migrations  (build the migrations image)\n"
            "  docker compose run --rm migrations python3 manage.py makemigrations data_refinery_common"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.parse_args(argv)
    if (rc := require_drdb("common:make-migrations")) != 0:
        return rc
    if (rc := bake_target("migrations")) != 0:
        return rc
    return _compose_run_migrations("makemigrations", "data_refinery_common")


def cmd_common_apply(argv):
    p = argparse.ArgumentParser(
        prog="rbio common:apply",
        description="Apply Django migrations + create the cache table.",
        epilog=(
            "wraps:\n"
            "  docker buildx bake migrations\n"
            "  docker compose run --rm migrations python3 manage.py migrate\n"
            "  docker compose run --rm migrations python3 manage.py createcachetable"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.parse_args(argv)
    if (rc := require_drdb("common:apply")) != 0:
        return rc
    if (rc := bake_target("migrations")) != 0:
        return rc
    if (rc := _compose_run_migrations("migrate")) != 0:
        return rc
    return _compose_run_migrations("createcachetable")


def cmd_common_update(argv):
    p = argparse.ArgumentParser(
        prog="rbio common:update",
        description=(
            "Full after-models-changed workflow: regenerate migrations, apply "
            "them to the database, and rebuild the common sdist so downstream "
            "images pick up the updated package."
        ),
        epilog=(
            "wraps:\n"
            "  rbio common:make-migrations\n"
            "  rbio common:apply\n"
            "  cd common && python3 setup.py sdist  (host)"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.parse_args(argv)
    if (rc := require_drdb("common:update")) != 0:
        return rc
    if (rc := bake_target("migrations")) != 0:
        return rc
    if (rc := _compose_run_migrations("makemigrations", "data_refinery_common")) != 0:
        return rc
    if (rc := _compose_run_migrations("migrate")) != 0:
        return rc
    if (rc := _compose_run_migrations("createcachetable")) != 0:
        return rc
    # Rebuild common sdist on the host so other images (built by bake) pick it
    # up via Dockerfile COPY common/dist/...
    common_dir = REPO_ROOT / "common"
    if not Globals.dry_run:
        # setup.py reads the version from this file
        (common_dir / "version").write_text(os.environ.get("SYSTEM_VERSION", "local") + "\n")
        for old in (common_dir / "dist").glob("*"):
            old.unlink()
    return run(["python3", "setup.py", "sdist"], cwd=common_dir)


COMMANDS = [
    ("common:make-migrations", cmd_common_make_migrations, "generate Django migration files"),
    ("common:apply", cmd_common_apply, "apply migrations to the database"),
    ("common:update", cmd_common_update, "full chain: make + apply + rebuild common sdist"),
]
