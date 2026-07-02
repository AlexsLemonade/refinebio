"""common:* — data_refinery_common schema + sdist workflow."""

import argparse
import os
from pathlib import Path

from lib._docker import bake_target, require_drdb
from lib._requirements import epilog_from, requires
from lib._runtime import REPO_ROOT, Globals, run


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


def _sdist_python():
    # `setup.py sdist` needs setuptools. System python3 on PEP-668 distros
    # (Ubuntu 23.04+, Debian 12+) won't accept casual pip-installs into the
    # system site, so we prefer the venv created by `rbio install:venv`.
    # Falls back to `python3` when the venv isn't present.
    venv_python = REPO_ROOT / "dr_env" / "bin" / "python3"
    if venv_python.exists():
        return str(venv_python)
    return "python3"


def _sdist_needs_rebuild(common_dir, sdist_file):
    sdist_mtime = sdist_file.stat().st_mtime
    sources = list((common_dir / "data_refinery_common").rglob("*.py"))
    sources += [common_dir / "MANIFEST.in", common_dir / "setup.py"]
    return any(f.stat().st_mtime > sdist_mtime for f in sources if f.exists())


def rebuild_common_sdist():
    """Rebuild the common sdist tarball if any source file is newer than the existing one.

    Used by `rbio common:build-sdist` (standalone), `rbio common:update`
    (as the final step), and `rbio build` (as a pre-bake step so downstream
    images pick up an updated package).

    The mtime gate is the cache-correctness primitive: when source is
    unchanged, we skip `setup.py sdist` entirely and reuse the existing
    tarball. Docker's COPY layer hash holds → downstream `pip install
    common` layers stay cached. When source changes, a fresh sdist is
    produced and the cache busts correctly.
    """
    common_dir = REPO_ROOT / "common"
    dist_dir = common_dir / "dist"
    existing = list(dist_dir.glob("data*refinery*common-*.tar.gz"))

    if existing and not _sdist_needs_rebuild(common_dir, existing[0]):
        return 0

    if not Globals.dry_run:
        # setup.py reads the version from this file
        (common_dir / "version").write_text(os.environ.get("SYSTEM_VERSION", "local") + "\n")
        for old in existing:
            old.unlink()

    return run([_sdist_python(), "setup.py", "sdist"], cwd=common_dir)


@requires(tools=["docker"])
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


@requires(tools=["docker"])
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


@requires(tools=["docker"])
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
            "  rbio common:build-sdist"
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
    return rebuild_common_sdist()


@requires(tools=["python3"])
def cmd_common_build_sdist(argv):
    p = argparse.ArgumentParser(
        prog="rbio common:build-sdist",
        description="Rebuild the data_refinery_common sdist if any source file is newer than the existing tarball.",
        epilog=epilog_from(
            cmd_common_build_sdist,
            "wraps:\n"
            "  cd common && python3 setup.py sdist\n"
            "\n"
            "Skips the rebuild entirely when no source file is newer than the\n"
            "current sdist — keeps docker layer cache valid for downstream images.\n"
            "Uses dr_env/bin/python3 when the venv exists (system python3 lacks\n"
            "setuptools on PEP-668 distros).",
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.parse_args(argv)
    return rebuild_common_sdist()


COMMANDS = [
    ("common:make-migrations", cmd_common_make_migrations, "generate Django migration files"),
    ("common:apply", cmd_common_apply, "apply migrations to the database"),
    ("common:update", cmd_common_update, "full chain: make + apply + rebuild common sdist"),
    (
        "common:build-sdist",
        cmd_common_build_sdist,
        "rebuild data_refinery_common sdist (mtime-gated)",
    ),
]
