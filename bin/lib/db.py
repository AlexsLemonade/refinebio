"""db:* — state operations against the running drdb container."""

import argparse
import subprocess

from lib import common
from lib._docker import require_drdb
from lib._runtime import Globals, run, stderr


def _db_exec(*sql_args, db=None):
    """`docker exec drdb psql -U postgres [-d <db>] <args>`."""
    cmd = ["docker", "exec", "drdb", "psql", "-U", "postgres"]
    if db:
        cmd.extend(["-d", db])
    cmd.extend(sql_args)
    return run(cmd)


def _db_exists(name):
    """True if a postgres database named `name` exists in drdb."""
    result = subprocess.run(
        [
            "docker",
            "exec",
            "drdb",
            "psql",
            "-U",
            "postgres",
            "-tAc",
            f"SELECT 1 FROM pg_database WHERE datname='{name}'",
        ],
        capture_output=True,
        text=True,
    )
    return result.returncode == 0 and result.stdout.strip() == "1"


def cmd_db_init(argv):
    p = argparse.ArgumentParser(
        prog="rbio db:init",
        description="Bootstrap the data_refinery database, role, and hstore extension.",
        epilog="wraps: docker exec drdb psql -U postgres ...",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.parse_args(argv)
    if (rc := require_drdb("db:init")) != 0:
        return rc
    if not Globals.dry_run and _db_exists("data_refinery"):
        stderr("rbio db:init: data_refinery database already exists.")
        stderr("to drop and recreate run:  rbio db:reset")
        return 1
    statements = [
        "CREATE DATABASE data_refinery",
        "CREATE ROLE data_refinery_user WITH LOGIN PASSWORD 'data_refinery_password'",
        "GRANT ALL PRIVILEGES ON DATABASE data_refinery TO data_refinery_user",
        "ALTER USER data_refinery_user CREATEDB",
        "ALTER ROLE data_refinery_user superuser",
    ]
    for stmt in statements:
        if (rc := _db_exec("-c", stmt)) != 0:
            return rc
    return _db_exec("-c", "CREATE EXTENSION IF NOT EXISTS hstore", db="data_refinery")


def cmd_db_psql(argv):
    p = argparse.ArgumentParser(
        prog="rbio db:psql",
        description="Open an interactive psql shell against the data_refinery database.",
        epilog="wraps: docker exec -it drdb psql -U postgres -d data_refinery",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.parse_args(argv)
    if (rc := require_drdb("db:psql")) != 0:
        return rc
    return run(["docker", "exec", "-it", "drdb", "psql", "-U", "postgres", "-d", "data_refinery"])


def cmd_db_reset(argv):
    p = argparse.ArgumentParser(
        prog="rbio db:reset",
        description="Drop and re-create the data_refinery database, then re-run migrations.",
        epilog=(
            "wraps:\n"
            "  docker exec drdb psql -U postgres -c 'DROP DATABASE/ROLE ...'\n"
            "  rbio db:init\n"
            "  rbio common:apply"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.parse_args(argv)
    if (rc := require_drdb("db:reset")) != 0:
        return rc
    # Drop in reverse order of init; ignore failures (idempotent reset).
    _db_exec("-c", "DROP DATABASE IF EXISTS data_refinery")
    _db_exec("-c", "DROP ROLE IF EXISTS data_refinery_user")
    if (rc := cmd_db_init([])) != 0:
        return rc
    return common.cmd_common_apply([])


COMMANDS = [
    ("db:init", cmd_db_init, "bootstrap the data_refinery database + role"),
    ("db:psql", cmd_db_psql, "open a psql shell against data_refinery"),
    ("db:reset", cmd_db_reset, "drop + recreate data_refinery, then re-apply migrations"),
]
