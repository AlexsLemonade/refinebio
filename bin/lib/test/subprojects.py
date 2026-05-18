"""test:api / test:common / test:foreman — Django subproject test runners."""

import argparse

from lib._docker import bake_target, require_services
from lib._runtime import REPO_ROOT, run
from lib.common import require_common_fresh
from lib.test._shared import coverage_command

# Per-subproject knobs for `rbio test:<name>`: which bake target supplies
# the test image, whether the suite needs elasticsearch reachable, and any
# default flags forwarded to manage.py test.
TEST_SUBPROJECTS = {
    "api": {
        "image": "api_local",
        "needs_es": True,
        "default_args": [],
    },
    "common": {
        "image": "common_tests",
        "needs_es": True,
        "default_args": ["--parallel"],
    },
    "foreman": {
        "image": "foreman",
        "needs_es": False,
        "default_args": ["--exclude-tag=manual"],
    },
}


def _test_subproject(name, argv):
    cfg = TEST_SUBPROJECTS[name]
    p = argparse.ArgumentParser(
        prog=f"rbio test:{name}",
        description=(
            f"Run the {name} test suite via docker compose. Any extra args "
            "are forwarded to manage.py test (e.g. -t TAG, --keepdb)."
        ),
        epilog=(
            "wraps:\n"
            f"  docker buildx bake {cfg['image']}\n"
            f"  docker compose run --rm test_{name} bash -c '<coverage runner>'"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--allow-stale-common",
        action="store_true",
        help="skip the pre-flight check that common/dist is up to date with common/ source",
    )
    args, forwarded = p.parse_known_args(argv)

    services = {"postgres": "drdb"}
    if cfg["needs_es"]:
        services["elasticsearch"] = "dres"
    if (rc := require_services(f"test:{name}", services)) != 0:
        return rc
    if not args.allow_stale_common and (rc := require_common_fresh(f"test:{name}")) != 0:
        return rc
    if (rc := bake_target(cfg["image"])) != 0:
        return rc

    # Pre-create the bind-mount target so dockerd doesn't create it as root
    # and prevent writing coverage.xml.
    (REPO_ROOT / "test_volume").mkdir(exist_ok=True)

    extra = cfg["default_args"] + forwarded
    return run(
        [
            "docker",
            "compose",
            "run",
            "--rm",
            f"test_{name}",
            "bash",
            "-c",
            coverage_command(extra),
        ]
    )


def cmd_test_api(argv):
    return _test_subproject("api", argv)


def cmd_test_common(argv):
    return _test_subproject("common", argv)


def cmd_test_foreman(argv):
    return _test_subproject("foreman", argv)


COMMANDS = [
    ("test:api", cmd_test_api, "run the api test suite"),
    ("test:common", cmd_test_common, "run the common test suite"),
    ("test:foreman", cmd_test_foreman, "run the foreman test suite"),
]
