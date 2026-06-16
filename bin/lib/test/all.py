"""test:all — orchestrate every project test suite sequentially."""

import argparse

from lib._docker import require_services
from lib._requirements import epilog_from, requires
from lib._runtime import REPO_ROOT
from lib.common import cmd_common_update
from lib.test.subprojects import cmd_test_api, cmd_test_common, cmd_test_foreman
from lib.test.workers import cmd_test_workers


@requires(tools=["docker"])
def cmd_test_all(argv):
    p = argparse.ArgumentParser(
        prog="rbio test:all",
        description="Run every project test suite (api, common, foreman, workers).",
        epilog=epilog_from(
            cmd_test_all,
            "wraps:\n"
            "  rbio common:update           (apply migrations, rebuild common sdist)\n"
            "  rbio test:api\n"
            "  rbio test:common\n"
            "  rbio test:foreman\n"
            "  rbio test:workers\n"
            "\n"
            "Suites run sequentially against the same local postgres + ES.\n"
            "First failure short-circuits the chain.",
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.parse_args(argv)

    if (rc := require_services("test:all", {"postgres": "drdb", "elasticsearch": "dres"})) != 0:
        return rc

    test_volume = REPO_ROOT / "test_volume"
    test_volume.mkdir(exist_ok=True)
    test_volume.chmod(0o777)

    for fn in [
        cmd_common_update,
        cmd_test_api,
        cmd_test_common,
        cmd_test_foreman,
        cmd_test_workers,
    ]:
        if (rc := fn([])) != 0:
            return rc
    return 0


COMMANDS = [
    ("test:all", cmd_test_all, "run every test suite sequentially"),
]
