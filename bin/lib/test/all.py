"""test:all — clean-state entry point: common:update + every test:<subproject> suite."""

import argparse

from lib import common
from lib._docker import require_drdb
from lib.test.subprojects import cmd_test_api, cmd_test_common, cmd_test_foreman
from lib.test.workers import cmd_test_workers


def cmd_test_all(argv):
    p = argparse.ArgumentParser(
        prog="rbio test:all",
        description=(
            "Run the full test suite from a clean state: common:update "
            "(rebuild migrations + sdist), then test:api, test:common, "
            "test:foreman, test:workers. Aborts on the first failure. Each "
            "suite is invoked with --allow-stale-common because we just "
            "refreshed."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.parse_args(argv)

    if (rc := require_drdb("test:all")) != 0:
        return rc
    if (rc := common.cmd_common_update([])) != 0:
        return rc

    for cmd_fn in (cmd_test_api, cmd_test_common, cmd_test_foreman, cmd_test_workers):
        if (rc := cmd_fn(["--allow-stale-common"])) != 0:
            return rc
    return 0


COMMANDS = [
    ("test:all", cmd_test_all, "full chain: common:update + every test:* suite"),
]
