"""install:* — one-time host setup actions. Thin Python wrappers over install/*.sh."""

import argparse

from lib._requirements import epilog_from, requires
from lib._runtime import REPO_ROOT, run


@requires(tools=[])
def cmd_install_deps(argv):
    p = argparse.ArgumentParser(
        prog="rbio install:deps",
        description="Install missing OS-level dev dependencies (idempotent).",
        epilog=epilog_from(
            cmd_install_deps,
            "wraps: ./install/install_deps.sh\n"
            "\n"
            "Run `rbio debug:deps` first to see what's missing. The script\n"
            "is idempotent — already-present tools are skipped.",
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("-y", "--yes", action="store_true", help="skip confirmation prompts")
    args = p.parse_args(argv)
    cmd = ["./install/install_deps.sh"]
    if args.yes:
        cmd.append("--yes")
    return run(cmd, cwd=str(REPO_ROOT))


@requires(tools=["python3"])
def cmd_install_venv(argv):
    p = argparse.ArgumentParser(
        prog="rbio install:venv",
        description=(
            "Create dr_env/ with pip-tools + every subproject's requirements + "
            "the data_refinery_common package as editable. Used by IDEs for "
            "imports/autocomplete/jump-to-def, and by `rbio common:build-sdist` "
            "on PEP-668 distros where system python3 lacks setuptools."
        ),
        epilog=epilog_from(
            cmd_install_venv,
            "wraps: ./install/setup_venv.sh\n"
            "\n"
            "After this completes, point your IDE at dr_env/bin/python3.",
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.parse_args(argv)
    return run(["./install/setup_venv.sh"], cwd=str(REPO_ROOT))


COMMANDS = [
    ("install:deps", cmd_install_deps, "install missing OS-level dev dependencies"),
    ("install:venv", cmd_install_venv, "set up dr_env/ for IDE + pip-compile"),
]
