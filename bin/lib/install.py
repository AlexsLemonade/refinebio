"""install:* — one-time host setup actions. Thin Python wrappers over install/*.sh."""

import argparse

from lib._requirements import epilog_from, requires
from lib._runtime import REPO_ROOT, Globals, run, stderr

# Pip-compile runs inside this image so the lock resolves against the same
# Python the service images use. Avoids host-python mismatches like a Mac
# Homebrew python3.14 that can't build older sdists in the deps.
PIPTOOLS_IMAGE = "python:3.10-slim"

# Every subproject that has pip-tools-managed requirements. Kept here so
# install:compile has a single source of truth; add new subprojects by
# appending their .in path.
REQUIREMENTS_INS = [
    "api/requirements.in",
    "common/requirements.in",
    "foreman/requirements.in",
    "workers/data_refinery_workers/processors/requirements.in",
    "workers/data_refinery_workers/downloaders/requirements.in",
]


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


@requires(tools=["docker"])
def cmd_install_compile(argv):
    p = argparse.ArgumentParser(
        prog="rbio install:compile",
        description=(
            "Regenerate every subproject's requirements.txt from its "
            "requirements.in via pip-compile. Runs inside the "
            f"{PIPTOOLS_IMAGE} container so the lock resolves against the "
            "same Python the service images use — independent of the host "
            "Python, and independent of whether dr_env is working."
        ),
        epilog=epilog_from(
            cmd_install_compile,
            "subprojects compiled (one pip-compile each):\n"
            + "".join(f"  {p}\n" for p in REQUIREMENTS_INS)
            + "\n"
            "No flags: pip-compile only resolves new/changed .in lines; every\n"
            "other pin stays put. Use --package NAME to upgrade a single\n"
            "package (repeatable). `--dry-run` prints what would run.\n"
            "\n"
            "Any extra args after `--` are forwarded to pip-compile. Bare\n"
            "`-- --upgrade` upgrades every transitive dep and requires --yes.",
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "subproject",
        nargs="*",
        help=(
            "Subproject .in paths to compile (default: all). Accepts the "
            "directory (`api`) or the file (`api/requirements.in`)."
        ),
    )
    p.add_argument(
        "-p",
        "--package",
        action="append",
        default=[],
        metavar="NAME",
        help="Upgrade one package (repeatable); maps to --upgrade-package.",
    )
    p.add_argument(
        "-y",
        "--yes",
        action="store_true",
        help="Acknowledge a blanket `-- --upgrade` (otherwise refused).",
    )
    # Split on `--` before argparse runs: argparse would otherwise consume
    # the separator and treat the pip-compile flags as positional args.
    if "--" in argv:
        sep = argv.index("--")
        own_args, passthrough = argv[:sep], argv[sep + 1 :]
    else:
        own_args, passthrough = argv, []
    args = p.parse_args(own_args)

    if "--upgrade" in passthrough and not args.yes and not Globals.dry_run:
        stderr(
            "rbio install:compile: `--upgrade` upgrades every transitive dependency\n"
            "across the touched subprojects, which usually isn't what you want.\n"
            "Prefer `--package NAME` (repeatable) for targeted upgrades.\n"
            "If you really want a blanket upgrade, rerun with `--yes`."
        )
        return 1

    targets = _resolve_targets(args.subproject)
    if not targets:
        return 1

    pip_compile_args = list(passthrough)
    for name in args.package:
        pip_compile_args.extend(["--upgrade-package", name])

    passthrough_shell = " ".join(_shell_quote(a) for a in pip_compile_args)
    for in_path in targets:
        rel_in = str(in_path.relative_to(REPO_ROOT))
        rel_out = str(in_path.with_suffix(".txt").relative_to(REPO_ROOT))
        shell = (
            "set -e; "
            "pip install --quiet --root-user-action=ignore pip-tools && "
            # --annotation-style=line matches the format the historic
            # requirements.txt files were generated with (`# via X, Y`
            # on the same line as the package). Dropping this makes
            # every lock-file regeneration churn thousands of lines.
            f"pip-compile --annotation-style=line "
            f"--output-file {_shell_quote(rel_out)} "
            f"{passthrough_shell} {_shell_quote(rel_in)}"
        )
        rc = run(
            [
                "docker",
                "run",
                "--rm",
                "-v",
                f"{REPO_ROOT}:/work",
                "-w",
                "/work",
                # rpy2 (in common/requirements.in) probes for R at
                # metadata-generation time; ABI mode skips the probe so the
                # slim image can resolve the lock without R installed.
                "-e",
                "RPY2_CFFI_MODE=ABI",
                PIPTOOLS_IMAGE,
                "bash",
                "-c",
                shell,
            ],
            cwd=str(REPO_ROOT),
        )
        if rc != 0:
            stderr(f"rbio install:compile: pip-compile failed for {rel_in}")
            return rc
    return 0


def _resolve_targets(selectors):
    """Return absolute Paths to the .in files to compile.

    No selectors → every entry in REQUIREMENTS_INS. Otherwise each selector
    is matched against the known list by directory prefix or exact path.
    """
    known = {REPO_ROOT / rel: rel for rel in REQUIREMENTS_INS}
    if not selectors:
        return list(known.keys())

    picked = []
    for sel in selectors:
        hits = [
            p for p, rel in known.items() if rel == sel or rel.startswith(sel.rstrip("/") + "/")
        ]
        if not hits:
            stderr(f"rbio install:compile: no requirements.in matches '{sel}'")
            stderr("known subprojects:\n" + "".join(f"  {rel}\n" for rel in REQUIREMENTS_INS))
            return []
        picked.extend(hits)
    # Preserve order and dedupe.
    seen = set()
    return [p for p in picked if not (p in seen or seen.add(p))]


def _shell_quote(arg):
    """Minimal shell quoting for arguments embedded in a `bash -c` string."""
    import shlex

    return shlex.quote(arg)


COMMANDS = [
    ("install:deps", cmd_install_deps, "install missing OS-level dev dependencies"),
    ("install:venv", cmd_install_venv, "set up dr_env/ for IDE + pip-compile"),
    ("install:compile", cmd_install_compile, "regenerate requirements.txt from requirements.in"),
]
