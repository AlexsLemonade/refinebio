"""debug:* — read-only diagnostics for env vars, external tools, and venv state.

Reads each command's @requires metadata and aggregates / scopes the view.
Never installs or modifies anything; use install:* commands for that.
"""

import argparse
import json
import os
import shutil
import subprocess

from lib._requirements import requires
from lib._runtime import REPO_ROOT, stderr

_EMPTY_REQS = {"env": [], "env_optional": {}, "env_alternatives": [], "tools": []}


def _all_commands():
    """Return [(name, fn, description, requirements)] across every namespace.

    Late-imported to avoid circular dependencies — debug.py needs to see
    every other namespace module's COMMANDS.
    """
    from lib import build, common, compose, db, deploy, dev, es, install, ops
    from lib.test import all as test_all, subprojects as test_subprojects, workers as test_workers

    namespaces = [
        build,
        common,
        compose,
        db,
        deploy,
        dev,
        es,
        install,
        ops,
        test_subprojects,
        test_workers,
        test_all,
    ]
    out = []
    for ns in namespaces:
        for name, fn, desc in ns.COMMANDS:
            req = getattr(fn, "_requirements", _EMPTY_REQS)
            out.append((name, fn, desc, req))
    return out


def _mask(name, val):
    """Mask all but the last 4 chars for secret-named vars; pass others through."""
    if not val:
        return ""
    if any(s in name.upper() for s in ("PASSWORD", "SECRET", "TOKEN", "KEY")):
        if len(val) <= 4:
            return "*" * len(val)
        return "*" * (len(val) - 4) + val[-4:]
    return val


def _alt_status(option):
    """Return ('ok'|'partial'|'empty', satisfied_count, total_count) for an alternative option."""
    set_count = sum(1 for v in option if os.environ.get(v))
    if set_count == len(option):
        return ("ok", set_count, len(option))
    if set_count > 0:
        return ("partial", set_count, len(option))
    return ("empty", 0, len(option))


@requires()
def cmd_debug_env(argv):
    p = argparse.ArgumentParser(
        prog="rbio debug:env",
        description="Show env var configuration. Reads @requires metadata across all commands.",
        epilog=(
            "Without --for, shows the flat aggregate of all required + optional\n"
            "env vars + alternative groups across every command.\n"
            "\n"
            "Pass --for <command> to scope to one command's requirements.\n"
            "\n"
            "Secret-named vars (containing PASSWORD/SECRET/TOKEN/KEY) are masked:\n"
            "all but the last 4 chars become *."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--for",
        dest="command",
        metavar="COMMAND",
        help="show requirements for a single command (e.g. 'ops:kill-jobs')",
    )
    args = p.parse_args(argv)

    all_cmds = _all_commands()
    if args.command:
        match = [c for c in all_cmds if c[0] == args.command]
        if not match:
            stderr(f"debug:env: unknown command '{args.command}'")
            return 1
        return _render_command_env(match[0])
    return _render_aggregate_env(all_cmds)


def _render_aggregate_env(all_cmds):
    required_by_cmd = {}
    optional_by_cmd = {}
    alt_groups = {}

    for name, _, _, req in all_cmds:
        for v in req["env"]:
            required_by_cmd.setdefault(v, []).append(name)
        for v, d in req["env_optional"].items():
            entry = optional_by_cmd.setdefault(v, (d, []))
            entry[1].append(name)
        for label, options in req["env_alternatives"]:
            entry = alt_groups.setdefault(label, (options, []))
            entry[1].append(name)

    if required_by_cmd:
        print("Required env vars (across all commands):")
        max_name = max(len(v) for v in required_by_cmd)
        for var in sorted(required_by_cmd):
            actual = os.environ.get(var)
            marker = "✓" if actual else "✗"
            val = _mask(var, actual) if actual else "unset"
            cmds = ", ".join(sorted(set(required_by_cmd[var])))
            print(f"  {marker} {var:<{max_name}}    {val:<30}  (needed by: {cmds})")
        print()

    for label in sorted(alt_groups):
        options, cmds = alt_groups[label]
        _render_alt_group(label, options, sorted(set(cmds)))
        print()

    if optional_by_cmd:
        print("Optional env vars (default in parens):")
        max_name = max(len(v) for v in optional_by_cmd)
        for var in sorted(optional_by_cmd):
            default, cmds = optional_by_cmd[var]
            actual = os.environ.get(var)
            marker = "✓" if actual else "·"
            val = _mask(var, actual) if actual else f"(default: {default or 'unset'})"
            cmds_str = ", ".join(sorted(set(cmds)))
            print(f"  {marker} {var:<{max_name}}    {val:<30}  (used by: {cmds_str})")

    return 0


def _render_command_env(entry):
    name, _, desc, req = entry
    print(f"# {name}  —  {desc}")
    print()
    if req["env"]:
        print("Required:")
        for v in req["env"]:
            actual = os.environ.get(v)
            marker = "✓" if actual else "✗"
            val = _mask(v, actual) if actual else "UNSET"
            print(f"  {marker} {v:<30}  {val}")
        print()
    for label, options in req["env_alternatives"]:
        _render_alt_group(label, options, cmds=None)
        print()
    if req["env_optional"]:
        print("Optional (defaults shown):")
        for v, d in req["env_optional"].items():
            actual = os.environ.get(v)
            marker = "✓" if actual else "·"
            val = _mask(v, actual) if actual else f"(default: {d or 'unset'})"
            print(f"  {marker} {v:<30}  {val}")
    return 0


def _render_alt_group(label, options, cmds):
    """Render a single alternatives group with status line."""
    header = f"{label} (need one)"
    if cmds:
        header += f"            (needed by: {', '.join(cmds)})"
    print(header)

    satisfied_option = None
    for option in options:
        status, _, _ = _alt_status(option)
        marker = {"ok": "✓", "partial": "✗", "empty": "·"}[status]
        opt_label = " + ".join(option)

        if len(option) == 1:
            v = option[0]
            actual = os.environ.get(v)
            val = _mask(v, actual) if actual else "unset"
            print(f"  {marker} {v:<40}    {val}")
        else:
            suffix = " incomplete" if status == "partial" else ""
            print(f"  {marker} {opt_label}{suffix}")
            for v in option:
                actual = os.environ.get(v)
                val = _mask(v, actual) if actual else "unset"
                extra = ""
                if status == "partial" and not actual:
                    extra = " (required for this option)"
                print(f"      {v}: {val}{extra}")

        if status == "ok":
            satisfied_option = opt_label

    if satisfied_option:
        print(f"  → status: OK (via {satisfied_option})")
    else:
        opt_strs = " or ".join(" + ".join(opt) for opt in options)
        print(f"  → status: MISSING — set {opt_strs}")


@requires()
def cmd_debug_deps(argv):
    p = argparse.ArgumentParser(
        prog="rbio debug:deps",
        description="Show external tool availability. Reads @requires metadata.",
        epilog=(
            "Without --for, shows the union of every command's tool deps.\n"
            "Pass --for <command> to scope to one command.\n"
            "\n"
            "Exits non-zero if any tool is missing — run `rbio install:deps` to install."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--for",
        dest="command",
        metavar="COMMAND",
        help="show tools required by a single command",
    )
    args = p.parse_args(argv)

    all_cmds = _all_commands()
    if args.command:
        match = [c for c in all_cmds if c[0] == args.command]
        if not match:
            stderr(f"debug:deps: unknown command '{args.command}'")
            return 1
        tools = sorted(set(match[0][3]["tools"]))
    else:
        tools = sorted({t for _, _, _, req in all_cmds for t in req["tools"]})

    if not tools:
        print("No external tools required.")
        return 0

    missing = []
    longest = max(len(t) for t in tools)
    print("External tools:")
    for tool in tools:
        path = shutil.which(tool)
        marker = "✓" if path else "✗"
        loc = path or "not found"
        print(f"  {marker} {tool:<{longest}}    {loc}")
        if not path:
            missing.append(tool)
    if missing:
        print(f"\n  → MISSING: {', '.join(missing)}.  Run `rbio install:deps`.")
        return 1
    return 0


_VENV_REQUIREMENT_FILES = [
    "api/requirements.txt",
    "foreman/requirements.txt",
    "common/requirements.txt",
    "workers/data_refinery_workers/processors/requirements.txt",
    "workers/data_refinery_workers/downloaders/requirements.txt",
]


@requires()
def cmd_debug_venv(argv):
    p = argparse.ArgumentParser(
        prog="rbio debug:venv",
        description=(
            "Check dr_env/ presence and whether installed packages match each "
            "subproject's requirements.txt."
        ),
        epilog="If anything is out of sync, run `rbio install:venv` to refresh.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.parse_args(argv)

    venv_path = REPO_ROOT / "dr_env"
    pip = venv_path / "bin" / "pip"
    if not pip.exists():
        print("dr_env/  ✗ not present. Run `rbio install:venv`.")
        return 1

    any_stale = False
    for rel in _VENV_REQUIREMENT_FILES:
        req_file = REPO_ROOT / rel
        result = subprocess.run(
            [
                str(pip),
                "install",
                "--dry-run",
                "--quiet",
                "--report",
                "-",
                "-r",
                str(req_file),
            ],
            capture_output=True,
            text=True,
            check=False,
        )
        if result.returncode != 0:
            stderr(f"debug:venv: pip install --dry-run failed for {rel}")
            stderr(result.stderr.strip())
            return 1
        try:
            report = json.loads(result.stdout)
        except json.JSONDecodeError:
            stderr(f"debug:venv: couldn't parse pip output for {rel}")
            return 1
        would_change = [
            f"{i['metadata']['name']}=={i['metadata']['version']}"
            for i in report.get("install", [])
        ]
        if would_change:
            any_stale = True
            print(f"  ✗ {rel}  ({len(would_change)} package(s) need install/upgrade)")
            for pkg in would_change[:5]:
                print(f"      {pkg}")
            if len(would_change) > 5:
                print(f"      ... and {len(would_change) - 5} more")
        else:
            print(f"  ✓ {rel}")

    if any_stale:
        print("\n→ Run `rbio install:venv` to refresh.")
        return 1
    print("\ndr_env/  ✓ all subproject requirements satisfied.")
    return 0


COMMANDS = [
    ("debug:env", cmd_debug_env, "show env var configuration"),
    ("debug:deps", cmd_debug_deps, "show external tool availability"),
    ("debug:venv", cmd_debug_venv, "check dr_env/ presence and freshness"),
]
