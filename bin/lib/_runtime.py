"""Process-level runtime: globals, subprocess wrapper, repo paths."""

import shlex
import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent.parent


class Globals:
    verbose = False
    quiet = False
    dry_run = False
    no_color = False


def stderr(msg):
    print(msg, file=sys.stderr)


def run(cmd, cwd=None, env=None):
    """Run a subprocess; honor --verbose / --dry-run."""
    rendered = " ".join(shlex.quote(c) for c in cmd)
    if Globals.dry_run:
        print(f"DRY RUN: {rendered}", file=sys.stderr)
        return 0
    if Globals.verbose:
        print(f"+ {rendered}", file=sys.stderr)
    return subprocess.call(cmd, cwd=str(cwd or REPO_ROOT), env=env)
