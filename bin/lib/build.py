"""build — docker image builds via buildx bake."""

import argparse
import json
import subprocess
import time

from lib._requirements import requires
from lib._runtime import REPO_ROOT, run, stderr
from lib.common import rebuild_common_sdist

RETRY_BACKOFF_SECONDS = 5


@requires(tools=["docker"])
def cmd_build(argv):
    needs_help = any(a in ("-h", "--help") for a in argv)
    epilog = (
        "wraps:\n"
        "  rbio common:build-sdist  (mtime-gated; skipped when common source unchanged)\n"
        "  docker buildx bake -f docker-bake.hcl <targets>"
    )
    if needs_help:
        section = render_bake_targets_section()
        if section:
            epilog = section + "\n\n" + epilog

    p = argparse.ArgumentParser(
        prog="rbio build",
        description="Build one or more docker images.",
        epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "targets",
        nargs="*",
        default=["default"],
        metavar="TARGET",
        help="image(s) or group(s) to build (see list below). default: 'default' group.",
    )
    p.add_argument(
        "--push", action="store_true", help="push to registry after build (requires creds)"
    )
    p.add_argument(
        "--load",
        action="store_true",
        help="load into local docker (default for single-target builds)",
    )
    p.add_argument("--no-cache", action="store_true", help="rebuild from scratch")
    p.add_argument("--pull", action="store_true", help="pull latest base images before building")
    p.add_argument(
        "--print",
        dest="print_only",
        action="store_true",
        help="resolve and print HCL; build nothing",
    )
    p.add_argument(
        "--retries",
        type=int,
        default=0,
        metavar="N",
        help=(
            "retry up to N times after a failed bake (useful in CI for "
            "transient registry push errors). default: 0."
        ),
    )
    args = p.parse_args(argv)

    # Bake's worker/foreman/api targets COPY common/dist/data?refinery?common-*
    # at image-build time. Make sure that tarball reflects current source
    # before bake reads it. The helper short-circuits when nothing changed,
    # preserving docker's COPY layer cache for downstream images.
    if not args.print_only:
        if (rc := rebuild_common_sdist()) != 0:
            return rc

    cmd = ["docker", "buildx", "bake", "-f", "docker-bake.hcl"]
    if args.push:
        cmd.append("--push")
    if args.load:
        cmd.append("--load")
    if args.no_cache:
        cmd.append("--no-cache")
    if args.pull:
        cmd.append("--pull")
    if args.print_only:
        cmd.append("--print")
    cmd.extend(args.targets)

    attempts = args.retries + 1
    for attempt in range(1, attempts + 1):
        rc = run(cmd)
        if rc == 0:
            return 0
        if attempt < attempts:
            stderr(
                f"rbio build: attempt {attempt}/{attempts} failed (exit {rc}), "
                f"retrying in {RETRY_BACKOFF_SECONDS}s..."
            )
            time.sleep(RETRY_BACKOFF_SECONDS)
    if attempts > 1:
        stderr(f"rbio build: failed after {attempts} attempts")
    return rc


# Keep in sync with docker-bake.hcl group definitions.
GROUP_DESCRIPTIONS = {
    "all": "every image",
    "default": "all images except affymetrix (long build time)",
    "deploy": "production images only (skips the local-dev variants)",
    "workers": "only the worker images",
}


def render_bake_targets_section():
    """Return formatted groups + targets sections as a string. Empty on failure."""
    try:
        result = subprocess.run(
            [
                "docker",
                "buildx",
                "bake",
                "-f",
                "docker-bake.hcl",
                "--list",
                "type=targets,format=json",
            ],
            capture_output=True,
            text=True,
            cwd=str(REPO_ROOT),
        )
        if result.returncode != 0:
            return ""
        entries = json.loads(result.stdout)
    except (FileNotFoundError, json.JSONDecodeError):
        return ""

    groups = sorted(e["name"] for e in entries if e.get("group"))
    targets = sorted(e["name"] for e in entries if not e.get("group"))
    if not (groups or targets):
        return ""

    lines = []

    if groups:
        lines.append("groups:")
        col = max(len(g) for g in groups) + 2
        for g in groups:
            lines.append(f"  {g.ljust(col)}{GROUP_DESCRIPTIONS.get(g, '')}")

    if targets:
        if groups:
            lines.append("")
        lines.append("images:")
        col_width = max(len(t) for t in targets) + 2
        n_cols = max(1, 78 // col_width)
        n_rows = (len(targets) + n_cols - 1) // n_cols
        for r in range(n_rows):
            row = []
            for c in range(n_cols):
                idx = c * n_rows + r
                if idx < len(targets):
                    row.append(targets[idx].ljust(col_width))
            lines.append("  " + "".join(row).rstrip())

    return "\n".join(lines)


COMMANDS = [
    ("build", cmd_build, "build one or more docker images"),
]
