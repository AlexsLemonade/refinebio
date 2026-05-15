"""deploy:* — deploy-time templating + ops invoked from infrastructure/."""

import argparse
import os
import re
from pathlib import Path

from lib._runtime import REPO_ROOT, Globals, stderr

# Per-template RAM fan-out for the workers project. A template not listed here
# uses DEFAULT_WORKER_RAMS; one in NO_RAM_WORKER_TEMPLATES is rendered once
# without a RAM postfix.
WORKER_RAM_VARIANTS = {
    "downloader.json": [1024, 4096, 16384],
    "create_compendia.json": [30000, 950000],
    "create_quantpendia.json": [30000, 131000],
}
NO_RAM_WORKER_TEMPLATES = {
    "create_qn_target.json",
    "smasher.json",
    "tximport.json",
    "qn_dispatcher.json",
    "janitor.json",
}
DEFAULT_WORKER_RAMS = [2048, 4096, 8192, 12288, 16384, 32768, 65536]
SURVEYOR_RAMS = [1024, 4096, 16384]

# Image-env defaults applied if the caller didn't set them. Mirrors the
# `if [ -z "$X" ]; then export X=...; fi` block at the top of the old script.
DEFAULT_IMAGES = {
    "FOREMAN_DOCKER_IMAGE": "dr_foreman",
    "DOWNLOADERS_DOCKER_IMAGE": "dr_downloaders",
    "TRANSCRIPTOME_DOCKER_IMAGE": "dr_transcriptome",
    "SALMON_DOCKER_IMAGE": "dr_salmon",
    "SMASHER_DOCKER_IMAGE": "dr_smasher",
    "AFFYMETRIX_DOCKER_IMAGE": "dr_affymetrix",
    "ILLUMINA_DOCKER_IMAGE": "dr_illumina",
    "NO_OP_DOCKER_IMAGE": "dr_no_op",
    "COMPENDIA_DOCKER_IMAGE": "dr_compendia",
}

TEMPLATE_VAR_RE = re.compile(r"\$\{\{([^}]+)\}\}")


def cmd_deploy_format_batch(argv):
    args = _parse_args(argv)
    env = _build_env(args)
    output_dir = _resolve_output_dir(args)
    if Globals.dry_run:
        stderr(f"DRY RUN: would render {args.project} templates to {output_dir}")
        return 0
    output_dir.mkdir(parents=True, exist_ok=True)
    return _render_project(args.project, env, output_dir)


def _resolve_output_dir(args):
    if args.output_dir:
        return Path(args.output_dir)
    return _source_dir(args.project) / "batch-job-specs"


def _parse_args(argv):
    p = argparse.ArgumentParser(
        prog="rbio deploy:format-batch",
        description=(
            "Render AWS Batch job specifications / environment files for a "
            "project, substituting ${{NAME}} placeholders with values from "
            "infrastructure/prod_env + the current environment."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "-p",
        "--project",
        required=True,
        choices=["api", "workers", "surveyor", "foreman"],
        help="project whose templates to render",
    )
    p.add_argument(
        "-e",
        "--env",
        default="local",
        help="deploy environment (informational; passed through to templates). default: local",
    )
    p.add_argument(
        "-o",
        "--output-dir",
        help=(
            "where to write rendered files. default: <project>/batch-job-specs for "
            "workers/surveyor, or the project directory for api/foreman"
        ),
    )
    p.add_argument(
        "-v",
        "--system-version",
        default=os.environ.get("system_version", "latest"),
        help="image tag for default *_DOCKER_IMAGE values. default: latest",
    )
    return p.parse_args(argv)


def _source_dir(project):
    # surveyor templates live under foreman/.
    return REPO_ROOT / ("foreman" if project == "surveyor" else project)


def _build_env(args):
    """Layer the substitution env from lowest priority to highest. Later writes win."""
    env = dict(os.environ)

    # Image-tag defaults — only fill in if the caller hasn't already exported them.
    env.setdefault("DOCKERHUB_REPO", "ccdlstaging")
    for key, image in DEFAULT_IMAGES.items():
        env.setdefault(key, f"{image}:{args.system_version}")

    # CLI args bound to the names templates reference (`${{env}}`, `${{project}}`, ...).
    env["env"] = args.env
    env["project"] = args.project
    env["system_version"] = args.system_version

    # Mount path for the worker EBS volume — referenced by every worker template.
    env["VOLUME_DIR"] = "/var/ebs"

    # Highest priority: deploy-time values from terraform.
    _load_env_file(REPO_ROOT / "infrastructure" / "prod_env", env)

    return env


def _load_env_file(path, env):
    """Merge KEY=VALUE lines from `path` into `env` (skip blanks and #-comments)."""
    if not path.exists():
        return
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        key, _, value = line.partition("=")
        if key:
            env[key] = value


def _render_project(project, env, output_dir):
    """Dispatch by project: workers/surveyor fan out a template dir; api/foreman do one environment.tpl."""
    if project in ("workers", "surveyor"):
        templates_dir = _source_dir(project) / "batch-job-templates"
        if not templates_dir.is_dir():
            stderr(f"rbio deploy:format-batch: missing templates dir: {templates_dir}")
            return 1
        policy = _worker_policy if project == "workers" else _surveyor_policy
        for tpl in sorted(templates_dir.glob("*.tpl.json")):
            _render_template(tpl, output_dir, env, policy)
        return 0

    # api / foreman: a single environment.tpl -> <out>/environment
    tpl_path = _source_dir(project) / "environment.tpl"
    if not tpl_path.is_file():
        stderr(f"rbio deploy:format-batch: missing {tpl_path}")
        return 1
    _render_file(tpl_path, output_dir / "environment", env)
    return 0


def _render_template(template_path, output_dir, env, policy):
    """Render one template, fanning out by RAM size if the policy returns a non-empty list."""
    base = template_path.name[: -len(".tpl.json")]
    rams = policy(f"{base}.json")
    if not rams:
        _render_file(template_path, output_dir / f"{base}.json", env)
        return
    for ram in rams:
        scoped = {**env, "RAM": str(ram), "RAM_POSTFIX": f"_{ram}"}
        _render_file(template_path, output_dir / f"{base}_{ram}.json", scoped)


def _worker_policy(out_name):
    """Return the RAM variants for a workers template, or an empty list for single-render templates."""
    if out_name in NO_RAM_WORKER_TEMPLATES:
        return []
    return WORKER_RAM_VARIANTS.get(out_name, DEFAULT_WORKER_RAMS)


def _surveyor_policy(out_name):
    if out_name == "surveyor_dispatcher.json":
        return []
    return SURVEYOR_RAMS


def _render_file(src, dst, env):
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_text(_render(src.read_text(), env))
    print(f"Made {dst}")


def _render(text, env):
    """Substitute ${{NAME}} -> env[NAME]. Leave the placeholder unchanged if NAME is unset."""
    return TEMPLATE_VAR_RE.sub(lambda m: env.get(m.group(1), m.group(0)), text)


COMMANDS = [
    ("deploy:format-batch", cmd_deploy_format_batch, "render Batch job specs/env from templates"),
]
