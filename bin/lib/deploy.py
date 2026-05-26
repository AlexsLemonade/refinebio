"""deploy:* — env-driven deploy invocation."""

import argparse
import os

from lib._requirements import epilog_from, requires
from lib._runtime import run, stderr

_VALID_ENVS = ("dev", "staging", "prod")


@requires(
    env=["DEPLOY_ENV", "DEPLOY_TAG", "DEPLOY_USER"],
    env_optional={"BATCH_USE_ON_DEMAND_INSTANCES": ""},
    tools=["tfenv"],
)
def cmd_deploy_up(argv):
    p = argparse.ArgumentParser(
        prog="rbio deploy:up",
        description="Deploy refinebio. Reads target/version from environment.",
        epilog=epilog_from(
            cmd_deploy_up,
            "Locally: set DEPLOY_ENV=dev, DEPLOY_TAG=<git short SHA or label>, DEPLOY_USER=$USER\n"
            "CI:      set by remote_deploy.sh via env_vars file (from 1Password secrets)\n"
            "\n"
            "wraps: cd infrastructure && tfenv install && ./deploy.sh -e <env> -v <tag> -u <user> [-i <on_demand>]",
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.parse_args(argv)

    env_name = os.environ.get("DEPLOY_ENV")
    tag = os.environ.get("DEPLOY_TAG")
    user = os.environ.get("DEPLOY_USER")
    on_demand = os.environ.get("BATCH_USE_ON_DEMAND_INSTANCES")

    if env_name not in _VALID_ENVS:
        stderr(f"rbio deploy:up: DEPLOY_ENV must be one of {_VALID_ENVS} (got: {env_name!r})")
        return 1
    if not tag or not user:
        stderr("rbio deploy:up: DEPLOY_TAG and DEPLOY_USER must be set")
        return 1

    print(f"Deploying tag {tag} to {env_name} (user={user})")

    # tfenv reads infrastructure/.terraform-version (added in this PR).
    if (rc := run(["tfenv", "install"], cwd="infrastructure")) != 0:
        return rc

    cmd = ["./deploy.sh", "-e", env_name, "-v", tag, "-u", user]
    if on_demand:
        cmd.extend(["-i", on_demand])
    return run(cmd, cwd="infrastructure")


COMMANDS = [
    ("deploy:up", cmd_deploy_up, "deploy refinebio (env-driven)"),
]
