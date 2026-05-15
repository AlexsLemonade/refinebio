"""install:* — host bootstrap (system packages + venv + initial DB / ES)."""

import argparse
import os
import platform
import shutil
from dataclasses import dataclass

from lib._runtime import REPO_ROOT, Globals, run, stderr

TERRAFORM_VERSION = "0.13.5"

# Returned by an installer that succeeded partially but needs the user to log
# out and back in before the rest of the bootstrap can run (e.g. docker group
# membership). Translated to a clean exit-0 by cmd_install_all so the user
# isn't told the install "failed". EX_TEMPFAIL.
_RELOGIN_REQUIRED = 75


@dataclass
class PackageManager:
    """Resolved package-manager handle: how to install a regular package, optionally a cask."""

    install: list  # e.g. ["sudo", "apt-get", "install", "--assume-yes"]
    cask: list  # e.g. ["brew", "install", "--cask"]; None if N/A
    kind: str  # "brew" | "apt" | "custom"

    @property
    def is_brew(self):
        return self.kind == "brew"

    @property
    def is_apt(self):
        return self.kind == "apt"


def cmd_install_all(argv):
    p = argparse.ArgumentParser(
        prog="rbio install:all",
        description=(
            "Install the host dependencies needed for refine.bio dev "
            "(docker, terraform, pre-commit, jq, ...) then bring up the "
            "local postgres + elasticsearch and prime the database."
        ),
        epilog=(
            "env vars:\n"
            "  INSTALL_CMD   override the package install command (skips OS detection)\n"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.parse_args(argv)

    pm = _detect_package_manager()
    if pm is None:
        return 1
    rc = _install_host_deps(pm)
    if rc == _RELOGIN_REQUIRED:
        # Docker permissions need a fresh login before we can continue.
        # User already saw the explanation from _install_docker.
        return 0
    if rc != 0:
        return rc
    return _bootstrap_dev_stack()


def _detect_package_manager():
    """Pick a package-manager shape based on host. No side effects."""
    custom = os.environ.get("INSTALL_CMD")
    if custom:
        return PackageManager(install=custom.split(), cask=None, kind="custom")

    system = platform.system()
    if system == "Darwin":
        return PackageManager(
            install=["brew", "install"],
            cask=["brew", "install", "--cask"],
            kind="brew",
        )
    if system == "Linux":
        if not _which("apt-get"):
            stderr("rbio install:all: only apt-based Linux distros are auto-supported")
            stderr("set INSTALL_CMD='<your install command>' to override")
            return None
        return PackageManager(
            install=["sudo", "apt-get", "install", "--assume-yes"],
            cask=None,
            kind="apt",
        )

    stderr(f"rbio install:all: unsupported operating system: {system}")
    return None


def _install_host_deps(pm):
    # Bootstrap the package manager itself if needed, then refresh its index.
    if pm.is_brew and not _which("brew"):
        if not _install_homebrew():
            return 1
    if pm.is_apt:
        if (rc := run(["sudo", "apt-get", "update"])) != 0:
            return rc

    # Docker is special — adding the user to the docker group requires a
    # fresh login before the rest of the bootstrap can proceed.
    rc = _install_docker(pm)
    if rc == _RELOGIN_REQUIRED:
        return rc
    if rc != 0:
        stderr("rbio install:all: failed to install docker")
        return rc

    for name, installer in [
        ("python/pip", _install_python),
        ("terraform", _install_terraform),
        ("pre-commit", _install_pre_commit),
        ("jq", _install_jq),
        ("ip", _install_iproute),
    ]:
        if (rc := installer(pm)) != 0:
            stderr(f"rbio install:all: failed to install {name}")
            return rc
    return 0


def _install_homebrew():
    if not _confirm("Would you like to install Homebrew?"):
        stderr("rbio install:all: Homebrew is required on macOS")
        return False
    return run(
        [
            "/bin/bash",
            "-c",
            'NONINTERACTIVE=1 /bin/bash -c "$(curl -fsSL '
            'https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"',
        ]
    ) == 0


def _install_docker(pm):
    if _which("docker"):
        return 0
    print("Installing Docker...")
    if pm.is_brew:
        return run([*pm.cask, "docker"])

    if (rc := run([*pm.install, "docker.io"])) != 0:
        stderr("You must manually install Docker")
        return rc

    user = os.environ.get("USER")
    if not user:
        stderr("Skipping docker group setup: USER is unset.")
        return 0
    print("Fixing Docker permissions...")
    run(["sudo", "groupadd", "-f", "docker"])
    run(["sudo", "usermod", "-aG", "docker", user])
    print()
    print(
        "Logout and log back in to apply the permissions changes, "
        "then re-run 'rbio install:all'."
    )
    # Signal to the caller (cmd_install_all) that we need a fresh login before
    # the rest of the bootstrap can proceed.
    return _RELOGIN_REQUIRED


def _install_python(pm):
    # macOS gets python via brew/system, so no-op there.
    if pm.is_brew or _which("pip3"):
        return 0
    print("Installing Python and pip...")
    rc = run([*pm.install, "python3-pip"])
    if rc != 0:
        stderr("You must manually install Python and pip")
    return rc


def _install_terraform(pm):
    if _which("terraform"):
        return 0
    print("Installing Terraform...")
    if pm.is_brew:
        return run([*pm.install, "terraform"])
    if not (pm.is_apt or _confirm("Automatically install Terraform for amd64 Linux?")):
        stderr("You need to manually install Terraform before continuing...")
        return 1
    if (rc := run([*pm.install, "unzip"])) != 0:
        return rc
    return _download_terraform_amd64()


def _download_terraform_amd64():
    zipname = f"terraform_{TERRAFORM_VERSION}_linux_amd64.zip"
    url = f"https://releases.hashicorp.com/terraform/{TERRAFORM_VERSION}/{zipname}"
    for step in (
        ["curl", "-0sL", url, "-o", zipname],
        ["sudo", "unzip", "-d", "/usr/bin", zipname],
        ["sudo", "chmod", "a+rx", "/usr/bin/terraform"],
    ):
        if (rc := run(step)) != 0:
            return rc
    if not Globals.dry_run:
        os.remove(zipname)
    return 0


def _install_pre_commit(pm):
    if _which("pre-commit"):
        return 0
    msg = (
        "Automatically install pre-commit? Note: This will install all the "
        "required dependencies (black, isort, etc) using an additional ~185MB "
        "of disk space."
    )
    if not (pm.is_apt or _confirm(msg)):
        print("Skipping installation of pre-commit.")
        return 0
    print("Installing pre-commit...")
    if (rc := run([*pm.install, "shellcheck"])) != 0:
        return rc
    if (rc := run(["pip3", "install", "pre-commit"])) != 0:
        return rc
    return run(["pre-commit", "install"], cwd=REPO_ROOT)


def _install_jq(pm):
    if _which("jq"):
        return 0
    print("Installing jq...")
    rc = run([*pm.install, "jq"])
    if rc != 0:
        stderr("You must manually install jq")
    return rc


def _install_iproute(pm):
    if _which("ip"):
        return 0
    pkg = "iproute2mac" if pm.is_brew else "iproute2"
    rc = run([*pm.install, pkg])
    if rc != 0:
        stderr(f"You must manually install {pkg}")
    return rc


def _bootstrap_dev_stack():
    """Bring up postgres + ES, init DB, build search index, create venv, sync common."""
    rbio = str(REPO_ROOT / "bin" / "rbio")
    print("Starting postgres and installing the database...")
    if (rc := run([rbio, "dev:up", "-s", "postgres"])) != 0:
        return rc
    if (rc := run([rbio, "db:init"])) != 0:
        return rc

    print("Starting elasticsearch and building the ES Indexes...")
    if (rc := run([rbio, "dev:up", "-s", "elasticsearch"])) != 0:
        return rc
    if (rc := run([rbio, "es:rebuild"])) != 0:
        return rc

    print("Creating virtual environment...")
    if (rc := cmd_install_venv([])) != 0:
        return rc
    print("Run `source dr_env/bin/activate` to activate the virtual environment.")

    print("Updating common dependencies...")
    venv_path = REPO_ROOT / "dr_env" / "bin"
    venv_env = {**os.environ, "PATH": f"{venv_path}:{os.environ.get('PATH', '')}"}
    return run([rbio, "common:update"], env=venv_env)


def cmd_install_venv(argv):
    p = argparse.ArgumentParser(
        prog="rbio install:venv",
        description="Create the dr_env virtualenv at the repo root and install pip-tools.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.parse_args(argv)

    if not _which("virtualenv"):
        if (rc := run(["pip3", "install", "virtualenv"])) != 0:
            return rc
    if (rc := run(["virtualenv", "-p", "python3", "dr_env"])) != 0:
        return rc
    pip = str(REPO_ROOT / "dr_env" / "bin" / "pip3")
    return run([pip, "install", "pip-tools"])


def _which(binary):
    return shutil.which(binary) is not None


def _confirm(prompt):
    """Interactive y/N. Returns True for 'y'. Falsy under --dry-run."""
    if Globals.dry_run:
        stderr(f"DRY RUN: would prompt: {prompt}")
        return False
    answer = input(f"{prompt} [y/N] ").strip().lower()
    return answer == "y"


COMMANDS = [
    ("install:all", cmd_install_all, "install host dependencies + bootstrap the dev stack"),
    ("install:venv", cmd_install_venv, "create the dr_env virtualenv"),
]
