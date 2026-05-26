#!/bin/sh
# install/install_deps.sh — idempotent OS-level dev dependency installer.
#
# Invoked by `rbio install:deps`. Detects macOS/Linux, dispatches to the
# right package manager, skips already-installed tools.

set -e

YES=0
if [ "${1:-}" = "--yes" ] || [ "${1:-}" = "-y" ]; then
    YES=1
fi

case "$(uname)" in
    Darwin) PM=brew ;;
    Linux)
        if command -v apt-get >/dev/null; then
            PM=apt
        else
            echo "Unsupported Linux distribution (no apt-get). Install dependencies manually." >&2
            exit 1
        fi
        ;;
    *)
        echo "Unsupported OS: $(uname). Install dependencies manually." >&2
        exit 1
        ;;
esac

confirm() {
    if [ "$YES" -eq 1 ]; then return 0; fi
    printf "%s [y/N] " "$1"
    read -r ans
    [ "$ans" = "y" ] || [ "$ans" = "yes" ]
}

# Tools list. install_simple installs the same name via brew/apt; the docker
# and tfenv installs need bespoke handling (see below).
install_simple() {
    cmd=$1; brew_pkg=$2; apt_pkg=$3
    if command -v "$cmd" >/dev/null; then return 0; fi
    if ! confirm "Install $cmd?"; then
        echo "Skipping $cmd."
        return 0
    fi
    case "$PM" in
        brew) brew install "$brew_pkg" ;;
        apt)  sudo apt-get update -qq && sudo apt-get install -y "$apt_pkg" ;;
    esac
}

install_simple jq         jq         jq
install_simple shellcheck shellcheck shellcheck
install_simple pre-commit pre-commit pre-commit

# Docker
if ! command -v docker >/dev/null; then
    if confirm "Install Docker?"; then
        case "$PM" in
            brew) brew install --cask docker ;;
            apt)
                sudo apt-get update -qq
                sudo apt-get install -y docker.io
                sudo groupadd -f docker
                sudo usermod -aG docker "$USER"
                echo
                echo "Log out and log back in for docker group membership to take effect,"
                echo "then re-run this script."
                exit 0
                ;;
        esac
    else
        echo "Skipping Docker. Most rbio commands will not work without it."
    fi
fi

# AWS CLI
if ! command -v aws >/dev/null; then
    if confirm "Install AWS CLI?"; then
        case "$PM" in
            brew) brew install awscli ;;
            apt)
                tmpdir=$(mktemp -d)
                curl -fsSL -o "$tmpdir/awscliv2.zip" \
                    "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip"
                unzip -q "$tmpdir/awscliv2.zip" -d "$tmpdir"
                sudo "$tmpdir/aws/install"
                rm -rf "$tmpdir"
                ;;
        esac
    else
        echo "Skipping AWS CLI. rbio ops:* and deploy paths will not work."
    fi
fi

# tfenv — brew has it natively; apt requires a git clone + PATH update.
if ! command -v tfenv >/dev/null; then
    if confirm "Install tfenv (terraform version manager)?"; then
        case "$PM" in
            brew) brew install tfenv ;;
            apt)
                if [ ! -d "$HOME/.tfenv" ]; then
                    git clone --depth 1 https://github.com/tfutils/tfenv "$HOME/.tfenv"
                fi
                if ! grep -q 'tfenv/bin' "$HOME/.profile" 2>/dev/null; then
                    # SC2016: literal string — $HOME and $PATH are expanded
                    # when ~/.profile is sourced, not when this script runs.
                    # shellcheck disable=SC2016
                    echo 'export PATH="$HOME/.tfenv/bin:$PATH"' >> "$HOME/.profile"
                    echo "Added tfenv to ~/.profile."
                fi
                echo "Reload your shell or run: export PATH=\"\$HOME/.tfenv/bin:\$PATH\""
                ;;
        esac
    else
        echo "Skipping tfenv. rbio deploy:up will not work without it."
    fi
fi

echo
echo "Done. Run \`rbio debug:deps\` to verify."
