#!/bin/sh

# Exit on error
set -e

# Config variables
TERRAFORM_VERSION="0.12.3"

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

print_description() {
    cat <<EOF
Automatically installs all of the required programs for refine.bio development.
EOF
}

print_usage() {
    cat <<EOF
Usage: [ENVIRONMENT_VARIABLE(S)] install_all.sh [OPTION(S)]

Options:
  -h		Print this help message
  -v		Verbose mode: show the output of package management commands

Environment variables:
  INSTALL_CMD   Use a custom package management command for package installations.
  		(Make sure it does not expect user input unless you run this script with -v)
                Note: It is pretty likely that at least one package is called something different
                in your repositories, so you may have to manually install some dependencies.
EOF
}

confirm() {
    printf "%s [y/N] " "$1"
    read -r confirmation
    if ! [ "$confirmation" = "y" ]; then
	echo "Confirmation failure" >&2
        exit 1
    fi
}

while getopts "hv" opt; do
    case $opt in
        h)
            print_description
            echo
            print_usage
            exit 0
            ;;
	v)
	    OUTPUT="/dev/stdout"
	    ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            print_usage >&2
            exit 1
            ;;
    esac
done

# Unless output was set to stdout by the verbose flag, set it to /dev/null
# to hide the stdout of package management commands
if [ -z "$OUTPUT" ]; then
    OUTPUT="/dev/null"
fi

if [ -z "$INSTALL_CMD" ]; then
    case "$(uname)" in
        "Darwin")
            if ! command -v brew >/dev/null; then
                confirm "Would you like to install Homebrew?"
                /usr/bin/ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)"
            fi

            INSTALL_CMD="brew install"
            INSTALL_CASK_CMD="brew cask install"
            BREW=true
            ;;
        "Linux")
            if command -v apt >/dev/null; then
                sudo apt-get update
                INSTALL_CMD="sudo apt-get install --assume-yes"
                APT=true
            else
                echo "Your Linux distribution is not officially supported," >&2
                echo "but it *should* be able to run the required services. You need to manually" >&2
                echo "install dependencies or give the command to install dependencies with \$INSTALL_CMD." >&2
                exit 1
            fi
            ;;
        *)
            echo "$(uname) is an unsupported operating system." >&2
            echo "You can try to provide a package manager command with \$INSTALL_CMD," >&2
            echo "but your mileage may vary." >&2
            exit 1
            ;;
    esac
fi

if ! command -v docker > /dev/null; then
    echo "Installing Docker..."

    # On macOS, install docker desktop with Homebrew cask
    if [ $BREW ]; then
        $INSTALL_CASK_CMD docker > $OUTPUT
    else
        $INSTALL_CMD docker.io > $OUTPUT || (echo "You must manually install docker" && exit 1)

        echo "Fixing docker permissions..."
        sudo groupadd -f docker
        sudo usermod -aG docker "$USER"

	echo
	echo "Logout and log back in to apply the permissions changes, then execute this script again."
	exit 0
    fi
fi

if ! command -v pip3 > /dev/null && ! [ $BREW ]; then # Don't reinstall python on macOS
    echo "Installing python and pip..."
    $INSTALL_CMD python3-pip > $OUTPUT || (echo "You must manually install python and pip" && exit 1)
fi

if ! command -v terraform > /dev/null; then
    echo "Installing terraform..."
    if [ $BREW ]; then
        $INSTALL_CMD terraform > $OUTPUT
    elif [ $APT ] || confirm "Would you like to automatically install Terraform for amd64 linux?"; then
        $INSTALL_CMD unzip > $OUTPUT
        curl -0s "https://releases.hashicorp.com/terraform/${TERRAFORM_VERSION}/terraform_${TERRAFORM_VERSION}_linux_amd64.zip" \
             > "terraform_${TERRAFORM_VERSION}_linux_amd64.zip"
        sudo unzip -d /usr/bin "terraform_${TERRAFORM_VERSION}_linux_amd64.zip"
        sudo chmod a+rx /usr/bin/terraform
	rm "terraform_${TERRAFORM_VERSION}_linux_amd64.zip"
    else
        echo "You need to manually install Terraform before continuing..." >&2
        exit 1
    fi
fi

if ! command -v nomad > /dev/null; then
    echo "Installing nomad..."
    if [ $BREW ]; then
        $INSTALL_CMD nomad > $OUTPUT
    else
        sudo ./install_nomad.sh
    fi
fi

if ! command -v pre-commit > /dev/null; then
    message="Would you like to automatically install pre-commit? \
Note: This will install all the required dependencies (black, isort, etc) \
using an additional ~185MB of disk space."
    if [ $APT ] || confirm $message; then
        echo "Installing pre-commit..."
        pip3 install pre-commit
        pre-commit install
    else
        echo "Skipping installation of pre-commit."
    fi
fi

if ! command -v git-crypt > /dev/null; then
    echo "Installing git-crypt..."
    $INSTALL_CMD git-crypt > $OUTPUT || (echo "You must manually install git-crypt" && exit 1)
fi
if ! command -v jq > /dev/null; then
    echo "Installing jq..."
    $INSTALL_CMD jq > $OUTPUT || (echo "You must manually install jq" && exit 1)
fi

if ! command -v ip > /dev/null; then
    if [ $BREW ]; then
        $INSTALL_CMD iproute2mac > $OUTPUT
    else
        $INSTALL_CMD iproute2 > $OUTPUT || (echo "You must manually install iproute2" && exit 1)
    fi
fi

if ! test -e ../.git/hooks/pre-commit > /dev/null; then
    if confirm "Would you like to configure a pre-commit hook for this project to run black (code formatter) prior to each commit?"; then
        echo "Installing pre-commit hook to auto-format code."
        cp hooks/autoformat.sh ../.git/hooks/pre-commit
        chmod +x ../.git/hooks/pre-commit
    else
        echo "Not installing pre-commit hook to auto-format code."
    fi
fi

echo "Starting postgres and installing the database..."
./run_postgres.sh > $OUTPUT
./install_db_docker.sh > $OUTPUT

echo "Starting elasticsearch and building the ES Indexes..."
./run_es.sh > $OUTPUT
./rebuild_es_index.sh > $OUTPUT

echo "Creating virtual environment..."
./create_virtualenv.sh > $OUTPUT
echo "Run \`source dr_env/bin/activate\` to activate the virtual environment."

echo "Updating common dependencies..."
# Source the virtual environment first
. ../dr_env/bin/activate
./update_models.sh > $OUTPUT
