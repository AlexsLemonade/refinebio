#!/bin/sh
# install/setup_venv.sh — create dr_env/ with all subproject requirements.
#
# Invoked by `rbio install:venv`. The venv exists for two purposes:
#   1. IDE language server / autocomplete / jump-to-def across the monorepo
#   2. `rbio common:build-sdist` on PEP-668 distros where system python3
#      lacks setuptools

set -e

cd "$(dirname "$0")/.."   # repo root

if [ ! -d dr_env ]; then
    python3 -m venv dr_env
fi

. dr_env/bin/activate

pip install --quiet --upgrade pip pip-tools

pip install --quiet -r api/requirements.txt
pip install --quiet -r foreman/requirements.txt
pip install --quiet -r common/requirements.txt
pip install --quiet -r workers/data_refinery_workers/processors/requirements.txt
pip install --quiet -r workers/data_refinery_workers/downloaders/requirements.txt

# Editable install of common so `import data_refinery_common` resolves to
# the source tree (IDEs follow imports correctly).
pip install --quiet -e ./common

echo
echo "Done. Point your IDE at $(pwd)/dr_env/bin/python3"
echo "Activate in your shell:  . dr_env/bin/activate"
