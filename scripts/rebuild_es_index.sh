#!/bin/sh

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(cd "$(dirname "$0")" || exit; pwd)"
cd "$script_directory" || exit

./run_manage.sh -i api_local -s api search_index --rebuild -f
