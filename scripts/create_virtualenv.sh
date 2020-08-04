#!/bin/sh

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# However in order to give virtualenv access to all the code we have to
# move up a level
cd ..

if ! command -v virtualenv >/dev/null; then
  pip3 install virtualenv
fi

virtualenv -p python3 dr_env

. dr_env/bin/activate

pip3 install pip-tools
