#!/bin/sh

# This script should always run from the context of the root directory of
# the project it is building.
script_directory="$(perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory/../.." || exit

set -e

files=`git diff --staged --name-only --diff-filter=d -- "*.py"`

seed-isort-config

for file in $files; do
  isort -y $file
  black --line-length 100 $file
  git add $file
done
