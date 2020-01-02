#!/bin/sh

set -e

files=`git diff --staged --name-only --diff-filter=d -- "*.py"`

for file in $files; do
  black --line-length 100 $file
  git add $file
done
