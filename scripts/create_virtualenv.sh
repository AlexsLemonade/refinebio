#!/bin/sh
exec "$(dirname "$0")/../bin/rbio" install:venv "$@"
