#!/bin/sh
exec "$(dirname "$0")/../bin/rbio" ops:kill-jobs "$@"
