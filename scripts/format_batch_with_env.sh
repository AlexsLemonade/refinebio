#!/bin/sh
exec "$(dirname "$0")/../bin/rbio" deploy:format-batch "$@"
