#!/bin/sh

set -e

REPO=$(echo "ghcr.io/$GITHUB_REPOSITORY" | tr '[:upper:]' '[:lower:]')
if [ -z "$IMAGES" ]; then
    echo "Error: must put images to pull in \$IMAGES" >&2
    exit 1
fi

for image in $IMAGES; do
    PACKAGE="$REPO/dr_$image"
    # Only pull the package if it already exists
    (docker pull "$PACKAGE" && docker tag "$PACKAGE" "ccdlstaging/dr_$image") || true
done
