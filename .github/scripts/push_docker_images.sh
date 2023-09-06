#!/bin/sh

set -e

if [ -z "$IMAGES" ]; then
    echo "Error: must put images to pull in \$IMAGES" >&2
    exit 1
fi

DOCKERHUB_REPO=$(echo "ghcr.io/$GITHUB_REPOSITORY" | tr '[:upper:]' '[:lower:]')

for image in $IMAGES; do
    PACKAGE="$DOCKERHUB_REPO/dr_$image"
    docker tag "ccdlstaging/dr_$image" "$PACKAGE"
    docker push "$PACKAGE"
done
