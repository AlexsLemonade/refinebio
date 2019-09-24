#!/bin/sh

# This script is executed once the container is created and the filesystem ready

pip install --editable common

# https://code.visualstudio.com/docs/remote/containers#_sharing-git-credentials-with-your-container
mkdir -p ~/.ssh \
    && cp -r ~/.ssh-localhost/* ~/.ssh \
    && chmod 700 ~/.ssh \
    && chmod 600 ~/.ssh/*

