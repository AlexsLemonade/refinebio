#!/bin/sh

# This script is executed once the container is created and the filesystem ready

pip install --editable common

# https://code.visualstudio.com/docs/remote/containers#_sharing-git-credentials-with-your-container
mkdir -p ~/.ssh \
    && cp -r ~/.ssh-localhost/* ~/.ssh \
    && chmod 700 ~/.ssh \
    && chmod 600 ~/.ssh/*

git config commit.gpgsign false

# https://click.palletsprojects.com/en/7.x/python3/#python-3-surrogate-handling
export LC_ALL=C.UTF-8
export LANG=C.UTF-8