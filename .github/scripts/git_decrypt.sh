#!/bin/bash

# Unlock encrypted files
# Temporary attempt to debug this.
pwd
cd
pwd
cd refinebio/.github/ || exit
git clean -f
openssl aes-256-cbc -md md5 -d -in git_crypt.key.enc -out git_crypt.key -k "$OPENSSL_KEY"
git-crypt unlock git_crypt.key
rm -f git_crypt.key
