#!/bin/bash -e

# This script shouldn't have to be run very often. It will be needed
# if we need to re-init git-crypt, as it generates an encrypted key
# file that can be added to the repo, along with an RSA key which can
# be used to decrypt it.

# Note that you will need to have access to git-crypt to run this.

openssl genrsa -out private.pem 4096

key=$(openssl rsa -in private.pem -outform PEM -pubout)

rm private.pem

git-crypt export-key git_crypt.key

openssl aes-256-cbc -md md5 -e -in git_crypt.key -out git_crypt.key.enc -k "$key"

rm git_crypt.key

echo "The encrypted git-crypt key file is: git_crypt.key.enc"
echo "The key for that file is:"
echo "$key"
echo ""
echo "It can be unlocked with:"
echo ""
echo "OPENSSL_KEY=\"$key\""
echo 'openssl aes-256-cbc -md md5 -pbkdf2 -iter 100000 -d -in git_crypt.key.enc -out git_crypt.key -k "$OPENSSL_KEY"'
