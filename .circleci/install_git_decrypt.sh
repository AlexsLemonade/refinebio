#!/bin/bash

# Install git-crypt
cd || exit
git clone https://github.com/AGWA/git-crypt.git
cd git-crypt || exit
make
sudo make install
