#!/bin/bash

# Install git-crypt
cd
git clone https://github.com/AGWA/git-crypt.git
cd git-crypt
make
sudo make install
