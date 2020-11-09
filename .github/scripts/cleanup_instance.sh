#!/bin/bash

rm -rf /usr/local/android-sdk-linux
rm -rf /usr/local/go
rm -rf /usr/local/android-ndk
rm -rf /home/ubuntu/nvm
rm -rf /home/ubuntu/.android
rm -rf /home/ubuntu/.phpenv
rm -rf /home/ubuntu/.rvm

# Suggestions from https://github.com/actions/virtual-environments/issues/709#issuecomment-615370473
sudo apt-get remove -y '^ghc-8.*'
sudo apt-get remove -y '^dotnet-.*'
sudo apt-get remove -y '^llvm-.*'
sudo apt-get remove -y 'php.*'
sudo apt-get remove -y azure-cli google-cloud-sdk hhvm google-chrome-stable firefox powershell mono-devel
sudo apt-get autoremove -y
sudo apt-get clean

rm -rf /usr/share/dotnet/
