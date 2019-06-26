#!/bin/sh

if ! type "virtualenv" > /dev/null; then
  pip3 install virtualenv
fi

virtualenv -p python3 dr_env

. dr_env/bin/activate

pip3 install pip-tools
