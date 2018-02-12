#!/bin/bash

if ! type "virtualenv" > /dev/null; then
  pip install virtualenv
fi

virtualenv dr_env

source dr_env/bin/activate

pip3 install pip-tools
