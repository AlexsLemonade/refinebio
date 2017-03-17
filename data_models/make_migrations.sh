#!/bin/bash

env $(cat environments/dev | xargs) ./manage.py makemigrations bioinformatics_mill_models
env $(cat environments/dev | xargs) ./manage.py migrate
