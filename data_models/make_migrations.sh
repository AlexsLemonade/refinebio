#!/bin/bash

env $(cat environments/dev | xargs) python3.5 manage.py makemigrations data_refinery_models
env $(cat environments/dev | xargs) python3.5 manage.py migrate
