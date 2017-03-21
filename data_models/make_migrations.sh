#!/bin/bash

env $(cat environments/dev | xargs) ./manage.py makemigrations data_refinery_models
env $(cat environments/dev | xargs) ./manage.py migrate
