#!/bin/bash
python manage.py collectstatic --noinput
uwsgi --ini uwsgi.ini --http 8081
