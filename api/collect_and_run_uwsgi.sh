#!/bin/sh
python manage.py collectstatic --noinput --clear
uwsgi --ini uwsgi.ini
