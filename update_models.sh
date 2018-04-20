#! /bin/bash
# Makes migrations and re-installs so Docker images update locally
./common/make_migrations.sh && cd common && python setup.py sdist
