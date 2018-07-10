#!/usr/bin/env bash
../run_postgres.sh
../common/install_db_docker.sh
../common/make_migrations.sh