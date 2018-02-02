#!/bin/bash

# Script for executing Django PyUnit tests within a Docker container.

# This script should always run as if it were being called from
# the directory it lives in.
script_directory=`perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0"`
cd $script_directory

# However in order to give Docker access to all the code we have to
# move up a level
cd ..

# Set up the test data volume directory if it does not already exist
volume_directory="$script_directory/test_volume"
if [ ! -d "$volume_directory" ]; then
    mkdir $volume_directory
    chmod 775 $volume_directory
fi

# In order to get the size of certain test data within the limit of
# Git LFS we have to use xz compression instead of gz. Therefore
# before we can run the tests we need to make sure the files are
# decompressed.
index_dir="$volume_directory/processed/TEST/TRANSCRIPTOME_INDEX/"
gz_index_path="$index_dir/Homo_sapiens_short.tar.gz"
if [ ! -e "$gz_index_path" ]; then
    # xz_index_path="$index_dir/Homo_sapiens_short.tar.xz"
    # unxz -k --stdout $xz_index_path | gzip > $gz_index_path
    mkdir -p workers/test_volume/processed/TEST/TRANSCRIPTOME_INDEX
    wget -O workers/test_volume/processed/TEST/TRANSCRIPTOME_INDEX/Homo_sapiens_short.tar.gz \
         https://s3.amazonaws.com/data-refinery-test-assets/Homo_sapiens_short.tar.gz
fi

docker build -t dr_worker_tests -f workers/Dockerfile.scan_tests .

source common.sh
DB_HOST_IP=$(get_docker_db_ip_address)
HOST_IP=$(get_ip_address)
NOMAD_HOST_IP=$(get_docker_nomad_ip_address)

docker run \
       --add-host=database:$DB_HOST_IP \
       --add-host=nomad:$NOMAD_HOST_IP \
       --env-file workers/environments/test \
       --volume $volume_directory:/home/user/data_store \
       --link drdb:postgres \
       --link nomad:nomad \
        -i dr_worker_tests python3 manage.py test --no-input --tag=scan "$@"  --verbosity=3  # This runs everything
       # -i dr_worker_tests python3 manage.py test data_refinery_workers.processors.test_salmon.SalmonTestCase.test_success --no-input "$@" # This runs a specific test
       # Can also be called like ./workers/run_tests.sh data_refinery_workers.downloaders.test_sra.DownloadSraTestCase.test_aspera_downloader 
