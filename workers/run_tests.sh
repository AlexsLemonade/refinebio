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

# Make sure test Transcriptome Index is downloaded from S3.
index_dir="$volume_directory/processed/TEST/TRANSCRIPTOME_INDEX/"
gz_index_path="$index_dir/Homo_sapiens_short.tar.gz"
if [ ! -e "$gz_index_path" ]; then
    mkdir -p $index_dir
    wget -O $gz_index_path \
         https://s3.amazonaws.com/data-refinery-test-assets/Homo_sapiens_short.tar.gz
fi

# Make sure data for Salmon test is downloaded from S3.
rna_seq_test_raw_dir="$volume_directory/raw/TEST/SALMON/"
read_1_name="ERR003000_1.fastq.gz"
read_2_name="ERR003000_2.fastq.gz"
rna_seq_test_data_1="$rna_seq_test_raw_dir/$read_1_name"
rna_seq_test_data_2="$rna_seq_test_raw_dir/$read_2_name"
if [ ! -e "$rna_seq_test_data_data_1" ]; then
    mkdir -p $rna_seq_test_raw_dir
    wget -O $rna_seq_test_data_1 \
         https://s3.amazonaws.com/data-refinery-test-assets/$read_1_name
    wget -O $rna_seq_test_data_2 \
         https://s3.amazonaws.com/data-refinery-test-assets/$read_2_name
fi

docker build -t dr_worker_tests -f workers/Dockerfile.tests .

source common.sh
HOST_IP=$(get_ip_address)
DB_HOST_IP=$(get_docker_db_ip_address)
NOMAD_HOST_IP=$(get_docker_nomad_ip_address)
NOMAD_LINK=$(get_nomad_link_option)

docker run \
       --add-host=database:$DB_HOST_IP \
       --add-host=nomad:$NOMAD_HOST_IP \
       --env-file workers/environments/test \
       --volume $volume_directory:/home/user/data_store \
       --link drdb:postgres $NOMAD_LINK \
        -i dr_worker_tests python3 manage.py test --no-input "$@" # This runs everything
       # -i dr_worker_tests python3 manage.py test data_refinery_workers.processors.test_salmon.SalmonTestCase.test_success --no-input "$@" # This runs a specific test
       # Can also be called like ./workers/run_tests.sh data_refinery_workers.downloaders.test_sra.DownloadSraTestCase.test_aspera_downloader
