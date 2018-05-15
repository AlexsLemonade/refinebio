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
    chmod -R a+rwX $volume_directory
fi

test_data_repo="https://s3.amazonaws.com/data-refinery-test-assets"

# Make sure test Transcriptome Index is downloaded from S3 for salmon tests.
index_dir="$volume_directory/processed/TEST/TRANSCRIPTOME_INDEX"
index_tarball="Homo_sapiens_short.tar.gz"
gz_index_path="$index_dir/$index_tarball"
if [ ! -e "$gz_index_path" ]; then
    mkdir -p $index_dir
    echo "Downloading Salmon index for Salmon tests."
    wget -q -O $gz_index_path \
         "$test_data_repo/$index_tarball"
fi

# Make sure data for Salmon test is downloaded from S3.
rna_seq_test_raw_dir="$volume_directory/raw/TEST/SALMON"
read_1_name="ERR003000_1.fastq.gz"
read_2_name="ERR003000_2.fastq.gz"
rna_seq_test_data_1="$rna_seq_test_raw_dir/$read_1_name"
rna_seq_test_data_2="$rna_seq_test_raw_dir/$read_2_name"
if [ ! -e "$rna_seq_test_data_1" ]; then
    mkdir -p $rna_seq_test_raw_dir
    echo "Downloading ERR003000_1.fastq.gz for Salmon tests."
    wget -q -O $rna_seq_test_data_1 \
         "$test_data_repo/$read_1_name"
    echo "Downloading ERR003000_2.fastq.gz for Salmon tests."
    wget -q -O $rna_seq_test_data_2 \
         "$test_data_repo/$read_2_name"
fi

# Make sure pcl for test is downloaded from S3
pcl_name="GSM1426071_CD_colon_active_1.pcl"
pcl_name2="GSM45588.pcl"
pcl_test_raw_dir="$volume_directory/raw/TEST/pcl"
pcl_test_data_1="$pcl_test_raw_dir/$pcl_name"
pcl_test_data_2="$pcl_test_raw_dir/$pcl_name2"
if [ ! -e "$pcl_test_data_1" ]; then
    mkdir -p $pcl_test_raw_dir
    echo "Downloading pcl for tests."
    wget -q -O $pcl_test_data_1 \
         "$test_data_repo/$pcl_name"
fi
if [ ! -e "$pcl_test_data_2" ]; then
    echo "Downloading Non-Brainarray pcl for tests."
    wget -q -O $pcl_test_data_2 \
         "$test_data_repo/$pcl_name2"
fi

#
# Make sure PCVLs for test is downloaded from S3
pcl_name="GSE22427_non-normalized.PCL"
pcl_name2="GSM466597_95899_agilent.PCL"
pcl_test_raw_dir="$volume_directory/raw/TEST/PCL"
pcl_test_data_1="$pcl_test_raw_dir/$pcl_name"
pcl_test_data_2="$pcl_test_raw_dir/$pcl_name2"
if [ ! -e "$pcl_test_data_1" ]; then
    mkdir -p $pcl_test_raw_dir
    echo "Downloading pcl for tests."
    wget -q -O $pcl_test_data_1 \
         "$test_data_repo/$pcl_name"
fi
if [ ! -e "$pcl_test_data_2" ]; then
    echo "Downloading pcl for tests."
    wget -q -O $pcl_test_data_2 \
         "$test_data_repo/$pcl_name2"
fi

# Download salmontools test data
rm -rf $volume_directory/salmontools/
git clone git@github.com:dongbohu/salmontools_tests.git $volume_directory/salmontools

# Make sure data for Transcriptome Index tests is downloaded.
tx_index_test_raw_dir="$volume_directory/raw/TEST/TRANSCRIPTOME_INDEX"
fasta_file="aegilops_tauschii_short.fa.gz"
if [ ! -e "$tx_index_test_raw_dir/$fasta_file" ]; then
    echo "Downloading fasta file for Transcriptome Index tests."
    wget -q -O "$tx_index_test_raw_dir/$fasta_file" \
         "$test_data_repo/$fasta_file"
fi
gtf_file="aegilops_tauschii_short.gtf.gz"
if [ ! -e "$tx_index_test_raw_dir/$gtf_file" ]; then
    echo "Downloading GTF file for Transcriptome Index tests."
    wget -q -O "$tx_index_test_raw_dir/$gtf_file" \
         "$test_data_repo/$gtf_file"
fi

# Illumina test file
ilu_file="GSE22427_non-normalized.txt"
ilu_test_raw_dir="$volume_directory/raw/TEST/ILLUMINA"
if [ ! -e "$ilu_test_raw_dir/$ilu_file" ]; then
    mkdir -p $ilu_test_raw_dir
    echo "Downloading Illumina file for Illumina tests."
    wget -q -O "$ilu_test_raw_dir/$ilu_file" \
         "$test_data_repo/$ilu_file"
fi

# Agilnt Two Color test file
at_file="GSM466597_95899_agilent.txt"
at_test_raw_dir="$volume_directory/raw/TEST/AGILENT_TWOCOLOR"
if [ ! -e "$at_test_raw_dir/$at_file" ]; then
    mkdir -p $at_test_raw_dir
    echo "Downloading Agilent file for A2C tests."
    wget -q -O "$at_test_raw_dir/$at_file" \
         "$test_data_repo/$at_file"
fi

# Ensure permissions are set for everything within the test data directory.
chmod -R a+rwX $volume_directory

docker build -t dr_worker_tests -f workers/Dockerfile.tests .

source common.sh
HOST_IP=$(get_ip_address)
DB_HOST_IP=$(get_docker_db_ip_address)

docker run \
       --add-host=database:$DB_HOST_IP \
       --add-host=nomad:$HOST_IP \
       --env-file workers/environments/test \
       --volume $volume_directory:/home/user/data_store \
       --link drdb:postgres \
       -it dr_worker_tests bash -c "$(run_tests_with_coverage $@)"
