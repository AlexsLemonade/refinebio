#!/bin/bash -e

# Script for executing Django PyUnit tests within a Docker container.

while getopts ":t:h" opt; do
    case $opt in
        t)
            tag=$OPTARG
            ;;
        h)
            echo "Runs the workers tests. These tests require different Docker containers depending "
            echo "on which code will be tested."
            echo '- by default runs all tests in the workers project, use -t to specify a tag to pass in.'
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            exit 1
            ;;
    esac
done

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

if [[ -z $tag || $tag == "salmon" ]]; then
    # Download salmontools test data
    rm -rf $volume_directory/salmontools/
    git clone git@github.com:dongbohu/salmontools_tests.git $volume_directory/salmontools

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
fi

if [[ -z $tag || $tag == "affymetrix" || $tag == "no_op" ]]; then
    # Make sure CEL for test is downloaded from S3
    cel_name="GSM1426071_CD_colon_active_1.CEL"
    cel_name2="GSM45588.CEL"
    cel_test_raw_dir="$volume_directory/raw/TEST/CEL"
    cel_test_data_1="$cel_test_raw_dir/$cel_name"
    cel_test_data_2="$cel_test_raw_dir/$cel_name2"
    if [ ! -e "$cel_test_data_1" ]; then
        mkdir -p $cel_test_raw_dir
        echo "Downloading CEL for tests."
        wget -q -O $cel_test_data_1 \
             "$test_data_repo/$cel_name"
    fi
    if [ ! -e "$cel_test_data_2" ]; then
        echo "Downloading Non-Brainarray CEL for tests."
        wget -q -O $cel_test_data_2 \
             "$test_data_repo/$cel_name2"
    fi
fi

if [[ -z $tag || $tag == "transcriptome" ]]; then
    # Make sure data for Transcriptome Index tests is downloaded.
    tx_index_test_raw_dir="$volume_directory/raw/TEST/TRANSCRIPTOME_INDEX"
    fasta_file="aegilops_tauschii_short.fa.gz"
    if [ ! -e "$tx_index_test_raw_dir/$fasta_file" ]; then
        mkdir -p $tx_index_test_raw_dir
        echo "Downloading fasta file for Transcriptome Index tests."
        wget -q -O "$tx_index_test_raw_dir/$fasta_file" \
             "$test_data_repo/$fasta_file"
    fi
fi

if [[ -z $tag || $tag == "illumina" ]]; then
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
fi

source common.sh
HOST_IP=$(get_ip_address)
DB_HOST_IP=$(get_docker_db_ip_address)

# Ensure permissions are set for everything within the test data directory.
chmod -R a+rwX $volume_directory

worker_images=(affymetrix illumina salmon transcriptome no_op downloaders)

for image in ${worker_images[*]}; do
    if [[ -z $tag || $tag == $image ]]; then
        image_name=ccdl/dr_$image

        # Strip out tag argument
        tag_string="-t $tag"
        args_without_tag="$(echo $@ | sed "s/-t $tag//")"
        test_command="$(run_tests_with_coverage --tag=$image $args_without_tag)"

        echo "Running tests with the following command:"
        echo $test_command
        docker run \
               --add-host=database:$DB_HOST_IP \
               --add-host=nomad:$HOST_IP \
               --env-file workers/environments/test \
               --volume $volume_directory:/home/user/data_store \
               --link drdb:postgres \
               -it $image_name bash -c "$test_command"
    fi
done
