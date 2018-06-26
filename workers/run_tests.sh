#!/bin/bash -e

# Script for executing Django PyUnit tests within a Docker container.

print_description() {
    echo "Runs the tests for workers. These tests require different Docker containers"
    echo "depending on which code will be tested. By default all tests in the workers"
    echo "project are run."
}

print_options() {
    echo "Options:"
    echo "    -h       Prints the help message"
    echo "    -t TAG   Runs all tests that are tagged with \$TAG"
}

while getopts ":t:h" opt; do
    case $opt in
        t)
            tag=$OPTARG
            ;;
        h)
            print_description
            echo
            print_options
            exit 0
            ;;
        \?)
            echo "Invalid option: -$OPTARG" >&2
            print_options >&2
            exit 1
            ;;
        :)
            echo "Option -$OPTARG requires an argument." >&2
            print_options >&2
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

# Ensure that postgres is running
if ! [[ $(docker ps --filter name=drdb -q) ]]; then
    echo "You must start Postgres first with './run_postgres.sh'" >&2
    exit 1
fi

# Set up the test data volume directory if it does not already exist
volume_directory="$script_directory/test_volume"
if [ ! -d "$volume_directory" ]; then
    mkdir $volume_directory
    chmod -R a+rwX $volume_directory
fi

test_data_repo="https://s3.amazonaws.com/data-refinery-test-assets"

if [[ -z $tag || $tag == "salmon" ]]; then
    # Download "salmon quant" test data
    echo "Downloading 'salmon quant' test data..."
    rm -rf $volume_directory/salmon_tests/

    wget -q -O $volume_directory/salmon_tests.tar.gz $test_data_repo/salmon_tests.tar.gz
    tar xzf $volume_directory/salmon_tests.tar.gz -C $volume_directory

    # Download salmontools test data
    rm -rf $volume_directory/salmontools/
    git clone https://github.com/dongbohu/salmontools_tests.git $volume_directory/salmontools

    # Download tximport test data
    rm -rf $volume_directory/tximport_test/
    git clone https://github.com/dongbohu/tximport_test.git $volume_directory/tximport_test

    # Make sure test Transcriptome Index is downloaded from S3 for salmon tests.
    index_dir="$volume_directory/processed/TEST/TRANSCRIPTOME_INDEX"
    index_tarball="Caenorhabditis_elegans_short_1527089586.tar.gz"
    gz_index_path="$index_dir/$index_tarball"
    if [ ! -e "$gz_index_path" ]; then
        mkdir -p $index_dir
        echo "Downloading Salmon index for Salmon tests."
        wget -q -O $gz_index_path \
             "$test_data_repo/$index_tarball"
    fi

    # Make sure data for Salmon test is downloaded from S3.
    rna_seq_test_raw_dir="$volume_directory/raw/TEST/SALMON"
    read_1_name="ERR1562482_1.fastq.gz"
    read_2_name="ERR1562482_2.fastq.gz"
    rna_seq_test_data_1="$rna_seq_test_raw_dir/$read_1_name"
    rna_seq_test_data_2="$rna_seq_test_raw_dir/$read_2_name"
    if [ ! -e "$rna_seq_test_data_1" ]; then
        mkdir -p $rna_seq_test_raw_dir
        echo "Downloading $read_1_name for Salmon tests."
        wget -q -O $rna_seq_test_data_1 \
             "$test_data_repo/$read_1_name"
        echo "Downloading $read_2_name for Salmon tests."
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
fi

if [[ -z $tag || $tag == "agilent" ]]; then
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

if [[ -z $tag || $tag == "smasher" ]]; then
    # Make sure PCL for test is downloaded from S3
    pcl_name="GSM1237810_T09-1084.PCL"
    pcl_name2="GSM1237812_S97-PURE.PCL"
    pcl_name3="GSM1238108-tbl-1.txt"
    pcl_test_raw_dir="$volume_directory/PCL"
    pcl_test_data_1="$pcl_test_raw_dir/$pcl_name"
    pcl_test_data_2="$pcl_test_raw_dir/$pcl_name2"
    pcl_test_data_3="$pcl_test_raw_dir/$pcl_name3"
    if [ ! -e "$pcl_test_data_1" ]; then
        mkdir -p $pcl_test_raw_dir
        echo "Downloading PCL for tests."
        wget -q -O $pcl_test_data_1 \
             "$test_data_repo/$pcl_name"
    fi
    if [ ! -e "$pcl_test_data_2" ]; then
        echo "Downloading PCL2 for tests."
        wget -q -O $pcl_test_data_2 \
             "$test_data_repo/$pcl_name2"
    fi
    if [ ! -e "$pcl_test_data_3" ]; then
        echo "Downloading PCL3 for tests."
        wget -q -O $pcl_test_data_3 \
             "$test_data_repo/$pcl_name3"
    fi

    export AWS_ACCESS_KEY_ID=`~/bin/aws configure get default.aws_access_key_id`
    export AWS_SECRET_ACCESS_KEY=`~/bin/aws configure get default.aws_secret_access_key`
fi

source common.sh
HOST_IP=$(get_ip_address)
DB_HOST_IP=$(get_docker_db_ip_address)

# Ensure permissions are set for everything within the test data directory.
chmod -R a+rwX $volume_directory

worker_images=(affymetrix illumina salmon transcriptome no_op downloaders agilent smasher)

for image in ${worker_images[*]}; do
    if [[ -z $tag || $tag == $image ]]; then
        if [[ $image == "agilent" || $image == "affymetrix" ]]; then
            # Agilent uses the same docker image as Affymetrix
            ./prepare_image.sh -p -i affymetrix -s workers
            image_name=ccdlstaging/dr_affymetrix
        else
            ./prepare_image.sh -i $image -s workers
            image_name=ccdlstaging/dr_$image
        fi

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
               --env AWS_ACCESS_KEY_ID \
               --env AWS_SECRET_ACCESS_KEY \
               --volume $volume_directory:/home/user/data_store \
               --link drdb:postgres \
               -it $image_name bash -c "$test_command"
    fi
done
