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
    git clone https://github.com/dongbohu/salmontools_tests.git $volume_directory/salmontools

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
    pcl_test_raw_dir="$volume_directory/PCL"
    pcl_test_data_1="$pcl_test_raw_dir/$pcl_name"
    pcl_test_data_2="$pcl_test_raw_dir/$pcl_name2"
    if [ ! -e "$pcl_test_data_1" ]; then
        mkdir -p $pcl_test_raw_dir
        echo "Downloading PCL for tests."
        wget -q -O $pcl_test_data_1 \
             "$test_data_repo/$pcl_name"
    fi
    if [ ! -e "$pcl_test_data_2" ]; then
        echo "Downloading PCL for tests."
        wget -q -O $pcl_test_data_2 \
             "$test_data_repo/$pcl_name2"
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
        if [[ $image == "agilent" ]]; then
            # Agilent uses the same docker image as Affymetrix
            ./prepare_image.sh -p -i affymetrix -s workers
            image_name=ccdl/dr_affymetrix
        elif [[ $image == "affymetrix" ]]; then
            ./prepare_image.sh -p -i $image -s workers
            image_name=ccdl/dr_$image
        else
            ./prepare_image.sh -i $image -s workers
            image_name=ccdl/dr_$image
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
