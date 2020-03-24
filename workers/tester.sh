#!/bin/bash

# Script for running a django management command to test the worker.

while getopts "hi:" opt; do
    case $opt in
        i)
            image=$OPTARG
            ;;
        h)
            echo "Runs a downloader or processor job. The following arguments are supported:"
            echo "-h : Print this help message and exit."
            echo "-i <IMAGE_NAME> : The image to use. Options are:"
            echo "    downloaders (default)"
            echo "    salmon"
            echo "    transcriptome"
            echo "    no_op"
            echo "    downloaders"
            echo "    illumina"
            echo "    affymetrix"
            echo "<MANAGEMENT COMMAND> : What kind of job to run."
            echo "    Must be either 'run_downloader_job' or 'run_processor_job'."
            echo "--job-name=<JOB_NAME> : The type of job to run."
            echo "    For processor jobs, options are:"
            echo "        AFFY_TO_PCL"
            echo "        AGILENT_TWOCOLOR_TO_PCL"
            echo "        SALMON"
            echo "        ILLUMINA_TO_PCL"
            echo "        TRANSCRIPTOME_INDEX_LONG"
            echo "        TRANSCRIPTOME_INDEX_SHORT"
            echo "        NO_OP"
            echo "    For downloader jobs, options are:"
            echo "        ARRAY_EXPRESS"
            echo "        SRA"
            echo "        TRANSCRIPTOME_INDEX"
            echo "        GEO"
            echo "--job-id=<JOB_ID> : The id of the job you want to run. Must already exist in the database."
            echo ""
            echo "Note that the <IMAGE_NAME> must correspond to the <JOB_NAME>."
            echo "     AGILENT_TWOCOLOR_TO_PCL is a special case because it requires the 'affymetrix' image."
            echo ""
            echo "Examples:"
            echo "    ./workers/tester.sh run_downloader_job --job-name=SRA --job-id=12345"
            echo "    ./workers/tester.sh -i affymetrix run_processor_job --job-name=AGILENT_TWOCOLOR_TO_PCL --job-id=54321"
            exit 0
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

if [[ -z "$image" ]]; then
    image="downloaders"
fi

# This script should always run as if it were being called from
# the directory it lives in.
script_directory=`perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0"`
cd "$script_directory"

# However in order to give Docker access to all the code we have to
# move up a level
cd ..

# Agilent uses the same image as affymetrix
if [[ "$image" == "affymetrix" || "$image" == "agilent" ]]; then
    ./scripts/prepare_image.sh -p -i affymetrix
    image_name="ccdlstaging/dr_affymetrix"
else
    ./scripts/prepare_image.sh -i "$image"
    image_name="ccdlstaging/dr_$image"
fi

volume_directory="$script_directory/volume"

source scripts/common.sh
HOST_IP=$(get_ip_address)
DB_HOST_IP=$(get_docker_db_ip_address)

export AWS_ACCESS_KEY_ID=`~/bin/aws configure get default.aws_access_key_id`
export AWS_SECRET_ACCESS_KEY=`~/bin/aws configure get default.aws_secret_access_key`

docker run \
	   -it \
       --add-host=database:"$DB_HOST_IP" \
       --add-host=nomad:"$HOST_IP" \
       --env-file workers/environments/local \
       --env AWS_ACCESS_KEY_ID \
       --env AWS_SECRET_ACCESS_KEY \
       --entrypoint ./manage.py \
       --volume "$volume_directory":/home/user/data_store \
       --link drdb:postgres \
       "$image_name" "${@: -3}" "${@: -2}" "${@: -1}"
