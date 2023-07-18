#!/bin/bash

# Script for running a Django management command to test the worker.

while getopts "hi:" opt; do
    case $opt in
    i)
        IMAGE="$OPTARG"
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
        echo "    ./workers/run_job.sh run_downloader_job --job-name=SRA --job-id=12345"
        echo "    ./workers/run_job.sh -i affymetrix run_processor_job --job-name=AGILENT_TWOCOLOR_TO_PCL --job-id=54321"
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

if [[ -z "$IMAGE" ]]; then
    IMAGE="downloaders"
fi

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(
    cd "$(dirname "$0")" || exit
    pwd
)"
cd "$script_directory" || exit

# However in order to give Docker access to all the code we have to
# move up a level
cd ..

# Agilent uses the same image as affymetrix
if [[ "$IMAGE" == "affymetrix" || "$IMAGE" == "agilent" ]]; then
    ./scripts/prepare_image.sh -i affymetrix
    IMAGE_NAME="$DOCKERHUB_REPO/dr_affymetrix"
else
    ./scripts/prepare_image.sh -i "$IMAGE"
    IMAGE_NAME="$DOCKERHUB_REPO/dr_$IMAGE"
fi

volume_directory="$script_directory/volume"
if [ ! -d "$volume_directory" ]; then
    mkdir "$volume_directory"
fi
chmod -R a+rwX "$volume_directory"

. scripts/common.sh

DB_HOST_IP=$(get_docker_db_ip_address)

docker run \
    --add-host=database:"$DB_HOST_IP" \
    --entrypoint ./manage.py \
    --env AWS_ACCESS_KEY_ID \
    --env AWS_SECRET_ACCESS_KEY \
    --env-file workers/environments/local \
    --interactive \
    --link drdb:postgres \
    --platform linux/amd64 \
    --tty \
    --volume "$volume_directory":/home/user/data_store \
    "$IMAGE_NAME" \
    "${@: -3}" "${@: -2}" "${@: -1}"
