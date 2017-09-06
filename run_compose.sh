# This script should always run as if it were being called from
# the directory it lives in.
script_directory=`dirname "${BASH_SOURCE[0]}"  | xargs realpath`
cd $script_directory

export VOLUME_DIRECTORY="$script_directory/workers/volume"
export STATIC_DIRECTORY="$script_directory/end_to_end_tests"
export HOST_IP=$(ip route get 8.8.8.8 | awk '{print $NF; exit}')

export EXPERIMENT_LIST=single-assay-experiment.txt

docker-compose build

docker-compose up
