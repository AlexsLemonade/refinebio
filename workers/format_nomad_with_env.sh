#!/bin/bash

while getopts ":e:o:h" opt; do
    case $opt in
        e)
            env=$OPTARG
            ;;
        o)
            output_dir=$OPTARG
            ;;
        h)
            echo "Formats Nomad Job Specifications with the specified environment overlaid "
            echo "onto the current environment."
            echo '- "dev" is the default enviroment, use -e to specify "prod" or "test".'
            echo '- "workers/" is the default output directory, use -o to specify'
            echo '      an absolute path to a directory (trailing / must be included).'
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

if [[ -z $env ]]; then
    env="dev"
fi

# This script should always run as if it were being called from
# the directory it lives in.
script_directory=`perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0"`
cd $script_directory

# It's important that these are run first so they will be overwritten
# by environment variables.
source ../common.sh
export DB_HOST_IP=$(get_docker_db_ip_address)
export NOMAD_HOST_IP=$(get_ip_address)

# What to do for the "prod" env is TBD.
if [ $env == "test" ]; then
    export VOLUME_DIR=$script_directory/test_volume
elif [ $env == "prod" ]; then
    export VOLUME_DIR=/home/ubuntu/docker_volume
else
    export VOLUME_DIR=$script_directory/volume
fi

while read line; do
    is_comment=$(echo $line | grep "^#")
    if [[ -n $line ]] && [[ -z $is_comment ]]; then
        export $line
    fi
done < "environments/$env"

# There is a current outstanding Nomad issue for the ability to
# template environment variables into the job specifications. Until
# this issue is resolved, we are using perl to accomplish this. The
# issue can be found here:
# https://github.com/hashicorp/nomad/issues/1185

if [[ ! -z $env && ! -d "$output_dir" ]]; then
    mkdir $output_dir
fi

# Perl magic found here: https://stackoverflow.com/a/2916159/6095378
cat downloader.nomad.tpl \
    | perl -p -e 's/\$\{\{([^}]+)\}\}/defined $ENV{$1} ? $ENV{$1} : $&/eg' \
           > "$output_dir"downloader.nomad \
           2> /dev/null

cat processor.nomad.tpl \
    | perl -p -e 's/\$\{\{([^}]+)\}\}/defined $ENV{$1} ? $ENV{$1} : $&/eg' \
           > "$output_dir"processor.nomad \
           2> /dev/null
