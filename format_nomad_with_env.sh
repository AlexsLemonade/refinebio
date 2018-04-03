#!/bin/bash

##
# This script takes your environment variables and uses them to populate
# Nomad job specifications, as defined in each project. 
##

while getopts ":p:e:o:h" opt; do
    case $opt in
        p)
            if [[ $OPTARG != "workers" && $OPTARG != "foreman" && $OPTARG != "api" ]]; then
                echo 'Error: -p must specify either "api", workers" or "foreman".'
                exit 1
            fi
            project=$OPTARG
            ;;
        e)
            env=$OPTARG
            ;;
        o)
            output_dir=$OPTARG
            ;;
        h)
            echo "Formats Nomad Job Specifications with the specified environment overlaid "
            echo "onto the current environment."
            echo '-p specifies the project to format. Valid values are "api", workers" or "foreman".'
            echo '- "dev" is the default enviroment, use -e to specify "prod" or "test".'
            echo '- the project directory will be used as the default output directory, use -o to specify'
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

# This script should always run from the context of the directory of
# the project it is building.
script_directory=`perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0"`
cd $script_directory/$project

# It's important that these are run first so they will be overwritten
# by environment variables.
source ../common.sh
export DB_HOST_IP=$(get_docker_db_ip_address)
export NOMAD_HOST_IP=$(get_ip_address)

if [ $env == "test" ]; then
    export VOLUME_DIR=$script_directory/test_volume
elif [ $env == "prod" ]; then
    # In production we use EFS as the mount.
    export VOLUME_DIR=/var/efs
else
    export VOLUME_DIR=$script_directory/volume
fi

# We need to specify the database and Nomad hosts for development, but
# not for production because we just point directly at the RDS/Nomad
# instances.
if [ $env != "prod" ]; then
    # This is a multi-line env var so that it can be formatted into
    # development job specs.
    export EXTRA_HOSTS="
        extra_hosts = [\"database:$DB_HOST_IP\",
                       \"nomad:$NOMAD_HOST_IP\"]
"
else
    export EXTRA_HOSTS=""
fi

# Skip all comments (lines starting with '#')
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

if [[ ! -z $output_dir && ! -d "$output_dir" ]]; then
    mkdir $output_dir
fi

# This actually performs the templating using Perl's regex engine.
# Perl magic found here: https://stackoverflow.com/a/2916159/6095378
if [[ $project == "workers" ]]; then
    cat downloader.nomad.tpl \
        | perl -p -e 's/\$\{\{([^}]+)\}\}/defined $ENV{$1} ? $ENV{$1} : $&/eg' \
               > "$output_dir"downloader.nomad \
               2> /dev/null

    cat processor.nomad.tpl \
        | perl -p -e 's/\$\{\{([^}]+)\}\}/defined $ENV{$1} ? $ENV{$1} : $&/eg' \
               > "$output_dir"processor.nomad \
               2> /dev/null
elif [[ $project == "foreman" ]]; then
    cat surveyor.nomad.tpl \
        | perl -p -e 's/\$\{\{([^}]+)\}\}/defined $ENV{$1} ? $ENV{$1} : $&/eg' \
               > "$output_dir"surveyor.nomad \
               2> /dev/null
elif [[ $project == "api" ]]; then
    cat environment.tpl \
        | perl -p -e 's/\$\{\{([^}]+)\}\}/defined $ENV{$1} ? $ENV{$1} : $&/eg' \
               > "$output_dir"environment \
               2> /dev/null
fi
