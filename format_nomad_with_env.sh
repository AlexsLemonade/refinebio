#!/bin/bash

##
# This script takes your environment variables and uses them to populate
# Nomad job specifications, as defined in each project.
##

while getopts ":p:e:o:h" opt; do
    case $opt in
    p)
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
        echo '- "local" is the default enviroment, use -e to specify "prod" or "test".'
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

if [[ $project != "workers" && $project != "foreman" && $project != "api" ]]; then
    echo 'Error: must specify project as either "api", "workers", or "foreman" with -p.'
    exit 1
fi

if [[ -z $env ]]; then
    env="local"
fi

# Default docker repo.
# This should work for local and test environments, but we want to
# let this be set outside the script so only set it if it isn't
# already set.
if [[ -z $DOCKERHUB_REPO ]]; then
    export DOCKERHUB_REPO="ccdlstaging"
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
    # Prevent test Nomad job specifications from overwriting
    # existing Nomad job specifications.
    export TEST_POSTFIX="_test"
elif [[ $env == "prod" || $env == "staging" || $env == "dev" ]]; then
    # In production we use EFS as the mount.
    export VOLUME_DIR=/var/efs
else
    export VOLUME_DIR=$script_directory/volume
fi

# We need to specify the database and Nomad hosts for development, but
# not for production because we just point directly at the RDS/Nomad
# instances.
# Conversely, in prod we need AWS credentials and a logging config but
# not in development.
# We do these with multi-line environment variables so that they can
# be formatted into development job specs.
if [ $env != "prod" && $env != "staging" && $env != "dev"]; then
    export EXTRA_HOSTS="
    extra_hosts = [\"database:$DB_HOST_IP\",
               \"nomad:$NOMAD_HOST_IP\"]
"
    export AWS_CREDS=""
    export LOGGING_CONFIG=""
    environment_file="environments/$env"
else
    export EXTRA_HOSTS=""
    export AWS_CREDS="
    AWS_ACCESS_KEY_ID = \"$AWS_ACCESS_KEY_ID_WORKER\"
    AWS_SECRET_ACCESS_KEY = \"$AWS_SECRET_ACCESS_KEY_WORKER\""
    # When deploying prod we write the output of Terraform to a
    # temporary environment file.
    environment_file="$script_directory/infrastructure/prod_env"
fi

# Read all environment variables from the file for the appropriate
# project and environment we want to run.
while read line; do
    # Skip all comments (lines starting with '#')
    is_comment=$(echo $line | grep "^#")
    if [[ -n $line ]] && [[ -z $is_comment ]]; then
        export $line
    fi
done < $environment_file

# There is a current outstanding Nomad issue for the ability to
# template environment variables into the job specifications. Until
# this issue is resolved, we are using perl to accomplish this. The
# issue can be found here:
# https://github.com/hashicorp/nomad/issues/1185

# If output_dir wasn't specified then assume the same folder we're
# getting the templates from.
if [[ -z $output_dir ]]; then
    output_dir=nomad-job-specs
elif [[ ! -d "$output_dir" ]]; then
    mkdir $output_dir
fi

export_log_conf (){
    if [[ $env == 'prod' ]]; then
        export LOGGING_CONFIG="
        logging {
          type = \"awslogs\"
          config {
            awslogs-region = \"$REGION\",
            awslogs-group = \"data-refinery-log-group-$USER-$STAGE\",
            awslogs-stream = \"log-stream-$1-docker-$USER-$STAGE\"
          }
        }"
    else
        export LOGGING_CONFIG=""
    fi
}

# This actually performs the templating using Perl's regex engine.
# Perl magic found here: https://stackoverflow.com/a/2916159/6095378
if [[ $project == "workers" ]]; then
    # Iterate over all the template files in the directory.
    for template in $(ls -1 nomad-job-specs | grep \.tpl); do
        # Strip off the trailing .tpl for once we've formatted it.
        output_file=${template/.tpl/}

        # Downloader logs go to a separate log stream.
        if [ $output_file == "downloader.nomad" ]; then
            export_log_conf "downloader"
        else
            export_log_conf "processor"
        fi

        cat nomad-job-specs/$template \
            | perl -p -e 's/\$\{\{([^}]+)\}\}/defined $ENV{$1} ? $ENV{$1} : $&/eg' \
                   > "$output_dir/$output_file$TEST_POSTFIX" \
                   2> /dev/null
    done
elif [[ $project == "foreman" ]]; then
    export_log_conf "surveyor"
    cat nomad-job-specs/surveyor.nomad.tpl \
        | perl -p -e 's/\$\{\{([^}]+)\}\}/defined $ENV{$1} ? $ENV{$1} : $&/eg' \
               > "$output_dir"/surveyor.nomad"$TEST_POSTFIX" \
               2> /dev/null
elif [[ $project == "api" ]]; then
    export_log_conf "api"
    cat environment.tpl \
        | perl -p -e 's/\$\{\{([^}]+)\}\}/defined $ENV{$1} ? $ENV{$1} : $&/eg' \
               > "$output_dir"/environment"$TEST_POSTFIX" \
               2> /dev/null
fi
