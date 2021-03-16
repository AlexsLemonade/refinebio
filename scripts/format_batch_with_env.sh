#!/bin/sh

##
# This script takes your environment variables and uses them to populate
# Nomad job specifications, as defined in each project.
##

print_description() {
    echo "Formats Nomad Job Specifications with the specified environment overlaid "
    echo "onto the current environment."
}

print_options() {
    echo "Options:"
    echo "    -h               Prints the help message"
    echo "    -p PROJECT       The project to format."
    echo "                     Valid values are 'api', 'workers', 'surveyor', or 'foreman'."
    echo "    -e ENVIRONMENT   The environemnt to run. 'local' is the default."
    echo "                     Other valid values are 'prod' or 'staging'."
    echo "    -o OUTPUT_DIR    The output directory. The default directory is the"
    echo "                     project directory, but you can specify an absolute"
    echo "                     path to another directory (trailing / must be included)."
}

while getopts ":p:e:o:v:h" opt; do
    case $opt in
        p)
            export project=$OPTARG
            ;;
        e)
            export env=$OPTARG
            ;;
        o)
            export output_dir=$OPTARG
            ;;
        v)
            export system_version=$OPTARG
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

if [ "$project" != "workers" ] && [ "$project" != "surveyor" ] && [ "$project" != "foreman" ] && [ "$project" != "api" ]; then
    echo 'Error: must specify project as either "api", "workers", "surveyor", or "foreman" with -p.'
    exit 1
fi

if [ -z "$env" ]; then
    env="local"
fi

if [ -z "$system_version" ]; then
    system_version="latest"
fi

# Default docker repo.
# This are very sensible defaults, but we want to let these be set
# outside the script so only set them if they aren't already set.
if [ -z "$DOCKERHUB_REPO" ]; then
    export DOCKERHUB_REPO="ccdlstaging"
fi
if [ -z "$FOREMAN_DOCKER_IMAGE" ]; then
    export FOREMAN_DOCKER_IMAGE="dr_foreman:$system_version"
fi
if [ -z "$DOWNLOADERS_DOCKER_IMAGE" ]; then
    export DOWNLOADERS_DOCKER_IMAGE="dr_downloaders:$system_version"
fi
if [ -z "$TRANSCRIPTOME_DOCKER_IMAGE" ]; then
    export TRANSCRIPTOME_DOCKER_IMAGE="dr_transcriptome:$system_version"
fi
if [ -z "$SALMON_DOCKER_IMAGE" ]; then
    export SALMON_DOCKER_IMAGE="dr_salmon:$system_version"
fi
if [ -z "$SMASHER_DOCKER_IMAGE" ]; then
    export SMASHER_DOCKER_IMAGE="dr_smasher:$system_version"
fi
if [ -z "$AFFYMETRIX_DOCKER_IMAGE" ]; then
    export AFFYMETRIX_DOCKER_IMAGE="dr_affymetrix:$system_version"
fi
if [ -z "$ILLUMINA_DOCKER_IMAGE" ]; then
    export ILLUMINA_DOCKER_IMAGE="dr_illumina:$system_version"
fi
if [ -z "$NO_OP_DOCKER_IMAGE" ]; then
    export NO_OP_DOCKER_IMAGE="dr_no_op:$system_version"
fi
if [ -z "$COMPENDIA_DOCKER_IMAGE" ]; then
    export COMPENDIA_DOCKER_IMAGE="dr_compendia:$system_version"
fi


# This script should always run from the context of the directory of
# the project it is building.
script_directory="$(perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0")"

project_directory="$script_directory/.."

# Correct for foreman and surveyor being in same directory:
if [ "$project" = "surveyor" ]; then
    cd "$project_directory/foreman" || exit
else
    cd "$project_directory/$project" || exit
fi

# It's important that these are run first so they will be overwritten
# by environment variables.

# shellcheck disable=SC1091
. ../scripts/common.sh
DB_HOST_IP=$(get_docker_db_ip_address)
export DB_HOST_IP
NOMAD_HOST_IP=$(get_ip_address)
export NOMAD_HOST_IP

export VOLUME_DIR=/var/ebs

export EXTRA_HOSTS=""
export AWS_CREDS="{\"name\": \"AWS_ACCESS_KEY_ID\", \"value\": \"$AWS_ACCESS_KEY_ID\"},
            {\"name\": \"AWS_SECRET_ACCESS_KEY\", \"value\": \"$AWS_SECRET_ACCESS_KEY\"},"
# When deploying prod we write the output of Terraform to a
# temporary environment file.
environment_file="$project_directory/infrastructure/prod_env"

# Temporarily set the logging config to nothing so we can see the logs
# through nomad.
if [ "$env" = "staging" ]; then
    export LOGGING_CONFIG=""
fi

# Read all environment variables from the file for the appropriate
# project and environment we want to run.
while read -r line; do
    # Skip all comments (lines starting with '#')
    is_comment=$(echo "$line" | grep "^#")
    if [ -n "$line" ] && [ -z "$is_comment" ]; then
        # shellcheck disable=SC2163
        export "$line"
    fi
done < "$environment_file"

# There is a current outstanding Nomad issue for the ability to
# template environment variables into the job specifications. Until
# this issue is resolved, we are using perl to accomplish this. The
# issue can be found here:
# https://github.com/hashicorp/nomad/issues/1185

# If output_dir wasn't specified then assume the same folder we're
# getting the templates from.
if [ -z "$output_dir" ]; then
    output_dir=nomad-job-specs
fi

if [ ! -d "$output_dir" ]; then
    mkdir "$output_dir"
fi

export_log_conf (){
    if [ "$env" = 'prod' ] || [ "$env" = 'staging' ] || [ "$env" = 'dev' ]; then
        export LOGGING_CONFIG="
        logging {
          type = \"awslogs\"
          config {
            awslogs-region = \"$AWS_REGION\",
            awslogs-group = \"data-refinery-log-group-$USER-$STAGE\",
            awslogs-stream = \"log-stream-$1-docker-$USER-$STAGE\"
          }
        }"

        # Only constrain smasher jobs in the cloud so that
        # local/test can still run smasher jobs.
        export SMASHER_CONSTRAINT="
        constraint {
          attribute = \"\${meta.is_smasher}\"
          operator = \"=\"
          value = \"true\"
        }"
    else
        export LOGGING_CONFIG=""
        export SMASHER_CONSTRAINT=""
    fi
}

# Not quite sure how to deal with this just yet, so punt.
# export INDEX=0

# This actually performs the templating using Perl's regex engine.
# Perl magic found here: https://stackoverflow.com/a/2916159/6095378
if [ "$project" = "workers" ]; then
    # Iterate over all the template files in the directory.
    for template in batch-job-templates/*.tpl.json; do
	template="$(basename "$template")"
        # Strip off the trailing .tpl for once we've formatted it.
        OUTPUT_BASE="$(basename "$template" .tpl.json)"
        FILETYPE=".json"
        OUTPUT_FILE="$OUTPUT_BASE$FILETYPE"

        # Downloader logs go to a separate log stream.
        if [ "$OUTPUT_FILE" = "downloader.json" ]; then
            export_log_conf "downloader"
            rams="1024 4096 16384"
            for r in $rams
            do
                export RAM_POSTFIX="_$r"
                export RAM="$r"
                FILEPATH="$output_dir/$OUTPUT_FILE$RAM_POSTFIX$FILETYPE"
                perl -p -e 's/\$\{\{([^}]+)\}\}/defined $ENV{$1} ? $ENV{$1} : $&/eg' \
                     < "batch-job-templates/$template" \
                     > "$FILEPATH" \
                     2> /dev/null
                echo "Made $FILEPATH"
            done
            echo "Made $output_dir/$OUTPUT_FILE"
        elif [ "$OUTPUT_FILE" = "smasher.json" ] || [ "$OUTPUT_FILE" = "create_qn_target.json" ] || [ "$OUTPUT_FILE" = "create_compendia.json" ] || [ "$OUTPUT_FILE" = "create_quantpendia.json" ] || [ "$OUTPUT_FILE" = "tximport.json" ]; then
            export_log_conf "processor"
            perl -p -e 's/\$\{\{([^}]+)\}\}/defined $ENV{$1} ? $ENV{$1} : $&/eg' \
                 < "batch-job-templates/$template" \
                 > "$output_dir/$OUTPUT_FILE" \
                 2> /dev/null
            echo "Made $output_dir/$OUTPUT_FILE"
        else
            export_log_conf "processor"
            # rams="1024 2048 3072 4096 5120 6144 7168 8192 9216 10240 11264 12288 13312"
            rams="2048 3072 4096 8192 12288 16384 32768 65536"
            for r in $rams
            do
                export RAM_POSTFIX="_$r"
                export RAM="$r"
                FILEPATH="$output_dir/$OUTPUT_BASE$RAM_POSTFIX$FILETYPE"
                perl -p -e 's/\$\{\{([^}]+)\}\}/defined $ENV{$1} ? $ENV{$1} : $&/eg' \
                     < "batch-job-templates/$template" \
                     > "$FILEPATH" \
                     2> /dev/null
                echo "Made $FILEPATH"
            done
        fi
    done
elif [ "$project" = "surveyor" ]; then

    # Iterate over all the template files in the directory.
    for template in batch-job-templates/*.tpl.json; do
	template="$(basename "$template")"
        # Strip off the trailing .tpl for once we've formatted it.
        OUTPUT_BASE="$(basename "$template" .tpl.json)"
        FILETYPE=".json"
        OUTPUT_FILE="$OUTPUT_BASE$FILETYPE"

        # Downloader logs go to a separate log stream.
        if [ "$OUTPUT_FILE" = "surveyor_dispatcher.json" ]; then
            export_log_conf "surveyor_dispatcher"
            perl -p -e 's/\$\{\{([^}]+)\}\}/defined $ENV{$1} ? $ENV{$1} : $&/eg' \
                 < "batch-job-templates/$template" \
                 > "$output_dir/$OUTPUT_FILE" \
                 2> /dev/null
            echo "Made $output_dir/$OUTPUT_FILE"
        else
            export_log_conf "surveyor"
            rams="256 4096 16384"
            for r in $rams
            do
                export RAM_POSTFIX="_$r"
                export RAM="$r"
                FILEPATH="$output_dir/$OUTPUT_BASE$RAM_POSTFIX$FILETYPE"
                perl -p -e 's/\$\{\{([^}]+)\}\}/defined $ENV{$1} ? $ENV{$1} : $&/eg' \
                     < "batch-job-templates/$template" \
                     > "$FILEPATH" \
                     2> /dev/null
                echo "Made $FILEPATH"
            done
        fi
    done

elif [ "$project" = "foreman" ]; then
    # foreman sub-project
    export_log_conf "foreman"
    perl -p -e 's/\$\{\{([^}]+)\}\}/defined $ENV{$1} ? $ENV{$1} : $&/eg' \
         < environment.tpl \
         > "$output_dir/environment" \
         2> /dev/null
elif [ "$project" = "api" ]; then
    export_log_conf "api"
    perl -p -e 's/\$\{\{([^}]+)\}\}/defined $ENV{$1} ? $ENV{$1} : $&/eg' \
         < environment.tpl \
         > "$output_dir/environment" \
         2> /dev/null
fi
