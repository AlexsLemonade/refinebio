#!/bin/sh

# Exit on failure
set -e

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

# Get access to all of refinebio
cd ..


print_description() {
    echo 'This script can download dev data for developement,'
    echo "it creates or replaces: \'./volume\' and \'volumes_postgres\'"
}

print_options() {
    echo '-h for help'
}

while getopts ":h" opt; do
    case $opt in
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

VOLUME_FOLDER=./volume
POSTGRES_VOLUME_FOLDER=./volumes_postgres

download_dev_data() {

    # Download Archived Data

    curl -o dev-data.tar.gz https://refinebio-dev-data.s3.amazonaws.com/dev-data.tar.gz

    # remove existing data folder
    if [ -f "$VOLUME_FOLDER" ]; then
        echo "found existing volumes folder... removing before replacing"
        rm -r "$VOLUME_FOLDER"
    fi

    # remove existing db folder
    if [ -f "$POSTGRES_VOLUME_FOLDER" ]; then
        echo "Found existing volumes_postgres... removing before replacing"
        rm -r "$POSTGRES_VOLUME_FOLDER"
    fi
    # Unarchive Volume and Data
    echo "unarchiving dev-data.tar.gz"
    tar -xzf dev-data.tar.gz

    # Clean Up
    echo "cleaning up downloaded zip file"
    rm ./dev-data.tar.gz

}


# donwload_dev_data
echo "Done!"
exit 0
