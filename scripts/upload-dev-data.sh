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
    echo 'This script can save a local DB and files to s3.'
    echo 'You must login to aws-cli with "$ aws login"'
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

upload_data() {

    # Archive data and volume
    echo "Archiving data..."
    tar -zcvf dev-data.tar.gz volume volumes_postgres

    # Send to s3
    echo "Sending to S3..."
    s3 cp dev-data.tar.gz s3://refinebio-dev-data/dev-data.tar.gz

    # Clean up
    echo "Cleaning up..."
    rm dev-data.tar.gz

}

if ! [ -x "$(command -v aws)" ]; then
    echo 'Error: aws is not installed.' >&2
    exit 1
fi

upload_data

echo "Done!"
exit 0
