#!/bin/bash

# This is a script for combining local data with data downloaded from refine.bio

###
# Define utility functions
###
print_description() {
    echo "Runs the tests for workers. These tests require different Docker containers"
    echo "depending on which code will be tested. By default all tests in the workers"
    echo "project are run."
}

print_options() {
    echo "REQUIRED Options:"
    echo "    -I REFINEBIO_DATA   The path to the zipfile downloaded from refine.bio."
    echo "    -D DATA_DIR         The path to the directory holding the data for analysis."
    echo "    -O ORGANISM         The organism that the samples belong to."
    echo "                        Please input the scientific name in underscore-delimited"
    echo "                        uppercase, i.e. 'GALLUS_GALLUS'."
    echo
    echo "Other Options:"
    echo "    -h                  Prints the help message and exits."
    echo "    -T TRANSFORM        The transform to be applied after combination."
    echo "                        The default option is 'MINMAX'."
    echo "                        The other options are 'NONE', 'STANDARD', and 'ROBUST'."

}

# We get the scientific name in underscore-delimited uppercase for the API, but for the surveyor we
# need it space-delimited and every word capitalized, so this function will transform the name for
# the surveyor
surveyor_fmt() {
    perl -e "\$str = $1; \$str = join ' ', map {ucfirst lc} split '_', \$str; print \$str"
}

# We get the scientific name in underscore-delimited uppercase for the API, but for the surveyor we
# need it underline-delimited and capitalized, so this function will transform the name for the local
# copies of the index
local_organism_fmt() {
    perl -e "\$str = $1; \$str = join '_', map {lc} split '_', \$str; print ucfirst \$str"
}

while getopts ":hO:I:D:T:" opt; do
    case $opt in
        I)
            REFINEBIO_DATA=$OPTARG
            ;;
        D)
            DATA_DIR=$OPTARG
            ;;
        T)
            TRANSFORM=$OPTARG
            ;;
        O)
            ORGANISM=$OPTARG
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

###
# Verify the command-line arguments
###

# Check if the refine.bio data was not specified
if [[ -z $REFINEBIO_DATA ]]; then
    echo "You must specify the path to your refine.bio data." >&2
    print_options >&2
    exit 1
fi

# Check if the data directory was not specified
if [[ -z $DATA_DIR ]]; then
    echo "You must specify a directory holding your private data." >&2
    print_options >&2
    exit 1
fi

# Check if the organism was not specified
if [[ -z $ORGANISM ]]; then
    echo "You must specify the organism that your samples belong to." >&2
    print_options >&2
    exit 1
# Make sure $ORGANISM matches a REGEX for two uppercase words separated by an underscore
elif ! echo $ORGANISM | grep -E "^[A-Z]+_[A-Z]+$" > /dev/null; then
    echo "Invalid form for the organism name. You must specify the scientific name in"
    echo "underscore-delimited uppercase, i.e. 'GALLUS_GALLUS'"
    print_options >&2
    exit 1
fi

# Set the default transform if not specified
if [[ -z $TRANSFORM ]]; then
    TRANSFORM="MINMAX"
# Make sure the specified transform is one of the ones allowed in the backend
elif ! echo $TRANSFORM | grep -E "^(MINMAX|NONE|STANDARD|ROBUST)$" > /dev/null; then
    echo "Invalid option $TRANSFORM for transform."
    print_options >&2
    exit 1
fi


###
# Handle environment variables that can be changed for testing and development.
###

# The URL of the API to query. The default one is the refine.bio production API URL.
if [[ -z $API_URL ]]; then
    API_URL="https://api.refine.bio"
fi

###
# Setup script for execution
###

# This script should always run as if it were being called from
# the directory it lives in.
script_directory=`perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0"`
cd $script_directory

# However in order to give Docker access to all the code we have to
# move up a level
cd ..

# Set up the data volume directory if it does not already exist
volume_directory="$script_directory/volume"
if [ ! -d "$volume_directory" ]; then
    mkdir $volume_directory
    chmod -R a+rwX $volume_directory
fi

# Ensure that Nomad is running first
# Double w in `ps` will cause the columns to never be truncated regardless of environment.
if ! ps auxww | grep nomad | grep -v grep > /dev/null; then
    echo "You must start the nomad first with 'sudo -E ./run_nomad.sh'" >&2
    exit 1
# Then ensure postgres is running
elif ! [[ $(docker ps --filter name=drdb -q) ]]; then
    echo "You must start Postgres first with './run_postgres.sh'" >&2
    exit 1
fi

###
# Query the refine.bio API for the indices for the specified organism.
###

# Long index
INDEX_LONG_URL=`curl "$API_URL/transcriptome_indices?organism=$ORGANISM&length=long" 2> /dev/null\
| jq -r '.s3_url'`

# Short index
INDEX_SHORT_URL=`curl "$API_URL/transcriptome_indices?organism=$ORGANISM&length=short" 2> /dev/null\
| jq -r '.s3_url'`

# Put the indices in the correct places
if [[ $(ls $volume_directory/$ORGANISM\_long.gz 2> /dev/null) && \
    $(ls $volume_directory/$ORGANISM\_short.gz 2> /dev/null) ]]; then
    : # The index already exists from a previous run, so do nothing and move on
# If the s3 urls exist, download the indexes from s3
elif [[ $INDEX_LONG_URL != "null" && $INDEX_SHORT_URL != "null" ]]; then
    echo "Downloading remote index..."
    wget -q -O $volume_directory/$ORGANISM\_long.gz $INDEX_LONG_URL
    wget -q -O $volume_directory/$ORGANISM\_short.gz $INDEX_SHORT_URL
# Otherwise, get the files from the surveyor
else
    # Get the organism name into the format used by the surveyor
    LOCAL_ORGANISM=$(local_organism_fmt $ORGANISM)

    # Check if the surveyor already has an index for the organism, and if so just use that
    if [[ $(ls volume/$LOCAL_ORGANISM/$LOCAL_ORGANISM\_long*.gz 2> /dev/null) && \
    $(ls volume/$LOCAL_ORGANISM/$LOCAL_ORGANISM\_short*.gz 2> /dev/null) ]]; then
        echo "No remote index found, but cached local one exists."
    # Otherwise, just manually download and process the index
    else
        echo "There was no pre-built transcriptome index found."
        echo "Falling back to local trancriptome index build."
        ./foreman/run_surveyor.sh survey_transcriptome Ensembl 1 "$(surveyor_fmt $ORGANISM)" > /dev/null
        echo "Downloading and processing (this may take a while)..."

        # Wait on the processing to be done
        while ! [[ $(ls volume/$LOCAL_ORGANISM/$LOCAL_ORGANISM\_short*.gz 2> /dev/null) ]]; do
            sleep 10
        done
        while ! [[ $(ls volume/$LOCAL_ORGANISM/$LOCAL_ORGANISM\_long*.gz 2> /dev/null) ]]; do
            sleep 10
        done
    fi

    # Move the index from the surveyor directory into $volume_directory

    # Collect the globs of the index files into a list variable
    short_file=(volume/$LOCAL_ORGANISM/$LOCAL_ORGANISM\_short*.gz)
    long_file=(volume/$LOCAL_ORGANISM/$LOCAL_ORGANISM\_long*.gz)

    # Copy the first glob match into the $volume_directory
    cp ${short_file[0]} $volume_directory/$ORGANISM\_short.gz
    cp ${long_file[0]} $volume_directory/$ORGANISM\_long.gz
fi
