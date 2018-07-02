#!/bin/bash

print_description() {
    echo "Runs the tests for workers. These tests require different Docker containers"
    echo "depending on which code will be tested. By default all tests in the workers"
    echo "project are run."
}

print_options() {
    echo "Options:"
    echo "    -h       Prints the help message"
    echo "    -I REFINEBIO_DATA   The path to the zipfile downloaded from refine.bio."
    echo "    -D DATA_DIR         The path to the directory holding the data for analysis."
    echo "    -T TRANSFORM        The transform to be applied after combination."
    echo "                        The default option is 'MINMAX'."
    echo "                        The other options are 'NONE', 'STANDARD', and 'ROBUST'."
    echo "    -O ORGANISM         The organism that the samples belong to."
    echo "                        Please input the scientific name in underscore-delimited"
    echo "                        uppercase, i.e. 'GALLUS_GALLUS'."
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

# Check if not specified
if [[ -z $REFINEBIO_DATA ]]; then
    echo "You must specify the path to your refine.bio data." >&2
    print_options >&2
    exit 1
fi

# Check if not specified
if [[ -z $DATA_DIR ]]; then
    echo "You must specify a directory holding your private data." >&2
    print_options >&2
    exit 1
fi

# Check if not specified
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