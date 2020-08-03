#!/bin/bash -e

print_description() {
    echo 'This script can be used to destroy an infrastructure stack that was created with deploy.sh.'
}

print_options() {
    echo 'This script accepts the following arguments: -e, -u, and -h.'
    echo 'Neither -e or -u is optional unless TF_VAR_stage and TF_VAR_user is set.'
    echo '-h prints this help message and exits.'
    echo '-e specifies the environment you would like to destroy.'
    echo '-u specifies the username you used to spin up the stack.'
    echo 'Both arguments are needed to determine which stack to destroy.'
}

while getopts ":e:u:h" opt; do
    case $opt in
    e)
        TF_VAR_stage=$OPTARG
        ;;
    u)
        TF_VAR_user=$OPTARG
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

if [[ $TF_VAR_stage != "dev" && $TF_VAR_stage != "staging" && $TF_VAR_stage != "prod" ]]; then
    echo 'Error: must specify environment as either "dev", "staging", or "prod" with -e or TF_VAR_enviroment.'
    exit 1
fi

if [[ -z $TF_VAR_user ]]; then
    echo 'Error: must specify the username by either providing the -u argument or setting TF_VAR_user.'
    exit 1
fi

# If this file still exists, the previous deploy didn't succeed, so recover
if [ -f nomad-configuration/client-instance-user-data.tpl.sh.bak ]; then
    mv nomad-configuration/client-instance-user-data.tpl.sh.bak nomad-configuration/client-instance-user-data.tpl.sh
fi

terraform destroy

rm nomad-configuration/client-instance-user-data.sh
