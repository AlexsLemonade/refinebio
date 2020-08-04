#!/bin/sh

# On the nomad client instances, our instance user data script is responsible
# for a lot of setup. This includes creating a bunch of files on the client
# instance, so we tar them up and append them to the end of the script. This is
# a similar trick to the script that we use to install the Aspera cli.

# This script should always run as if it were being called from
# the directory it lives in.
script_directory="$(perl -e 'use File::Basename;
 use Cwd "abs_path";
 print dirname(abs_path(@ARGV[0]));' -- "$0")"
cd "$script_directory" || exit

if [ -f client-instance-user-data.tpl.sh.bak ]; then
    mv client-instance-user-data.tpl.sh.bak client-instance-user-data.tpl.sh
fi

cp client-instance-user-data.tpl.sh client-instance-user-data.tpl.sh.bak

# Remove all comments and empty lines except for the shebang
sed -i -e '/^#[^!]/d' -e '/^$/d' client-instance-user-data.tpl.sh

# Add this to the end of the file so we can stick things after it.
# The script will exit before #__EOF__ but untar everything after it.
printf "exit 0\n#__EOF__\n" >> client-instance-user-data.tpl.sh

# Tar up all of the files we need in the client instance
tar c nomad-client.service ../nomad-job-specs clean-nomad-jobs.sh | \
    # Gzip them (using --best because we only have 16k to work with)
    gzip --best | \
    # base64 encode the output because terraform requires the file to be valid UTF-8
    # then append to the end of client-instance-user-data.sh
    base64 >> client-instance-user-data.tpl.sh
