#! /bin/bash

# Connects to each host in the file 'hosts'
# 'hosts' can be populated with 'lists_hosts.sh' which contains usage instructions.

if [ $(cssh -v > /dev/null 2>&1; echo $?) == 0 ]; then
    cssh --options "-i data-refinery-key.pem" --username ubuntu $(<hosts)
elif [ $(csshX -v  > /dev/null 2>&1; echo $?) == 2 ]; then # MacOS, csshX has an exit code of 2 for -v...
    csshX --ssh_args "-i data-refinery-key.pem" --login ubuntu --hosts hosts
else
    echo "You must have either cssh (linux) or csshX (MacOSX) installed."
    echo 'cssh can be installed with: `sudo apt install clusterssh`'
    echo 'csshX can be installed with: `brew install csshX`'
fi
