#!/bin/sh

set -e

# Update the apt keys for dl.google.com
apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 78BD65473CB3BD13
# Update the apt keys for packagecloud.io
curl -L https://packagecloud.io/circleci/trusty/gpgkey | apt-key add -

# Update the local CA certificates
# adapted from: https://stackoverflow.com/questions/62107565/wget-and-curl-stopped-working-with-https-wrongly-complain-about-an-expired-cert
apt update && apt install ca-certificates
sed -i '/mozilla\/AddTrust_External/d' /etc/ca-certificates.conf
update-ca-certificates -f -v
