# Import Hashicorps' Key.
curl https://keybase.io/hashicorp/pgp_keys.asc | gpg --import


# Install terraform and nomad
cd
TERRAFORM_VERSION=0.11.8
wget https://releases.hashicorp.com/terraform/$TERRAFORM_VERSION/terraform_${TERRAFORM_VERSION}_linux_amd64.zip
wget https://releases.hashicorp.com/terraform/$TERRAFORM_VERSION/terraform_${TERRAFORM_VERSION}_SHA256SUMS
wget https://releases.hashicorp.com/terraform/$TERRAFORM_VERSION/terraform_${TERRAFORM_VERSION}_SHA256SUMS.sig


# Verify the signature file is untampered.
gpg_ok=$(gpg --verify terraform_${TERRAFORM_VERSION}_SHA256SUMS.sig terraform_${TERRAFORM_VERSION}_SHA256SUMS |& grep Good)
if [[ "$gpg_ok" == "" ]]; then
    echo "Could not verify the signature from HashiCorp Security <security@hashicorp.com>."
    exit 1
fi

# Verify the SHASUM matches the binary.
shasum_ok=$(sha256sum -c terraform_${TERRAFORM_VERSION}_SHA256SUMS |& grep OK)
if [[ "$shasum_ok" == "" ]]; then
    echo "Could not verify the Terraform checksum provided by Hashicorp."
    exit 1
fi

unzip terraform_${TERRAFORM_VERSION}_linux_amd64.zip
sudo mv terraform /usr/local/bin/

sudo apt-get update
sudo apt-get install lxc -y  # Install lxc, which is required by nomad

NOMAD_VERSION=0.8.3
wget https://releases.hashicorp.com/nomad/$NOMAD_VERSION/nomad_${NOMAD_VERSION}_linux_amd64.zip
wget https://releases.hashicorp.com/nomad/$NOMAD_VERSION/nomad_${NOMAD_VERSION}_SHA256SUMS
wget https://releases.hashicorp.com/nomad/$NOMAD_VERSION/nomad_${NOMAD_VERSION}_SHA256SUMS.sig


# Verify the signature file is untampered.
gpg_ok=$(gpg --verify nomad_${NOMAD_VERSION}_SHA256SUMS.sig nomad_${NOMAD_VERSION}_SHA256SUMS |& grep Good)
if [[ "$gpg_ok" == "" ]]; then
    echo "Could not verify the signature from HashiCorp Security <security@hashicorp.com>."
    exit 1
fi

# Verify the SHASUM matches the binary.
shasum_ok=$(sha256sum -c nomad_${NOMAD_VERSION}_SHA256SUMS |& grep OK)
if [[ "$shasum_ok" == "" ]]; then
    echo "Could not verify the Nomad checksum provided by Hashicorp."
    exit 1
fi

unzip nomad_${NOMAD_VERSION}_linux_amd64.zip
sudo mv nomad /usr/local/bin/

cd ~/refinebio/infrastructure

# Circle won't set the branch name for us, so do it ourselves.
source ~/refinebio/common.sh
branch=$(get_master_or_dev)

if [[ $branch == "master" ]]; then
    ENVIRONMENT=prod
    BUCKET_NAME="refinebio-tfstate-deploy-production"
elif [[ $branch == "dev" ]]; then
    ENVIRONMENT=staging
    BUCKET_NAME="refinebio-tfstate-deploy-staging"
else
    echo "Why in the world was run_terraform.sh called from a branch other than dev or master?!?!?"
    exit 1
fi

state_files=$(aws s3 ls s3://$BUCKET_NAME)
if [[ ! -z $state_files ]]; then
   # Download encrypted tfstate files from S3, if they exist
   aws s3 cp s3://$BUCKET_NAME/$TFSTATE.enc .
   aws s3 cp s3://$BUCKET_NAME/$TFSTATE_BAK.enc .

   # Decrypt tfstate files
   openssl aes-256-cbc -d -in $TFSTATE.enc -out $TFSTATE -k $OPENSSL_KEY
   openssl aes-256-cbc -d -in $TFSTATE_BAK.enc -out $TFSTATE_BAK -k $OPENSSL_KEY
fi

# New deployment
./deploy.sh -e $ENVIRONMENT -v $CIRCLE_TAG -u circleci
exit_code=$?

# Encrypt new tfstate files
openssl aes-256-cbc -e -in $TFSTATE -out $TFSTATE.enc -k $OPENSSL_KEY
openssl aes-256-cbc -e -in $TFSTATE_BAK -out $TFSTATE_BAK.enc -k $OPENSSL_KEY

# Upload encrypted tfstate files back to S3
aws s3 cp $TFSTATE.enc s3://$BUCKET_NAME/
aws s3 cp $TFSTATE_BAK.enc s3://$BUCKET_NAME/

exit $exit_code
