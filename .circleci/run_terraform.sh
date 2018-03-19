# Install terraform and nomad
cd
wget https://releases.hashicorp.com/terraform/0.11.3/terraform_0.11.3_linux_amd64.zip
unzip terraform_0.11.3_linux_amd64.zip
sudo mv terraform /usr/local/bin/

sudo apt-get update
sudo apt-get install lxc -y  # Install lxc, which is required by nomad
wget https://releases.hashicorp.com/nomad/0.7.1/nomad_0.7.1_linux_amd64-lxc.zip
unzip nomad_0.7.1_linux_amd64-lxc.zip
sudo mv nomad /usr/local/bin/

cd ~/refinebio/.circleci/s3_tfstate
BUCKET_NAME=`terraform output terraform_state_s3_bucket`

cd ~/refinebio/infrastructure
terraform init

# Download encrypted tfstate files from S3
aws s3 cp s3://$BUCKET_NAME/$TFSTATE.enc .
aws s3 cp s3://$BUCKET_NAME/$TFSTATE_BAK.enc .

# Decrypt tfstate files
openssl aes-256-cbc -d -in $TFSTATE.enc -out $TFSTATE -k $OPENSSL_KEY
openssl aes-256-cbc -d -in $TFSTATE_BAK.enc -out $TFSTATE_BAK -k $OPENSSL_KEY

# New deployment
TF_VAR_user=deploy TV_VAR_stage=production ./deploy.sh

# Encrypt new tfstate files
openssl aes-256-cbc -e -in $TFSTATE -out $TFSTATE.enc -k $OPENSSL_KEY
openssl aes-256-cbc -e -in $TFSTATE_BAK -out $TFSTATE_BAK.enc -k $OPENSSL_KEY

# Upload encrypted tfstate files back to S3
aws s3 cp $TFSTATE.enc s3://$BUCKET_NAME/
aws s3 cp $TFSTATE_BAK.enc s3://$BUCKET_NAME/
