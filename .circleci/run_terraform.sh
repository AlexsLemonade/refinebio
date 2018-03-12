# TODO: remove ".s3" suffix"

cd
wget https://releases.hashicorp.com/terraform/0.11.3/terraform_0.11.3_linux_amd64.zip
unzip terraform_0.11.3_linux_amd64.zip
sudo mv terraform /usr/local/bin/

cd ~/refineio/infrastructure

# Download encrypted tfstate files from S3
aws s3 cp s3://refinebio-tfstate/$TFSTATE.enc .
aws s3 cp s3://refinebio-tfstate/$TFSTATE_BAK.enc .

# Decrypt tfstate files
openssl aes-256-cbc -d -in $TFSTATE.enc -out $TFSTATE.s3 -k $OPENSSL_KEY
openssl aes-256-cbc -d -in $TFSTATE_BAK.enc -out $TFSTATE_BAK.s3 -k $OPENSSL_KEY

# New deployment
#terraform init
#terraform plan
#terraform apply -auto-approve

# Encrypt new tfstate files
openssl aes-256-cbc -e -in $TFSTATE.s3 -out $TFSTATE.enc -k $OPENSSL_KEY
openssl aes-256-cbc -e -in $TFSTATE_BAK.s3 -out $TFSTATE_BAK.enc -k $OPENSSL_KEY

# Upload encrypted tfstate files back to S3
aws s3 cp $TFSTATE.enc.s3 s3://refinebio-tfstate/
aws s3 cp $TFSTATE_BAK.enc.s3 s3://refinebio-tfstate/
