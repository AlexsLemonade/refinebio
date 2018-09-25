# Install git-crypt
cd
git clone https://github.com/AGWA/git-crypt.git
cd git-crypt
make
sudo make install

# Unlock encrypted files
cd ~/refinebio/.circleci
openssl aes-256-cbc -d -in exported.key.enc -out $KEY_FILENAME -k $OPENSSL_KEY
git-crypt unlock $KEY_FILENAME
rm -f $KEY_FILENAME
