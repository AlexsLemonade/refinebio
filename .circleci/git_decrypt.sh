# Unlock encrypted files
cd ~/refinebio/.circleci
openssl aes-256-cbc -md md5 -d -in exported.key.enc -out $KEY_FILENAME -k $OPENSSL_KEY
git-crypt unlock $KEY_FILENAME
rm -f $KEY_FILENAME
