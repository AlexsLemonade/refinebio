
./scripts/update_models.sh

./scripts/prepare_image.sh -i downloaders -d localhost:5000
docker push localhost:5000/dr_downloaders

./scripts/prepare_image.sh -i no_op -d localhost:5000
docker push localhost:5000/dr_no_op