#! /bin/sh

docker rm -f dres 2> /dev/null

# Check if a docker database named "dres" exists, and if so just start it
if [ "$(docker ps -a --filter name=dres -q)" ]; then
  docker start dres > /dev/null
# Otherwise, run it with `docker run`
else
    docker run --name dres -p 9200:9200 -p 9300:9300 -e "discovery.type=single-node" -e "indices.query.bool.max_clause_count=16384" -d docker.elastic.co/elasticsearch/elasticsearch:6.5.4
fi

echo "Started ElasticSearch."
