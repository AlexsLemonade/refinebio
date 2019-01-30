#! /bin/bash

docker rm -f dres 2> /dev/null

docker run --name dres -p 9200:9200 -p 9300:9300 -e "discovery.type=single-node" -d docker.elastic.co/elasticsearch/elasticsearch:6.5.4
