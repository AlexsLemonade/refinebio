#!/bin/bash

docker run -d --hostname rabbit-queue --name message-queue rabbitmq:3
