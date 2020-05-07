#!/bin/bash

if [ "$(docker ps -aq -f name=minimization)" != "" ]; then
    echo "removing exited container"
    docker rm -f minimization
fi

docker run -d \
--name minimization \
--restart unless-stopped \
-e ARGS="$*" \
minimization
