#!/bin/bash

wget https://mimuw.edu.pl/~bp209493/phasedancer_example.tar.gz
tar -xvzf phasedancer_example.tar.gz
docker run -d -v `pwd`/phasedancer_example:/phaseDancerData -it --name phasedancer phasedancer:1.0
docker exec phasedancer ./run_phaseDancer.sh flat2 all 1 10
docker exec phasedancer ./unload_index.sh flat2 1
docker stop phasedancer
docker rm phasedancer
