#!/bin/bash

curl https://mimuw.edu.pl/~bp209493/phasedancer_example.tar.gz -o phasedancer_example.tar.gz
tar -xvzf phasedancer_example.tar.gz
docker run -d -v `pwd`/phasedancer_example:/phaseDancerData -it --name phasedancer bposzewiecka/phasedancer:1.0
docker exec phasedancer ./run_phaseDancer.sh flat2 all 1 20 >  /dev/null
docker exec phasedancer ./unload_index.sh flat2 1
docker stop phasedancer
docker rm phasedancer
docker run -d -it -v `pwd`/phasedancer_example:/phaseDancerViewer_data -p 8000:8000 --name phasedancer_viewer bposzewiecka/phasedancer_viewer:1.0
