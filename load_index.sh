#!/bin/bash

if [  $# -ne 2 ]; then
  echo "Use: ./load_index.sh sample_name number_of_indicies"
  exit 1
fi

let TO=$2-1

for NO in  $(seq 0 ${TO});
do

    NAME=$PHASEDANCER_DATA_DIR/data/$1/index/$1_${NO}

    if [ ! -p ${NAME}_in.fifo ]; then

        mkfifo ${NAME}_in.fifo

        if [ $? -ne 0 ]; then
            exit 1
        fi

    fi

    if [ ! -p ${NAME}_out.fifo ]; then

        mkfifo ${NAME}_out.fifo

        if [ $? -ne 0 ]; then
            exit 1
        fi

    fi

    if [ ! -f ${NAME}.mmi ]; then

        echo "File: \"${NAME}.mmi\" does not exists."
        exit 1
    fi

    if [ ! -f ${NAME}.fasta ]; then
       echo "File: \"${NAME}.fasta\" does not exists."
       exit 1
    fi

    which unbuffer  >> /dev/null 2>&1

    if [ $? -ne 0 ]; then
        echo "unbuffer is not installed."
        exit 1
    fi

    PID_FILE=${NAME}.index_pids.txt

    { while true; do unbuffer cat ${NAME}_in.fifo; done > ${NAME}_out.fifo; } &

    PID=$!

    echo "PID of fifo queue parent process: " $PID
    echo $PID > $PID_FILE
    echo $NAME.fasta


    { unbuffer  minimap2 -p 0 -N 300000 -x map-pb -t 1 -K 1 --sam-hit-only ${NAME}.mmi <(unbuffer cat  ${NAME}_out.fifo) | unbuffer -p  python3 src/mem_index/filter.py "${NAME}.fasta" $NO  > output_file; } &

    PID=$!
    echo "PID of minimap2 index parent process: " $PID
    echo $PID >> $PID_FILE
done
