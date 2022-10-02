#!/bin/bash

list_descendants ()
{
  local children=$(ps -o pid= --ppid "$1")

  for pid in $children
  do
    list_descendants "$pid"
  done

  echo "$children"
}

if [  $# -ne 2 ]; then
  echo "Use: ./unload_index.sh sample_name number_of_indicies"
  exit 1
fi

let TO=$2-1

for NO in  $(seq 0 ${TO});
do

    NAME=$PHASEDANCER_DATA_DIR/data/$1/index/$1_${NO}

    PID_FILE=${NAME}.index_pids.txt

    while read -r PID; do

        COMMAND=$(ps -o cmd= --pid $PID)

        if [ -z "$COMMAND" ]; then
            echo "No process with $PID."
            exit 1
        fi

        if [[ "$COMMAND" != *load_index* ]]; then
            echo "Error in process" $DESC_PID ": Root process should be executing a bash command."
            exit 1
        fi

        DESC_PIDS=$(list_descendants $PID);

        for DESC_PID in $DESC_PIDS; do

            COMMAND=$(ps -o cmd= --pid $DESC_PID)
            COMMAND_OK=$(echo $COMMAND | grep -c -P "($NAME|bash)")

            if [ $COMMAND_OK -ne 1 ]; then
                ps -o cmd= --pid $DESC_PID
                echo "Error in process" $DESC_PID ": Command (" $COMMAND ") is not related to" ${NAME}.mmi "index."
                exit 1
            fi

        done

    done  < $PID_FILE

    while read -r PID; do

        DESC_PIDS=$(list_descendants $PID);
        kill -9 $PID $DESC_PIDS;
        echo "Processes with PIDs:" $PID $DESC_PIDS "killed.";

    done  < $PID_FILE
done
