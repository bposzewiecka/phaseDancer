#!/bin/bash

snakemake --snakefile=GenerateIndices.smk --config sample=$1 --cores 5
./load_index.sh $1 $3
snakemake --config sample=$1 contig=$2 --cores $4
