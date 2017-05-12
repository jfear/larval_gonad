#!/bin/bash

PROJ=/data/fearjm/Projects/larval_gonad
DATA_DIR=$PROJ/data/rnaseq_samples
SAMPLE=$PROJ/lcdb-wf/config/sampletable.tsv

while IFS=$'\t' read -r -a myArray
do
    if [[ "${myArray[0]}" != "sample" ]]; then
        sname=${myArray[0]}
        fname=${myArray[1]}
        if [ ! -e $DATA_DIR/$sname ]; then
            mkdir $DATA_DIR/$sname
        fi
        ln -s $DATA_DIR/raw/$fname $DATA_DIR/$sname/${sname}_R1.fastq.gz
    fi
done < $SAMPLE
