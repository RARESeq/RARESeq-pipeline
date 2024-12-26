#!/bin/bash
mv $2 ${2}.fastq
mv $3 ${3}.fastq
python $1 ${2}.fastq ${3}.fastq $4 $5
rm ${2}.fastq ${3}.fastq
