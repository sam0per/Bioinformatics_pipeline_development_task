#!/bin/bash

#$ -cwd

bed=$1
out=$2
var=$3

cat $bed | awk '{for(i=$2;i<$3+1;i++) print $1"\\t"i}' | sort -nk2 | uniq > $out

invar=$(comm -23 $out $var | wc -l)

printf 'TN\tSPECIFICITY\n\n\n20\t10' | paste results/qc/gatkC_samt_truth.tsv -