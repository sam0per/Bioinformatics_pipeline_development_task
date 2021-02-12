#!/bin/bash

#$ -cwd

# Count number of variant bases considering insertions and deletions using vcf files and bedops
# Example: bash ./scripts/count_var_position.sh ground_truth.vcf.gz variant_call.vcf.gz > variant_base_pos.txt

first=$1
second=$2

if [ -d "tmp" ] 
then
    rm -rf tmp
fi

mkdir tmp

for var in "$@"
do
    # Check file extension
    if [[ $var == *.gz ]]
    then
        gunzip -c $var > ${var%.*}
        var=$(echo ${var%.*})
    fi

    # Transform vcf to bed
    vcf2bed --insertions < $var | cut -f 1,2,3 > tmp/$(basename $var).insertions.bed
    vcf2bed --deletions < $var | cut -f 1,2,3 > tmp/$(basename $var).deletions.bed
    bedops --everything tmp/$(basename $var).insertions.bed \
    tmp/$(basename $var).deletions.bed > tmp/$(basename $var).indels.bed

    # # Expand base positions
    cat tmp/$(basename $var).indels.bed | awk '{for(i=$2;i<$3+1;i++) print $1"\t"i}' | \
    sort -nk2 | uniq > tmp/$(basename $var).expandedpos.txt
    vcf2bed --snvs < $var | cut -f 1,3 >> tmp/$(basename $var).expandedpos.txt
done

cat tmp/$(basename ${first%.*})*.expandedpos.txt tmp/$(basename ${second%.*})*.expandedpos.txt | sort -nk2 | uniq
