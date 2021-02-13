#!/bin/bash

#$ -cwd

# Calculate analytical performance statistics using bcftools isec.
# The script takes four arguments:
# The first argument is the raw vcf including also invariant sites
# The second argument is the true variant call
# The third argument is the path to the directory where the results from bcftools isec will be saved
# The fourth argument is the summary output from gatk Concordance
# Example: bash ./scripts/edit_concordance.sh raw_call.vcf.gz ground_truth.vcf.gz isec_dir concordance_summary.tsv

raw=$1
truth=$2
out=$3
tbl=$4

# pipeline: indel invariant sites
gunzip -c $raw | grep -v "INDEL" | envs/var_call_v1/bin/bcftools \
filter -i 'ALT="."' -Oz -o tmp/$(basename ${raw%.*}).indel.neg.gz
tabix -p vcf tmp/$(basename ${raw%.*}).indel.neg.gz

# pipeline: snp invariant sites
vcftools --gzvcf $raw --remove-indels --recode --recode-INFO-all --stdout | envs/var_call_v1/bin/bcftools \
filter -i 'ALT="."' -Oz -o tmp/$(basename ${raw%.*}).snp.neg.gz
tabix -p vcf tmp/$(basename ${raw%.*}).snp.neg.gz

# truth: snp variant sites
vcftools --gzvcf $truth --remove-indels --recode --recode-INFO-all \
--stdout | bgzip -c > tmp/$(basename ${truth%.*}).snp.gz
tabix -p vcf tmp/$(basename ${truth%.*}).snp.gz

# truth: truth variant sites
vcftools --gzvcf $truth --keep-only-indels --recode --recode-INFO-all \
--stdout | bgzip -c > tmp/$(basename ${truth%.*}).indel.gz
tabix -p vcf tmp/$(basename ${truth%.*}).indel.gz

# Get the number of invariant bases of the pipeline that are also invariant in the truth
# For indel and snp separately
for i in snp indel; do
    envs/var_call_v1/bin/bcftools isec -p ${out}_${i} -Oz tmp/$(basename ${raw%.*}).${i}.neg.gz \
    tmp/$(basename ${truth%.*}).${i}.gz

    tn=$(gunzip -c ${out}_${i}/0000.vcf.gz | grep -v "^#" | wc -l)

    if [ ${i} == "snp" ]
    then
        # Grab SNP FP from input table
        fp=$(cut -f 3 $tbl | tail -n +2 | head -1)
    else
        fp=$(cut -f 3 $tbl | tail -1)
    fi

    den_spec=$(( tn + fp ))
    specificity=$(echo "scale=4; $tn/$den_spec" | bc)

    printf "tTN\tSpecificity\n${tn}\t${specificity}\n" > ${tbl%.*}.${i}.tsv
done
