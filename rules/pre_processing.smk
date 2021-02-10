rule mark_duplicates:
    input: "../Files_needed_for_task/chr19.bam"
    output:
        bam="data/chr19_dedup.bam",
        metrics="results/qc/chr19_dedup.metrics.txt",
        sortbam="data/chr19_dedup_sort.bam"
    log: "logs/picard/chr19_dedup.log"
    shell:
        """
        picard MarkDuplicates I={input} O={output.bam} M={output.metrics} \
        REMOVE_DUPLICATES=true ASSUME_SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT \
        USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
        picard SortSam INPUT={output.bam} OUTPUT={output.sortbam} SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT \
        USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
        samtools index {output.sortbam}
        """

rule base_recal:
    input: "data/chr19_dedup_sort.bam"
    output:
        tbl="results/preprocess/recal_data.table",
        bam="data/chr19_recal.bam"
    shell:
        """
        # wget -t 100 -O data/common_all_20180423.vcf.gz https://ftp.ncbi.nih.gov/snp/organisms/human_9606/VCF/common_all_20180423.vcf.gz
        # gunzip data/common_all_20180423.vcf.gz
        # gatk IndexFeatureFile -I data/common_all_20180423.vcf
        # rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz ./data/
        # gunzip data/hg19.fa.gz
        # sed -i -e 's/^>chr/>/g' data/hg19.fa.gz
        # gatk BaseRecalibrator -I {input} -R data/hg19.fa --known-sites data/common_all_20180423.vcf \
        # -O {output.tbl}
        # gatk ApplyBQSR -R data/hg19.fa -I {input} --bqsr-recal-file {output.tbl} -O {output.bam}
        # gatk AnalyzeCovariates -bqsr {output.tbl} -plots figures/AnalyzeCovariates.pdf
        """

rule chr19_interval:
    input: "data/chr19_recal.bam"
    output: "results/1_read_coverage/chr19_interval.txt"
    shell: 
        """
        samtools view {input} | awk '{{print $3 "\\t" $4 "\\t" $4+length($10)-1}}' > {output}
        sort -nk2 {output} | awk 'NR==1{{print "min:" $2}} END{{print "max:" $3}}'
        """
