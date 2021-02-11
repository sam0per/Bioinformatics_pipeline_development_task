rule mark_duplicates:
    input: "../Files_needed_for_task/{sample}.bam".format(sample=config["sample"])
    output:
        bam="data/{sample}_dedup.bam",
        metrics="results/qc/{sample}_dedup.metrics.txt",
        sortbam="data/{sample}_dedup_sort.bam"
    log: "logs/picard/{sample}_dedup.log"
    shell:
        """
        picard MarkDuplicates I={input} O={output.bam} M={output.metrics} \
        REMOVE_DUPLICATES=true ASSUME_SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT \
        USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
        
        # QC-passed properly paired
        samtools view -b -f 2 -F 1804 {output.bam} > data/{wildcards.sample}_ppqc_nodup.bam
        picard SortSam INPUT=data/{wildcards.sample}_ppqc_nodup.bam OUTPUT={output.sortbam} SORT_ORDER=coordinate \
        VALIDATION_STRINGENCY=LENIENT USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
        samtools index {output.sortbam}
        """

rule fasta:
    input: "data/{sample}_dedup_sort.bam"
    output:
        fas="data/{sample}.fa",
        bam="data/{sample}_RG.bam"
    params:
        rel=config["ref"]["release"],
        ref=config["ref"]["site"]
    shell:
        """
        # Download sequence fasta
        # rsync -avzP {params.ref} ./data/
        gunzip -c data/{wildcards.sample}.fa.gz | sed -e 's/^>chr/>/g' > {output.fas}
        samtools faidx {output.fas}
        gatk CreateSequenceDictionary -R {output.fas}
        
        # Add read groups and build index
        picard AddOrReplaceReadGroups I={input} O={output.bam} \
        RGLB=lib1 RGPL=ILLUMINA RGPU=unit1 RGSM=XXX VALIDATION_STRINGENCY=LENIENT USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
        picard BuildBamIndex I={output.bam} VALIDATION_STRINGENCY=LENIENT USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
        """

rule interval:
    input: "data/{sample}_RG.bam"
    output: "results/1_read_coverage/{sample}_interval.txt"
    shell: 
        """
        samtools view {input} | awk '{{print $3 "\\t" $4 "\\t" $4+length($10)-1}}' > {output}
        sort -nk2 {output} | awk 'NR==1{{print "min:" $2}} END{{print "max:" $3}}'
        """

rule stats:
    input:
        bam="data/{sample}_RG.bam",
        fas="data/{sample}.fa"
    output: "results/qc/{sample}_samtools.metrics.txt"
    params:
        reg=config["ref"]["reg"]
    shell:
        """
        samtools stats -r {input.fas} {input.bam} {params.reg} > {output}
        """