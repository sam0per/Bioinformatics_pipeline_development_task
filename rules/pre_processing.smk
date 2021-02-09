rule mark_duplicates:
    input: "../Files_needed_for_task/chr19.bam"
    output:
        bam="../Files_needed_for_task/chr19_dedup.bam",
        metrics="results/qc/chr19_dedup.metrics.txt",
        sortbam="../Files_needed_for_task/chr19_dedup_sort.bam"
    log: "logs/picard/chr19_dedup.log"
    shell:
        """
        picard MarkDuplicates I={input} O={output.bam} M={output.metrics} \
        REMOVE_DUPLICATES=true ASSUME_SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT \
        USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
        picard SortSam INPUT={output.bam} OUTPUT={output.sortbam} SORT_ORDER=coordinate
        """

rule chr19_interval:
    input: "../Files_needed_for_task/chr19_dedup_sort.bam"
    output: "results/1_read_coverage/chr19_interval.txt"
    shell: 
        """
        samtools view {input} | awk '{{print $3 "\\t" $4 "\\t" $4+length($10)-1}}' > {output}
        sort -nk2 {output} | awk 'NR==1{{print "min:" $2}} END{{print "max:" $3}}'
        """
