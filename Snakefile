##### Target rules #####

rule all:
    input:
        "figures/1_read_coverage/chr19_coverage",
        "figures/1_read_coverage/chr19_coverage_target",
        "figures/1_read_coverage/chr19_gc.png",
        "figures/1_read_coverage/chr19_r_coverage.png"

##### Modules #####

include: "rules/pre_processing.smk"
include: "rules/read_coverage.smk"
include: "rules/gc_content.smk"
include: "rules/annotation.smk"







rule base_recal:
    input: "../Files_needed_for_task/chr19_dedup_sort.bam"
    output:
        tbl="../Files_needed_for_task/recal_data.table",
        bam="../Files_needed_for_task/chr19_recal.bam"
    shell:
        """
        gatk BaseRecalibrator -I {input} -R data/hg19.fa --known-sites ../Files_needed_for_task/ground_truth.vcf \
        -O {output.tbl}
        gatk ApplyBQSR -R data/hg19.fa -I {input} --bqsr-recal-file {output.tbl} -O {output.bam}
        gatk AnalyzeCovariates -bqsr {output.tbl} -plots figures/AnalyzeCovariates.pdf
        """