import os
configfile: "config.yaml"

callers=["gatk", "samtools"]

##### Target rules #####

rule all:
    input:
        "results/annotation/ground_truth.ann.vcf",
        expand("results/1_read_coverage/{sample}_interval.txt", sample=config["sample"]),
        expand("results/1_read_coverage/{sample}_depth.txt", sample=config["sample"]),
        expand("figures/1_read_coverage/{sample}_r_coverage.png", sample=config["sample"]),
        expand("figures/1_read_coverage/GC_{sample}_00.png", sample=config["sample"]),
        expand("results/qc/{sample}_filt.{caller}.variant_calling_detail_metrics", sample=config["sample"], caller=callers),
        expand("results/qc/{sample}_{caller}.concord.truth.tsv", sample=config["sample"], caller=callers),
        expand("results/3_performance/{sample}_{caller}.scores.tsv", sample=config["sample"], caller=callers),
        expand("qc/{sample}_multiqc.html", sample=config["sample"]),
        "html"

##### Modules #####

include: "rules/annotation.smk"
include: "rules/pre_processing.smk"
include: "rules/read_coverage.smk"
include: "rules/gc_content.smk"
include: "rules/calling.smk"
include: "rules/performance.smk"
include: "rules/qc.smk"
include: "rules/report.smk"
