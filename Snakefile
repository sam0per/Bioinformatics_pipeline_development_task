import os
import pandas as pd
configfile: "config.yaml"

# callers_tb = pd.read_table(config["caller_file"], dtype=str).set_index("caller", drop=False)
# caller=callers_tb.index.values

# callers_names = {"gatk": "GATK", "samt": "SAMtools"}
# callers=callers_names.keys()

callers=["gatk", "samtools"]

# callers="|".join(vcallers.index)

##### Target rules #####

rule all:
    input:
        expand("results/1_read_coverage/{sample}_interval.txt", sample=config["sample"]),
        expand("results/1_read_coverage/{sample}_depth.txt", sample=config["sample"]),
        # "results/annotation/ground_truth.ann.vcf",
        expand("figures/1_read_coverage/{sample}_r_coverage.png", sample=config["sample"]),
        expand("results/3_performance/{sample}_{caller}.scores.tsv", sample=config["sample"], caller=callers),
        expand("qc/{sample}_{caller}_multiqc.html", sample=config["sample"], caller=callers)
        # "figures/1_read_coverage/chr19_gc.png",
        # "html"

##### Modules #####

# include: "rules/annotation.smk"
include: "rules/pre_processing.smk"
include: "rules/read_coverage.smk"
include: "rules/gc_content.smk"
include: "rules/calling.smk"
include: "rules/performance.smk"
include: "rules/qc.smk"
# include: "rules/report.smk"
