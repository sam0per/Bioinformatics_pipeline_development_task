##### Target rules #####

rule all:
    input:
        "results/annotation/ground_truth.ann.vcf",
        "qc/multiqc.html",
        "results/1_read_coverage/chr19_interval.txt"
        # "figures/1_read_coverage/chr19_coverage",
        # "figures/1_read_coverage/chr19_coverage_target",
        # "figures/1_read_coverage/chr19_gc.png",
        # "figures/1_read_coverage/chr19_r_coverage.png",
        # "html"

##### Modules #####

include: "rules/annotation.smk"
include: "rules/pre_processing.smk"
# include: "rules/read_coverage.smk"
# include: "rules/gc_content.smk"
# include: "rules/calling.smk"
include: "rules/qc.smk"
# include: "rules/report.smk"
