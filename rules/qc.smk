import glob

def qcstats(wildcards):
    return glob.glob("results/qc/wildcards.sample*")

rule multiqc:
    input:
        "results/qc/ground_truth.ann.csv",
        "results/qc/{sample}_dedup.metrics.txt",
        "results/qc/{sample}_dedup.metrics.txt",
        "results/qc/deep_coverage_{sample}.txt",
        "results/qc/deep_{sample}_11_coverage.txt",
        "results/qc/deep_{sample}_00_coverage.txt",
        "results/qc/{sample}_hs_metrics.txt",
        "results/qc/{sample}_dedup.GC",
        "results/qc/{sample}_11_GC.summary.txt",
        "results/qc/GC_{sample}_corr.gcbias.txt",
        "results/qc/{sample}_bcfstats_callers.txt",
        "results/qc/{sample}_bcfstats_callers.txt"
    output: "qc/{sample}_multiqc.html"
    shell:
        """
        multiqc -n {output} results/qc
        """