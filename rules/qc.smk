rule multiqc:
    input:
        "results/qc/{sample}_dedup.metrics.txt",
        "results/qc/{sample}_samtools.metrics.txt",
        "results/qc/{sample}_hs_metrics.txt",
        "results/1_read_coverage/deep_coverage_{sample}.txt",
        "results/1_read_coverage/deep_{sample}_11_coverage.txt",
        "results/1_read_coverage/deep_{sample}_00_coverage.txt",
        "results/qc/{sample}_gc_freq.txt",
        "results/qc/GC_{sample}_11_freq.txt",
        "results/qc/GC_{sample}_00_freq.txt"
    output: "qc/{sample}_multiqc.html"
    shell:
        """
        multiqc -n {output} results/qc
        """