rule multiqc:
    input: "results/qc/{sample}_dedup.metrics.txt", "results/qc/{sample}_samtools.metrics.txt", "results/qc/{sample}_hs_metrics.txt"
    output: "qc/{sample}_multiqc.html"
    shell:
        """
        multiqc -n {output} results/qc
        """