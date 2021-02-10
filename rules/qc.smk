rule multiqc:
    input: "results/qc/chr19_dedup.metrics.txt"
    output: "qc/multiqc.html"
    shell:
        """
        multiqc -n {output} results/qc
        """