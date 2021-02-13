rule annotate_truth:
    input:
        vcf=config["truth"]
    output:
        ann="results/annotation/ground_truth.ann.vcf",
        csv="results/qc/ground_truth.ann.csv"
    params:
        ref=config["ref"]["release"],
        eff=config["modules"]["snpeff"]
    shell:
        """
        java -Xmx8g -jar {params.eff} -stats {output}.html -csv {output.csv} {params.ref} {input.vcf} > {output.ann}
        """