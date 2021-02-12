rule annotate_truth:
    input:
        vcf=config["truth"]
    output: "results/annotation/ground_truth.ann.vcf"
    params:
        ref=config["ref"]["release"],
        eff=config["modules"]["snpeff"]
    shell:
        """
        java -Xmx8g -jar {params.eff} -stats {output}.html {params.ref} {input.vcf} > {output}
        """