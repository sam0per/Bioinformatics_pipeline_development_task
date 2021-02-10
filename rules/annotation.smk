rule annotate_truth:
    input:
        vcf="../Files_needed_for_task/ground_truth.vcf"
    output: "results/annotation/ground_truth.ann.vcf"
    shell:
        """
        java -Xmx8g -jar $HOME/miniconda/envs/variant_calling1/share/snpeff-5.0-0/snpEff.jar \
        -stats results/annotation/ground_truth.html hg19 {input.vcf} > {output}
        """