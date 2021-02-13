rule mpileup:
    input:
        ref="data/{sample}.fa",
        bam="results/1_read_coverage/gc/GC_{sample}_11_corr.bam"
    output:
        raw="results/2_variant_calling/{sample}_11_raw.samtools.vcf.gz",
        fil="results/2_variant_calling/{sample}_filt.samtools.vcf.gz"
    params:
        bcf=config["modules"]["bcf"],
        exc=config["filters"]["bcf"]
    shell:
        """
        {params.bcf} mpileup -a 'FORMAT/DV,FORMAT/AD,FORMAT/DP' -f {input.ref} {input.bam} | \
        {params.bcf} call -m -Oz -g 10 -o {output.raw}
        {params.bcf} filter -e 'ALT="."' {output.raw} | {params.bcf} filter -g3 -G10 \
        -e '{params.exc}' -Oz -o {output.fil}
        tabix -p vcf {output.fil}
        """

rule gatk_H:
    input:
        ref="data/{sample}.fa",
        bam="results/1_read_coverage/gc/GC_{sample}_11_corr.bam",
        bed=config["target"]
    output: "results/2_variant_calling/{sample}_11.g.vcf.gz"
    params:
        rgsm=config["samplename"],
        reg=config["dbSNP"]["chr"]
    shell:
        """
        picard BedToIntervalList I={input.bed} O={input.bed}.interval_list SD=data/{wildcards.sample}.dict
        gatk --java-options "-Xmx8G" HaplotypeCaller -R {input.ref} -I {input.bam} \
        -O {output} -ERC GVCF --sample-name {params.rgsm} -L {input.bed}.interval_list \
        --output-mode EMIT_ALL_CONFIDENT_SITES
        """

rule gatk_G:
    input:
        ref="data/{sample}.fa",
        bam="results/2_variant_calling/{sample}_11.g.vcf.gz"
    output:
        raw="results/2_variant_calling/{sample}_11_raw.gatk.vcf.gz",
        fil="results/2_variant_calling/{sample}_filt.gatk.vcf.gz"
    params:
        bed=config["target"] + ".interval_list",
        bcf=config["modules"]["bcf"],
        exc=config["filters"]["bcf"]
    shell:
        """
        gatk --java-options "-Xmx8G" GenotypeGVCFs -R {input.ref} -V {input.bam} -O {output.raw} -L {params.bed} \
        --include-non-variant-sites true
        {params.bcf} filter -e 'ALT="."' {output.raw} | {params.bcf} filter -g3 -G10 \
        -e '{params.exc}' -Oz -o {output.fil}
        tabix -p vcf {output.fil}
        """