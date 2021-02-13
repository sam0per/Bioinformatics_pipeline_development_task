rule bcfstats:
    input:
        ref="data/{sample}.fa",
        samt="results/2_variant_calling/{sample}_filt.samtools.vcf.gz",
        gatk="results/2_variant_calling/{sample}_filt.gatk.vcf.gz"
    output: "results/qc/{sample}_bcfstats_callers.txt"
    params:
        bcf=config["modules"]["bcf"]
    shell:
        """
        {params.bcf} stats --depth 0,300,50 --fasta-ref {input.ref} {input.gatk} {input.samt} > {output}
        """

rule concordance:
    input:
        ref="data/{sample}.fa",
        tru=config["truth"],
        bed="results/2_variant_calling/{sample}_target.interval_list",
        eva="results/2_variant_calling/{sample}_filt.{caller}.vcf.gz"
    output:
        con="results/qc/{sample}_{caller}.concord.truth.tsv"
    shell:
        """
        gatk Concordance -R {input.ref} -summary {output.con} -eval {input.eva} --truth {input.tru} \
        -L {input.bed}
        """

rule specificty:
    input:
        tru=config["truth"],
        raw="results/2_variant_calling/{sample}_11_raw.{caller}.vcf.gz",
        con="results/qc/{sample}_{caller}.concord.truth.tsv"
    output:
        stat="results/3_performance/{sample}_{caller}.scores.tsv"
    params:
        edi="results/3_performance/{sample}_{caller}_neg/",
        tbl=lambda wildcards, input: os.path.splitext(input.con)[0]
        # bams=lambda wildcards, input: " ".join(input.con)
    shell:
        """
        bash ./scripts/edit_concordance.sh {input.raw} {input.tru} {params.edi} {input.con}
        cat <(printf "TN\tSPECIFICITY\n") <(awk 'FNR>1' {params.tbl}.snp.tsv {params.tbl}.indel.tsv) | \
        paste {input.con} - > {output.stat}
        """
