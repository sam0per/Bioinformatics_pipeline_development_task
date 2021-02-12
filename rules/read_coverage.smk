rule all_depth:
    input: "data/{sample}_RG.bam"
    output: "results/1_read_coverage/{sample}_depth.txt"
    params:
        reg=config["ref"]["reg"]
    shell:
        """
        samtools depth -a -q 20 -Q 30 -r {params.reg} {input} > {output}
        awk '{{ sum += $3 }} END {{ if (NR > 0) print "Overall mean coverage: " sum / NR "X" }}' {output}
        """

rule DeepT_coverage:
    input: "data/{sample}_RG.bam"
    output:
        # fig="figures/1_read_coverage/{sample}_coverage",
        dat="results/1_read_coverage/deep_coverage_{sample}.txt"
    params:
        reg=os.path.basename(config["ref"]["reg"]).replace("-", ":")
    shell:
        """
        plotCoverage -b {input} --ignoreDuplicates \
        --minMappingQuality 30 -r {params.reg} --outRawCounts {output.dat} -n 100000 \
        -p max/2
        """

rule intersect:
    input:
        bam="data/{sample}_RG.bam",
        bed=config["target"]
    output:
        inbam="results/1_read_coverage/{sample}_11.bam",
        ofbam="results/1_read_coverage/{sample}_00.bam",
        ofbed="results/1_read_coverage/{sample}_00.bed"
    shell:
        # awk -v OFS='\\t''$1="chr"$1' {input.bed} > ../Files_needed_for_task/target_regions_chr.bed
        """
        bedtools intersect -abam {input.bam} -b {input.bed} > {output.inbam}
        bedtools intersect -abam {input.bam} -b {input.bed} -v > {output.ofbam}
        bedtools intersect -abam {input.bam} -b {input.bed} -v -bed | \
        cut -f 1,2,3 > {output.ofbed}
        samtools index {output.inbam}
        samtools index {output.ofbam}
        """

rule DeepT_target:
    input:
        inbam="results/1_read_coverage/{sample}_11.bam",
        inbed=config["target"],
        ofbam="results/1_read_coverage/{sample}_00.bam"
    output:
        indat="results/1_read_coverage/deep_{sample}_11_coverage.txt",
        ofdat="results/1_read_coverage/deep_{sample}_00_coverage.txt"
    params:
        reg=os.path.basename(config["ref"]["reg"]).replace("-", ":")
    shell:
        """
        # In-target
        plotCoverage -b {input.inbam} --ignoreDuplicates \
        --minMappingQuality 30 -r {params.reg} --outRawCounts {output.indat} -n 100000 --BED {input.inbed} \
        -p max/2
        
        # Of-target
        plotCoverage -b {input.ofbam} --ignoreDuplicates \
        --minMappingQuality 30 -r {params.reg} --outRawCounts {output.ofdat} -n 100000 --BED {input.ofbam}.bed \
        -p max/2
        """

rule plot_coverage:
    input:
        who="results/1_read_coverage/deep_coverage_{sample}.txt",
        inr="results/1_read_coverage/deep_{sample}_11_coverage.txt",
        our="results/1_read_coverage/deep_{sample}_00_coverage.txt"
    output: "figures/1_read_coverage/{sample}_r_coverage.png"
    params:
        r=config["modules"]["r"]
    shell:
        """
        {params.r} scripts/chr19_1_plot_coverage.R -x {input.who} -t {input.inr} -f {input.our} -o {output}
        """

rule hs_metrics:
    input:
        bed=config["target"],
        bam="data/{sample}_RG.bam"
    output:
        ls="results/1_read_coverage/{sample}_target_regions.interval_list",
        hs="results/qc/{sample}_hs_metrics.txt"
    shell:
        """
        picard BedToIntervalList I={input.bed} O={output.ls} SD={input.bam} VALIDATION_STRINGENCY=LENIENT USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
        picard CollectHsMetrics I={input.bam} O={output.hs} BAIT_INTERVALS={output.ls} \
        TARGET_INTERVALS={output.ls} VALIDATION_STRINGENCY=LENIENT USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
        """