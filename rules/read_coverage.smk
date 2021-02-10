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

import os
rule DeepT_coverage:
    input: "data/{sample}_RG.bam"
    output:
        fig="figures/1_read_coverage/{sample}_coverage",
        dat="results/1_read_coverage/deep_{sample}_coverage.txt"
    params:
        reg=os.path.basename(config["ref"]["reg"]).replace("-", ":")
    shell:
        """
        plotCoverage -b {input} --plotFile {output.fig} --ignoreDuplicates \
        --minMappingQuality 30 -r {params.reg} --outRawCounts {output.dat} -n 100000
        """

rule intersect:
    input:
        bam="data/{sample}_RG.bam",
        bed="../Files_needed_for_task/target_regions.bed"
    output:
        inbam="results/1_read_coverage/{sample}_intarget.bam",
        inbed="results/1_read_coverage/{sample}_intarget.bed",
        oubed="results/1_read_coverage/{sample}_outarget.bed"
    shell:
        """
        # awk -v OFS='\\t''$1="chr"$1' {input.bed} > ../Files_needed_for_task/target_regions_chr.bed
        bedtools intersect -abam {input.bam} -b {input.bed} > {output.inbam}
        bedtools intersect -abam {input.bam} -b {input.bed} -bed > {output.inbed}
        bedtools intersect -abam {input.bam} -b {input.bed} -v > {output.oubed}
        
        # Mean depth for on- and off-target
        bedtools coverage -a {input.bed} -b {input.bam} -mean \
        > results/1_read_coverage/{wildcards.sample}_intarget_mean.bedgraph
        bedtools coverage -a {output.oubed} -b {input.bam} -mean \
        > results/1_read_coverage/{wildcards.sample}_outarget_mean.bedgraph
        """

rule plot_coverage:
    input:
        who="results/1_read_coverage/deep_{sample}_coverage.txt",
        tar="results/1_read_coverage/{sample}_intarget.bed"
    output: "figures/1_read_coverage/{sample}_r_coverage.png"
    params:
        r=config["modules"]["r"]
    shell:
        """
        {params.r} scripts/chr19_1_plot_coverage.R -x {input.who} -t {input.tar} -o {output}
        """

rule hs_metrics:
    input:
        bed="../Files_needed_for_task/target_regions.bed",
        bam="data/{sample}_RG.bam"
    output:
        ls="results/1_read_coverage/{sample}_target_regions.interval_list",
        hs="results/qc/{sample}_hs_metrics.txt"
    shell:
        """
        picard BedToIntervalList I={input.bed} O={output.ls} SD={input.bam}
        picard CollectHsMetrics I={input.bam} O={output.hs}
        """