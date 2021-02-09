rule chr19_depth:
    input: "../Files_needed_for_task/chr19_dedup_sort.bam"
    output: "results/1_read_coverage/chr19_depth.txt"
    shell:
        """
        samtools index {input}
        samtools depth -a -q 20 -Q 30 -r 19:60004-14992498 {input} > {output}
        """

rule chr19_coverage:
    input: "results/1_read_coverage/chr19_depth.txt"
    output: "results/1_read_coverage/chr19_cov_redundancy.txt"
    shell:
        """
        cut -f 3 {input} | sort | uniq -c > {output}
        awk '{{ sum += $3 }} END {{ if (NR > 0) print "Mean coverage: " sum / NR "X" }}' {input}
        """

rule DeepT_coverage:
    input: "../Files_needed_for_task/chr19_dedup_sort.bam"
    output:
        bam="../Files_needed_for_task/chr19_rehead.bam",
        fig="figures/1_read_coverage/chr19_coverage",
        dat="results/1_read_coverage/chr19_coverage.txt"
    shell:
        """
        bash scripts/samtools_rehead.sh {input} {output.bam}
        samtools index {output.bam}
        plotCoverage -b {output.bam} --plotFile {output.fig} --ignoreDuplicates \
        --minMappingQuality 30 -r chr19:60004:14992498 --outRawCounts {output.dat}
        """

rule intersect:
    input:
        bam="../Files_needed_for_task/chr19_rehead.bam",
        bed="../Files_needed_for_task/target_regions.bed"
    output:
        intar="../Files_needed_for_task/chr19_intarget.bam",
        outar="../Files_needed_for_task/chr19_outarget.bam"
    shell:
        """
        awk -v OFS='\\t''$1="chr"$1' {input.bed} > ../Files_needed_for_task/target_regions_chr.bed
        bedtools intersect -abam {input.bam} -b ../Files_needed_for_task/target_regions_chr.bed > {output.intar}
        bedtools intersect -abam {input.bam} -b ../Files_needed_for_task/target_regions_chr.bed -v > {output.outar}
        """

rule DeepT_coverage_target:
    input:
        inbam="../Files_needed_for_task/chr19_intarget.bam",
        oubam="../Files_needed_for_task/chr19_outarget.bam"
    output:
        fig="figures/1_read_coverage/chr19_coverage_target",
        dat="results/1_read_coverage/chr19_coverage_target.txt"
    shell:
        """
        samtools index {input.inbam}
        samtools index {input.oubam}
        plotCoverage -b {input.inbam} {input.oubam} --plotFile {output.fig} --ignoreDuplicates \
        --minMappingQuality 30 --outRawCounts {output.dat}
        """

rule plot_coverage:
    input:
        xxx="results/1_read_coverage/chr19_coverage.txt",
        tar="results/1_read_coverage/chr19_coverage_target.txt"
    output: "figures/1_read_coverage/chr19_r_coverage"
    shell:
        """
        conda activate r_env
        Rscript scripts/chr19_1_plot_coverage.R -x {input.xxx} -t {input.tar} -o {output}
        conda deactivate
        """

rule hs_metrics:
    input:
        bed="../Files_needed_for_task/target_regions_chr.bed",
        bam="../Files_needed_for_task/chr19_rehead.bam"
    output: "results/1_read_coverage/target_regions_chr.interval_list"
    shell:
        """
        picard BedToIntervalList I={input.bed} O={output} SD={input.bam}
        """