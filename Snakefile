rule chr19_interval:
    input: "../Files_needed_for_task/chr19.bam"
    output: "results/1_read_coverage/chr19_interval.txt"
    shell: 
        """
        samtools view {input} | awk '{{print $3 "\\t" $4 "\\t" $4+length($10)-1}}' > {output}
        sort -nk2 {output} | awk 'NR==1{{print "min:" $2}} END{{print "max:" $3}}'
        """

rule chr19_depth:
    input: "../Files_needed_for_task/chr19.bam"
    output: "results/1_read_coverage/chr19_depth.txt"
    shell: "samtools depth {input} > {output}"

rule chr19_coverage:
    input: "results/1_read_coverage/chr19_depth.txt"
    output: "results/1_read_coverage/chr19_cov_redundancy.txt"
    shell: "cut -f 3 {input} | sort | uniq -c > {output}"

rule plot_coverage:
    input: "../Files_needed_for_task/chr19_rehead.bam"
    output:
        fig="figures/1_read_coverage/chr19_coverage",
        dat="results/1_read_coverage/chr19_coverage.txt"
    shell:
        """
        bash scripts/samtools_rehead.sh
        samtools index {input}
        plotCoverage -b {input} --plotFile {output.fig} --ignoreDuplicates --minMappingQuality 10 -r chr19:60004:14992498 \
        --outRawCounts {output.dat}
        """

rule plot_coverage_target:
    input:
        bam="../Files_needed_for_task/chr19_rehead.bam",
        bed="../Files_needed_for_task/target_regions.bed"
    output:
        fig="figures/1_read_coverage/chr19_coverage_target",
        dat="results/1_read_coverage/chr19_coverage_target.txt"
    shell:
        """
        plotCoverage -b {input.bam} --plotFile {output.fig} --BED {input.bed} --ignoreDuplicates \
        --minMappingQuality 10 -r chr19:60004:14992498 \
        --outRawCounts {output.dat}
        """

rule chr19_gc:
    input: "../Files_needed_for_task/chr19_rehead.bam"
    output:
        fig="figures/1_read_coverage/chr19_gc.png",
        dat="results/1_read_coverage/GC_freq.txt"
    shell:
        """
        computeGCBias -b {input} --effectiveGenomeSize 2864785220 \
        -g data/hg19.2bit --GCbiasFrequenciesFile {output.dat} \
        -r chr19:60004:14992498 --biasPlot {output.fig}
        """