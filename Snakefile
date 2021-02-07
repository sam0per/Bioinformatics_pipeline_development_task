rule all:
    input:
        "figures/1_read_coverage/chr19_coverage",
        "figures/1_read_coverage/chr19_coverage_target",
        "figures/1_read_coverage/chr19_gc.png"


rule mark_duplicates:
    input: "../Files_needed_for_task/chr19.bam"
    output:
        bam="../Files_needed_for_task/chr19_dedup.bam",
        metrics="results/qc/chr19_dedup.metrics.txt",
        sortbam="../Files_needed_for_task/chr19_dedup_sort.bam"
    log: "logs/picard/chr19_dedup.log"
    shell:
        """
        picard MarkDuplicates I={input} O={output.bam} M={output.metrics} \
        REMOVE_DUPLICATES=true ASSUME_SORT_ORDER=coordinate VALIDATION_STRINGENCY=LENIENT \
        USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
        picard SortSam INPUT={output.bam} OUTPUT={output.sortbam} SORT_ORDER=coordinate
        """

rule chr19_interval:
    input: "../Files_needed_for_task/chr19_dedup_sort.bam"
    output: "results/1_read_coverage/chr19_interval.txt"
    shell: 
        """
        samtools view {input} | awk '{{print $3 "\\t" $4 "\\t" $4+length($10)-1}}' > {output}
        sort -nk2 {output} | awk 'NR==1{{print "min:" $2}} END{{print "max:" $3}}'
        """

rule chr19_depth:
    input: "../Files_needed_for_task/chr19_dedup_sort.bam"
    output: "results/1_read_coverage/chr19_depth.txt"
    shell:
        """
        samtools index {input}
        samtools depth -a -q 10 -Q 10 -r 19:60004-14992498 {input} > {output}
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
        --minMappingQuality 10 -r chr19:60004:14992498 --outRawCounts {output.dat}
        """

rule DeepT_coverage_target:
    input:
        bam="../Files_needed_for_task/chr19_rehead.bam",
        bed="../Files_needed_for_task/target_regions.bed"
    output:
        fig="figures/1_read_coverage/chr19_coverage_target",
        dat="results/1_read_coverage/chr19_coverage_target.txt"
    shell:
        """
        plotCoverage -b {input.bam} --plotFile {output.fig} --BED {input.bed} --ignoreDuplicates \
        --minMappingQuality 10 -r chr19:60004:14992498 --outRawCounts {output.dat}
        """

rule plot_coverage:
    input:
        red="results/1_read_coverage/chr19_cov_redundancy.txt",
        cov="results/1_read_coverage/chr19_coverage.txt"
    output: "figures/1_read_coverage/chr19_r_coverage.png"
    shell:
        """
        conda activate r_env
        Rscript scripts/chr19_1_plot_coverage.R -r {input.red} -b {input.cov} -o {output}
        conda deactivate
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

rule base_recal:
    input: "../Files_needed_for_task/chr19_dedup_sort.bam"
    output:
        tbl="../Files_needed_for_task/recal_data.table",
        bam="../Files_needed_for_task/chr19_recal.bam"
    shell:
        """
        gatk BaseRecalibrator -I {input} -R data/hg19.fa --known-sites ../Files_needed_for_task/ground_truth.vcf \
        -O {output.tbl}
        gatk ApplyBQSR -R data/hg19.fa -I {input} --bqsr-recal-file {output.tbl} -O {output.bam}
        gatk AnalyzeCovariates -bqsr {output.tbl} -plots figures/AnalyzeCovariates.pdf
        """