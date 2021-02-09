rule chr19_gc:
    input: "../Files_needed_for_task/chr19_rehead.bam"
    output:
        fig="figures/1_read_coverage/chr19_gc.png",
        dat="results/1_read_coverage/chr19_gc_freq.txt"
    shell:
        """
        computeGCBias -b {input} --effectiveGenomeSize 2864785220 \
        -g data/hg19.2bit --GCbiasFrequenciesFile {output.dat} \
        -r chr19:60004:14992498 --biasPlot {output.fig}
        """

rule target_gc:
    input: "../Files_needed_for_task/chr19_intarget.bam"
    output:
        fig="figures/1_read_coverage/intarget_gc.png",
        dat="results/1_read_coverage/intarget_gc_freq.txt"
    shell:
        """
        computeGCBias -b {input} --effectiveGenomeSize 2864785220 -g data/hg19.2bit \
        --extraSampling ../Files_needed_for_task/target_regions.bed --GCbiasFrequenciesFile {output.dat} \
        -r chr19:60004:14992498 --biasPlot {output.fig}
        """