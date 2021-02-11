rule gc:
    input:
        bam="data/{sample}_RG.bam",
        ref="chrom/{sample}_hg19.2bit".format(sample=config["sample"])
    output:
        fig="figures/1_read_coverage/{sample}_gc.png",
        dat="results/qc/{sample}_gc_freq.txt",
        obam="results/1_read_coverage/gc/{sample}_GC_corr.bam"
    params:
        # bedtools genomecov -ibam data/chr19_RG.bam -bga | awk '{s+=$4}END{print s}' 
        # gen=len(open("results/1_read_coverage/{sample}_depth.txt").readlines())
        reg=os.path.basename(config["ref"]["reg"]).replace("-", ":")
    shell:
        """
        computeGCBias -b {input.bam} --effectiveGenomeSize 125954048 \
        -g {input.ref} --GCbiasFrequenciesFile {output.dat} \
        -r {params.reg} --biasPlot {output.fig} -p max/2

        correctGCBias -b {input.bam} --effectiveGenomeSize 125954048 -g {input.ref} \
        --GCbiasFrequenciesFile {output.dat} -o {output.obam} -p max/2
        """

rule target_gc:
    input:
        inbam="results/1_read_coverage/{sample}_11.bam",
        ref="chrom/{sample}_hg19.2bit",
        bed="../Files_needed_for_task/target_regions.bed"
    output:
        fig="figures/1_read_coverage/GC_{sample}_11.png",
        dat="results/qc/GC_{sample}_11_freq.txt",
        obam="results/1_read_coverage/gc/GC_{sample}_11_corr.bam"
    params:
        reg=os.path.basename(config["ref"]["reg"]).replace("-", ":")
    shell:
        """
        computeGCBias -b {input.inbam} --effectiveGenomeSize 111525929 -g {input.ref} \
        --extraSampling {input.bed} --GCbiasFrequenciesFile {output.dat} \
        -r {params.reg} --biasPlot {output.fig} -p max/2

        correctGCBias -b {input.inbam} --effectiveGenomeSize 111525929 -g {input.ref} \
        --GCbiasFrequenciesFile {output.dat} -o {output.obam} -p max/2
        """

rule oftarget_gc:
    input:
        ofbam="results/1_read_coverage/{sample}_00.bam",
        ref="chrom/{sample}_hg19.2bit",
        bed="results/1_read_coverage/{sample}_00.bam.bed"
    output:
        fig="figures/1_read_coverage/GC_{sample}_00.png",
        dat="results/qc/GC_{sample}_00_freq.txt",
        obam="results/1_read_coverage/gc/GC_{sample}_00_corr.bam"
    params:
        reg=os.path.basename(config["ref"]["reg"]).replace("-", ":")
    shell:
        """
        computeGCBias -b {input.ofbam} --effectiveGenomeSize 7817187 -g {input.ref} \
        --extraSampling {input.bed} --GCbiasFrequenciesFile {output.dat} \
        -r {params.reg} --biasPlot {output.fig} -p max/2

        correctGCBias -b {input.ofbam} --effectiveGenomeSize 7817187 -g {input.ref} \
        --GCbiasFrequenciesFile {output.dat} -o {output.obam} -p max/2
        """