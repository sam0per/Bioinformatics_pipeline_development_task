rule gc:
    input: "data/{sample}_RG.bam"
    output: "results/qc/{sample}_dedup.GC"
    params:
        ref=config["ref"]["complete"],
        fas="data/" + config["ref"]["release"] + ".fa.gz"
    shell:
        """
        rsync -avzP {params.ref} ./data/
        gatk CreateSequenceDictionary -R {params.fas}

        picard CollectGcBiasMetrics I={input} CHART={output}.pdf O={output} \
        S={output}.summary.txt R={params.fas} VALIDATION_STRINGENCY=LENIENT USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
        """

rule target_gc:
    input:
        inbam="results/1_read_coverage/{sample}_11.bam",
        ref="chrom/{sample}_hg19.2bit",
        bed=config["target"]
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

rule GC_metrics:
    input:
        gcbam="results/1_read_coverage/{sample}_11.bam"
    output:
        mtx="results/qc/{sample}_11_GC.bias.txt",
        summ="results/qc/{sample}_11_GC.summary.txt"
    params:
        fas="data/" + config["ref"]["release"] + ".fa.gz"
    shell:
        """
        picard CollectGcBiasMetrics I={input.gcbam} CHART=gc_bias_metrics.pdf O={output.mtx} \
        S={output.summ} R={params.fas} VALIDATION_STRINGENCY=LENIENT USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
        """

rule GCcorr_metrics:
    input:
        nogc="results/1_read_coverage/gc/GC_{sample}_11_corr.bam"
    output:
        mtx="results/qc/GC_{sample}_corr.metrics.txt",
        summ="results/qc/GC_{sample}_corr.gcbias.txt"
    params:
        fas="data/" + config["ref"]["release"] + ".fa.gz"
    shell:
        """
        picard CollectGcBiasMetrics I={input.nogc} CHART=gc_bias_metrics.pdf O={output.mtx} \
        S={output.summ} R={params.fas} VALIDATION_STRINGENCY=LENIENT USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
        """

rule oftarget_gc:
    input:
        ofbam="results/1_read_coverage/{sample}_00.bam",
        ref="chrom/{sample}_hg19.2bit",
        bed="results/1_read_coverage/{sample}_00.bed"
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