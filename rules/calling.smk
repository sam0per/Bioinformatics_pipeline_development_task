rule mpileup:
    input: 
    output: 
    shell:
        """
        conda activate variant_calling1
        bcftools mpileup -f data/chr19.fa ../Files_needed_for_task/chr19.bam | bcftools call -mv -Ob -o test/calls.bcf
        conda deactivate
        """

rule name:
    input: 
    output: 
    run: 