rule chr19_interval:
    input: "../Files_needed_for_task/chr19.bam"
    output: "../Read_coverage/chr19_interval.txt"
    shell: 
        """
        samtools view {input} | awk '{{print $3 "\\t" $4 "\\t" $4+length($10)-1}}' > {output}
        sort -nk2 {output} | awk 'NR==1{{print "min:" $2}} END{{print "max:" $3}}'
        """

rule chr19_depth:
    input: "../Files_needed_for_task/chr19.bam"
    output: "1_read_coverage/chr19_depth.txt"
    shell: "samtools depth {input} > {output}"

rule chr19_coverage:
    input: "1_read_coverage/chr19_depth.txt"
    output: "1_read_coverage/chr19_cov_redundancy.txt"
    shell: "cut -f 3 {input} | sort | uniq -c > {output}"

rule chr19_gc:
    input: 
    output: 
    run: 