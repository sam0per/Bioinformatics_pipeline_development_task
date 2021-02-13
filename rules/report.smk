rule create_html:
    input:
        rmd="Report_pipeline_development_Samuel.Rmd",
        tbl=expand("results/3_performance/{sample}_{caller}.scores.tsv", sample=config["sample"], caller=callers),
        fig=expand("figures/1_read_coverage/GC_{sample}_00.png", sample=config["sample"])
    output: "html"
    shell: "bash scripts/render_rmd.sh {input.rmd} {output}"