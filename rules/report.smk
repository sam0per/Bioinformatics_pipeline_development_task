rule create_html:
    input:
        rmd="Report_pipeline_development_Samuel.Rmd",
        tbl=expand("results/3_performance/{sample}_{caller}.scores.tsv", sample=config["sample"], caller=callers)
    output: "html"
    shell: "bash scripts/render_rmd.sh {input.rmd} {output}"