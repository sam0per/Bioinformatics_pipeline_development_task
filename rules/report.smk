rule create_html:
    input: "Report_pipeline_development_Samuel.Rmd"
    output: "html"
    shell: "bash scripts/render_rmd.sh {input} {output}"