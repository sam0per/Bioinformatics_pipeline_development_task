# Analysis Assignment – Pipeline Development

Author: Samuel Perini

With this assignment, I aim to demonstrate my exploratory analysis skills and methodology in three main subjects:

__1. Read coverage.__

__2. Variant calling analysis.__

__3. Analytical performance.__

A detailed description of each of the three subjects is presented below and it is based on the questions that were given for this assignment. Additional points were included to clarify analytical steps and results that might be less accessible to those that are less familiar with DNA sequencing and variant calling.  
The pipeline was built using [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html), a Python-based workflow management system, and after [cloning](https://help.github.com/en/articles/cloning-a-repository) this repository, the same bioinformatics steps can be performed by typing the following command line codes:

1. [Install Miniconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

2. Create an environment called "pipeline_development" with the required software listed inside the YAML file and activate it
```
    cd Bioinformatics_pipeline_development_task
    conda env create --name pipeline_development --file environment.yaml
    conda activate pipeline_development
```

3. Download human assembly hg19 and chromosome 19 fasta file
```
    rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit .
    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr19.fa.gz
    gunzip chr19.fa.gz
```

4. Execute the workflow locally using 1 core (change the number of cores as you wish)
```
    snakemake --cores 1
```

5. Output a self-contained interactive HTML report with all results
```
    snakemake --report report.html
```
