# Analysis Assignment – Pipeline Development

Author: Samuel Perini

With this assignment, I aim to demonstrate my exploratory analysis skills and methodology in three main subjects:

__1. Read coverage.__

__2. Variant calling analysis.__

__3. Analytical performance.__

This repository contains the scripting files that are necessary for running the pipeline but it does not contain the input data that were given for the assignment. The input data should be downloaded and placed in a directory called `Files_needed_for_task` and this directory should be in the same folder as the repository folder (not inside the repository). For example, both `Files_needed_for_task` and `Bioinformatics_pipeline_development_task` (the repository folder name) directories can be placed inside a folder called `Assignment_Pipeline_Development`.  
The pipeline was built using [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html), a Python-based workflow management system, and after [cloning](https://help.github.com/en/articles/cloning-a-repository) this repository, the same bioinformatics steps can be performed by typing the following command line codes:

1. [Install Miniconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

1. Create an environment called "pipeline_development" with the required software listed inside the YAML file and activate it
```
    cd Bioinformatics_pipeline_development_task
    conda env create --name pipeline_development --file environment.yaml
    conda activate pipeline_development
```

1. Create a second environment called "r_env" with the required software listed inside a different YAML file. This environment will be activated inside the pipeline to avoid conflicting packages.
```
    conda env create --name r_env --file r-environment.yaml
```

1. Download human assembly hg19, chromosome 19 FASTA file and chromosome 19 dbSNP BED file
```
    mkdir data
    cd data
    rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit .
    rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz .
    gunzip hg19.fa.gz
    wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr19.fa.gz
    gunzip chr19.fa.gz
    wget ftp://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/BED/bed_chr_19.bed.gz
    cd ../
```

1. Execute the workflow locally using 1 core (change the number of cores as you wish)
```
    snakemake --cores 1
```

1. Output a self-contained interactive HTML report with all results
```
    snakemake --report report.html
```
