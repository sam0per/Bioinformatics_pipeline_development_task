# Analysis Assignment – Pipeline Development

Author: Samuel Perini

With this assignment, I aim to demonstrate my exploratory analysis skills and methodology in three main subjects:

__1. Read coverage.__

__2. Variant calling analysis.__

__3. Analytical performance.__

This repository contains the scripting files that are necessary for running the pipeline but it does not contain the input data that were given for the assignment. The input data should be downloaded and placed in a directory called `Files_needed_for_task` and this directory should be in the same folder as the repository folder (not inside the repository). For example, both `Files_needed_for_task` and `Bioinformatics_pipeline_development_task` (the repository folder name) directories can be placed inside a folder called `Assignment_Pipeline_Development`.

```
Assignment_Pipeline_Development/
├── Bioinformatics_pipeline_development_task
└── Files_needed_for_task
```

The pipeline was built using [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html), a Python-based workflow management system, and after [cloning](https://help.github.com/en/articles/cloning-a-repository) this repository inside the directory `Assignment_Pipeline_Development`, the same bioinformatics steps can be performed by typing the following command line codes:

1. [Install Miniconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)

2. Create an environment called "pipeline_development" with the required software listed inside the YAML file and activate it
```
    cd Bioinformatics_pipeline_development_task
    conda env create --name pipeline_development --file environment.yaml
    conda activate pipeline_development
```

3. Create other environments with the required up-to-date software listed inside different YAML files. The packages of these new environments will be activated inside the pipeline to avoid conflicts with the main environment. Just run the commands below to create them but do not activate these new environments.
```
    conda env create --prefix ./envs/r-envir --file envs/r-environment.yaml
    conda env create --prefix ./envs/var_call_v1 --file ./envs/var-call-env1.yaml
    conda env create --prefix ./envs/var_call_v2 --file ./envs/var-call-env2.yaml
```

4. Download human assembly hg19, chromosome 19 FASTA file and chromosome 19 dbSNP BED file
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

5. Execute the workflow locally printing the commands `-p` and using 1 core (change the number of cores as you wish)
```
    snakemake -p --cores 1
```