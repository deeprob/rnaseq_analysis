# RNASeq analysis steps
Steps: | Trim | Align | Count | DE
|---|---|---|---|---|
Algorithms: | Trimmomatic | STAR | HTSeq | DESeq2

**Algorithms were chosen following the works of Corchete et. al., Systematic comparison and assessment of RNA‑seq procedures for gene expression quantitative analysis**

# Conda environment
```bash
foo@bar:~$ conda create -n rnaseq -c conda-forge -c anaconda -c bioconda python=3 trimmomatic star htseq bioconductor-deseq2 -y
```
