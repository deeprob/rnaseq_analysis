# RNASeq analysis steps
Steps: | Trim | Align | Count | DE
|---|---|---|---|---|
Algorithms: | Trimmomatic | STAR | HTSeq | DESeq2

**Algorithms were chosen following the works of Corchete et. al., Systematic comparison and assessment of RNA‑seq procedures for gene expression quantitative analysis**

# Conda environment for reads to differential expression
```bash
foo@bar:~$ conda create -n rnaseq -c conda-forge -c anaconda -c bioconda python=3 trimmomatic star htseq bioconductor-deseq2 pandas -y
```

# Conda environment for quality control 
```bash
foo@bar:~$ conda create -n rnaseq_qc -c conda-forge -c anaconda -c bioconda python=3 pandas scikit-learn seaborn -y
```

# Conda environment for enrichment analysis 
```bash
foo@bar$ conda create -n gsea -c bioconda -c conda-forge bioconductor-clusterprofiler bioconductor-biomart bioconductor-org.hs.eg.db bioconductor-enrichplot
```
