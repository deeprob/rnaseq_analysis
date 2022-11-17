# RNASeq analysis steps
Steps: | Trim | Align | Count | DE
|---|---|---|---|---|
Algorithms: | Trimmomatic | STAR | HTSeq | DESeq2

**Algorithms were chosen following the works of Corchete et. al., Systematic comparison and assessment of RNA‑seq procedures for gene expression quantitative analysis**

# Conda environment for reads to differential expression
```bash
foo@bar:~$ conda create -n rnaseq -c conda-forge -c anaconda -c bioconda python=3 trimmomatic star htseq bioconductor-deseq2 pandas bioconductor-clusterprofiler==3.6.0 bioconductor-biomart==2.34.2 bioconductor-org.hs.eg.db bioinfokit scikit-learn seaborn -y
```

# Conda environment for quality control 
```bash
foo@bar:~$ conda create -n rnaseq_qc -c conda-forge -c anaconda -c bioconda python=3 pandas bioinfokit scikit-learn seaborn -y
```

# Conda environment for enrichment analysis 
```bash
foo@bar$ conda create -n gsea -c bioconda -c conda-forge bioconductor-clusterprofiler==3.6.0 bioconductor-biomart==2.34.2 bioconductor-org.hs.eg.db
```

# Conda environment for enrichment visualization 
```bash
foo@bar$ conda create -n gseaviz -c bioconda -c conda-forge bioconductor-enrichplot
```