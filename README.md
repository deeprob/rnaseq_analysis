# RNASeq analysis steps
Steps: | Trim | Align | Count | DE | QC | Enrich |
|---|---|---|---|---|---|---|
Tools: | Trimmomatic | STAR | HTSeq | DESeq2 | sklearn & seaborn | clusterprofiler & enrichplot |

# Steps

- Create metadata file
- Create conda environments
- Run 0_reads_to_counts.py
- Run 1_counts_to_de.py
- Run 2_analyze.py

## Create metadata file

## Create conda environments

### Conda environment for reads to differential expression
```bash
foo@bar:~$ conda create -n rnaseq -c conda-forge -c anaconda -c bioconda python=3 trimmomatic star htseq bioconductor-deseq2 pandas scikit-learn seaborn -y
```

### Conda environment for enrichment analysis 
```bash
foo@bar$ conda create -n gsea -c bioconda -c conda-forge bioconductor-clusterprofiler bioconductor-biomart bioconductor-org.hs.eg.db bioconductor-enrichplot
```

## Run 0_reads_to_counts.py

## Run 1_counts_to_de.py

## Run 2_analyze.py

**Most tools were chosen following the works of Corchete et. al., Systematic comparison and assessment of RNA‑seq procedures for gene expression quantitative analysis**
