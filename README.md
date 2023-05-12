# RNASeq analysis steps
Steps: | Trim | Align | Count | DE | QC | Enrich |
|---|---|---|---|---|---|---|
Tools: | Trimmomatic | STAR | HTSeq | DESeq2 | sklearn & seaborn | clusterprofiler & enrichplot |

# TODO:
1. Divide rnaseq counts and de completely (glrnacounts and glrnade)
2. No more conda env - use docker envs


# Steps

- Create metadata file
- Create conda environments
- Run 0_reads_to_counts.py
- Run 1_counts_to_de.py
- Run 2_analyze.py

**Most tools were chosen following the works of Corchete et. al., Systematic comparison and assessment of RNA‑seq procedures for gene expression quantitative analysis**
