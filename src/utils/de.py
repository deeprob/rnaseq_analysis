import os
import subprocess

deseq2_script_path = os.path.join(os.path.dirname(__file__), "deseq2.R")


def run_deseq2(counts_file, counts_cols, designfile, design_formula, contrast, de_file):
	counts_cols_to_str = ",".join(counts_cols)
	cmd = ["Rscript", deseq2_script_path, counts_file, counts_cols_to_str, designfile, design_formula, contrast, de_file]
	subprocess.run(cmd)
	return
