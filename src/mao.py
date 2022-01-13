import os
import subprocess

# create a metadict with the primary sample and fwd reads and rev reads for seven reps
meta_samples = ("Homo", "Hetz", "WT")
meta_dict = {k:[] for k in meta_samples}

root_dir = "/data5/deepro/mao_lab/embryonic_brains_gyf2_ko/data"
sample_reps = [f.name for f in os.scandir(root_dir) if os.path.isdir(f.path)]

for rep in sample_reps:
    sample = rep.split("-")[0]
    meta_dict[sample].append(rep)

# sort metadict
meta_dict = {k:sorted(v) for k,v in meta_dict.items()}

# for each sample, generate fwd and rev read files, then run the rnaseq pipeline
ROOT_DIR="/data5/deepro/mao_lab/embryonic_brains_gyf2_ko/results"
GENOME="/data5/deepro/genomes/mm39/mm39.fa"
GTF="/data5/deepro/genomes/mm39/refGene.gtf"


def run_rnaseq(sample_name, fwd_samples, rev_samples):
    command = ["python", "run.py", ROOT_DIR, sample_name]
    command += fwd_samples
    if rev_samples:
        command += ["-r"]
        command += rev_samples
    command += ["-g", GENOME, "-e", GTF]
    print(command)
    subprocess.run(command)
    return

for sample, sample_reps in meta_dict.items():
    sample_rep_fwd = [os.path.join(root_dir, sr, sr+"_1.fq.gz") for sr in sample_reps]
    sample_rep_rev = [os.path.join(root_dir, sr, sr+"_2.fq.gz") for sr in sample_reps]

    run_rnaseq(sample, sample_rep_fwd, sample_rep_rev)





