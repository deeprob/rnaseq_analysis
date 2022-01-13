import subprocess


def create_star_index(genome_file, gtf_file, storage_dir, nthreads=64):
    # create index for star aligner
    subprocess.call(["STAR", "--runMode", "genomeGenerate",  "--genomeDir", f"{storage_dir}",
                    "--genomeFastaFiles", f"{genome_file}", "--sjdbGTFfile", f"{gtf_file}", "--sjdbOverhang", "100", 
                    "--runThreadN", f"{nthreads}"])
    return   


def align(read_files, output_prefix, starindex_dir, nthreads=64):
    # align using star
    read_files_list = read_files.split()
    command = ["STAR", "--genomeDir", f"{starindex_dir}", "--runThreadN", f"{nthreads}", "--outFileNamePrefix", f"{output_prefix}", "--readFilesIn"]
    command += read_files_list
    command += ["--outSAMtype", "BAM", "SortedByCoordinate", "--outSAMunmapped", "Within", "--outSAMattributes", "Standard", "--readFilesCommand", "zcat"]
    subprocess.call(command)
    return
