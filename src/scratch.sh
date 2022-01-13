#!/bin/bash
set -ue


# File to run rnaseq for starrseq data
ROOT_DIR="/data5/deepro/starrseq/rnaseq"
REP_DIR="all"
READ1="/data5/deepro/starrseq/rnaseq/raw/all/16P12_1.fastq.gz"
READ2=""
GENOME="/data5/deepro/genomes/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta"
GTF="/data5/deepro/genomes/gencode.v38.annotation.gtf"

if [ -z "$READ2" ]
then 
    time python run.py ${ROOT_DIR} ${REP_DIR} ${READ1} "" ${GENOME} ${GTF}
else 
    time python run.py ${ROOT_DIR} ${REP_DIR} ${READ1} ${READ2} ${GENOME} ${GTF}
fi
