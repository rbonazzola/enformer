#!/bin/bash
#SBATCH --job-name=convert
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=0:10:00
#SBATCH --array=0-130

FILE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" tfrecords.txt)
python extract_tfr_seq_to_fasta.py "$FILE"
