#!/bin/bash
#SBATCH --job-name=convert
#SBATCH --output=logs/alignment_%x_%j.out
#SBATCH --error=logs/alignment_%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=3:00:00
#SBATCH --array=0-130

echo $TMPDIR
cp data/datasets/mouse_mm10_db* ${TMPDIR}
echo "Finished copying file"

TFR_FILE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" tfrecords.txt)
FASTA_FILE=${TFR_FILE%.*}.fasta

blastn -query fasta_seq/${FASTA_FILE} \
       -db ${TMPDIR}/mouse_mm10_db \
       -out alignments/${FASTA_FILE%.*}_blast_results_subset.tsv \
       -evalue 1e-300 \
       -max_target_seqs 1 \
       -num_threads 16 \
       -outfmt '6 qseqid sseqid sstart send evalue bitscore'
