#!/bin/bash

#SBATCH --time=6:00:00
#SBATCH -N 1
#SBATCH -c 16
#SBATCH --mem=70G
#SBATCH --output=slurm.output_STAR_indices%A.txt
#SBATCH --error=slurm.error_STAR_indices%A.txt
#SBATCH --job-name=genome_indices


#### Time and memory requests were overkill, took about 15 minutes and less than 5 GB

genome_dir=/lustre/nobackup/WUR/ABGC/arida001/mapping/sscrofa_genome
indices_dir=${genome_dir}/genome_indices

eval "$(/home/WUR/arida001/miniforge3/bin/mamba shell hook --shell bash)"

mamba activate Mapping


STAR --runThreadN 16 --runMode genomeGenerate --genomeDir ${indices_dir} \
     --genomeFastaFiles ${genome_dir}/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa \
     --sjdbGTFfile ${genome_dir}/Sus_scrofa.Sscrofa11.1.115.gtf \
     --sjdbOverhang 149


sacct --format="CPUTime,MaxRSS"

echo "successfully done"
