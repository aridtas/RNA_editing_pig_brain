#!/bin/bash

#SBATCH --time=5:00
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH --output=slurm.output_vep_%A.txt
#SBATCH --error=slurm.error_vep_%A.txt
#SBATCH --job-name=vep

tool_dir=/home/WUR/arida001/ensembl-vep
genome_dir=/lustre/nobackup/WUR/ABGC/arida001/mapping/sscrofa_genome


${tool_dir}/vep -i $1 -o $2 --stats_file $3 --pick --offline --cache --species sus_scrofa
