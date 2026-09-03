#!/bin/bash

#SBATCH --time=10:00
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=6G
#SBATCH --output=slurm.output_multiqc_%A.txt
#SBATCH --error=slurm.error_multiqc_%A.txt
#SBATCH --job-name=multiqc_after

fastqc_dir=/lustre/nobackup/WUR/ABGC/arida001/quality_control/fastqc_output

eval "$(/home/WUR/arida001/miniforge3/bin/mamba shell hook --shell bash)"

mamba activate qualcon

multiqc ${fastqc_dir} -o ${fastqc_dir}

mv slurm.output_multiqc_${SLURM_JOB_ID}.txt $fastqc_dir
mv slurm.error_multiqc_${SLURM_JOB_ID}.txt $fastqc_dir

