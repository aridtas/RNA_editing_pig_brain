#!/bin/bash

#SBATCH --time=20:00
#SBATCH -N 1
#SBATCH -c 12
#SBATCH --mem=35Gallaksetarreeeeeee
#SBATCH --output=slurm.output_STAR_mapping_%A_%a.txt
#SBATCH --error=slurm.error_STAR_mapping_%A_%a.txt
#SBATCH --job-name=mapping
#SBATCH --array=1-3

reads_dir=/lustre/nobackup/WUR/ABGC/arida001/raw_reads
sample_file=${reads_dir}/samples_ID.txt
genome_indices_dir=/lustre/nobackup/WUR/ABGC/arida001/mapping/sscrofa_genome/genome_indices/
output_dir=/lustre/nobackup/WUR/ABGC/arida001/mapping

sname=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$sample_file")

read_files=( $(find ${reads_dir}/${sname} -name "*.fq.gz" | sort) )

mkdir -p "${output_dir}/${sname}"

eval "$(/home/WUR/arida001/miniforge3/bin/mamba shell hook --shell bash)"

mamba activate Mapping

STAR --runThreadN 12 \
     --genomeDir "${genome_indices_dir}" \
     --readFilesIn "${read_files[0]}" "${read_files[1]}" \
     --readFilesCommand zcat \
     --outSAMtype BAM Unsorted \
     --outReadsUnmapped Fastx \
     --outFileNamePrefix "$output_dir/${sname}/${sname}." \
     --quantMode GeneCounts


sacct --format="CPUTime,MaxRSS"

echo "process finished"
