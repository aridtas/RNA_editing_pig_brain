#!/bin/bash

#SBATCH --time=30:00
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=8G
#SBATCH --output=slurm.output_fastqc_%A_%a.txt
#SBATCH --error=slurm.error_fastqc_%A_%a.txt
#SBATCH --job-name=fastqc_parallel
#SBATCH --array=1-3


lustre_dir=/lustre/nobackup/WUR/ABGC/arida001
output_dir=${lustre_dir}/quality_control/fastqc_output/
input_dir=${lustre_dir}/raw_reads/
sample_file=${input_dir}/samples_ID.txt


sname=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$sample_file")


mkdir -p ${output_dir}/${sname}

sample_dir=${input_dir}/$sname

read_files=( $(find ${sample_dir} -name "*.fq.gz") )

eval "$(/home/WUR/arida001/miniforge3/bin/mamba shell hook --shell bash)"

mamba activate qualcon

fastqc -t 2 -o ${output_dir}/${name} ${read_files[@]}

sacct --format="CPUTime,MaxRSS"

mv slurm.output_fastqc_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} ${output_dir}/${sname}
mv slurm.error_fastqc_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID} ${output_dir}/${sname}

echo "successfully done"
