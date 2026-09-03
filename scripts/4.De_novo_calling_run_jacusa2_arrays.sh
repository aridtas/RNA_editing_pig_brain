#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH -N 1
#SBATCH -c 2
#SBATCH --mem=10G
#SBATCH --output=slurm.output_hypothalamus_jacusa_arrays_%A_%a.txt
#SBATCH --error=slurm.error_hypothalamus_jacusa_arrays_%A_%a.txt
#SBATCH --job-name=jac_hypo_arr
#SBATCH --array=1-7


output_dir=/lustre/nobackup/WUR/ABGC/arida001/de_novo_rna_edit/jacusa_output_files
bam_dir=/lustre/nobackup/WUR/ABGC/arida001/de_novo_rna_edit/compressed_bams
tool_dir=/home/WUR/arida001
snp_dir=/lustre/nobackup/WUR/ABGC/arida001/de_novo_rna_edit/snp_files/split_by_chromosome_beds
sample_list=/lustre/nobackup/WUR/ABGC/arida001/de_novo_rna_edit/compressed_bams/B_and_A_samples_comp.txt
region_dir=/lustre/nobackup/WUR/ABGC/arida001/de_novo_rna_edit/snp_files/chromosome_coords


chrom=${SLURM_ARRAY_TASK_ID}


if [ -z "$chrom" ]; then
    echo "Error: No chromosome found for Array Task ID $SLURM_ARRAY_TASK_ID"
    exit 1
fi

echo "Processing Array Task: $SLURM_ARRAY_TASK_ID"
echo "Target Chromosome: $chrom"


bam_list=$(tr '\n' ',' < "$sample_list")


java -Xmx8g -jar ${tool_dir}/jacusa2 call-1 \
	 -r ${output_dir}/hypothalamus_calls_chrom_${chrom}.out \
	 -a B,M,I:distance=10,S,Y,E:file=${snp_dir}/pig_SNPs_chrom_${chrom}.bed:type=BED \
	 -b ${region_dir}/chromosome_${chrom}.bed \
	 -f B \
	 -P RF-FIRSTSTRAND \
	 -p 2 \
	 -s \
	 -m 20 \
	 -q 24 \
	 -c 1 \
	 ${bam_list}
