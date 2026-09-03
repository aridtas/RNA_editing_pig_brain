#!/bin/bash

#SBATCH --time=5:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --job-name=dedup_calm
#SBATCH --array=1-42
#SBATCH --output=slurm.output_e_samples_dedup_and_calmd_%A_%a.txt
#SBATCH --error=slurm.error_e_samples_dedup_and_calmd_%A_%a.txt

sample_list=/lustre/nobackup/WUR/ABGC/arida001/de_novo_rna_edit/compressed_bams/B_and_A_samples_for_preprocessing.txt
out_dir=/lustre/nobackup/WUR/ABGC/arida001/de_novo_rna_edit/compressed_bams
genome=/lustre/nobackup/WUR/ABGC/arida001/mapping/sscrofa_genome/older_reference_104.fa


input_bam_path=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$sample_list")
sname=$(basename ${input_bam_path} .mapped.q30.bam.gz)
temp_bam=${out_dir}/${sname}_temp.bam
dedup_md_bam=${out_dir}/${sname}.dedup.md.bam

# --- Print metadata for the log file ---
echo "--------------------------------"
echo "Job ID: $SLURM_JOB_ID, Task ID: $SLURM_ARRAY_TASK_ID"
echo "Input:    $input_bam_path"
echo "Output:   $dedup_md_bam"
echo "--------------------------------"

gunzip -c "$input_bam_path" | \
samtools sort -n - | \
samtools fixmate -m - - | \
samtools sort - | \
samtools markdup -r - "$temp_bam"

echo "Deduplication finished, initiating MD tag calculation"

samtools calmd -b ${temp_bam} ${genome} > ${dedup_md_bam}

rm -f ${temp_bam}

samtools index ${dedup_md_bam}

echo "MD calculation finished"


#read_header=$(samtools view "$temp_bam" | head -n 1 | cut -f 1)
#flowcell=$(echo "$read_header" | cut -d ':' -f 3)
#lane=$(echo "$read_header" | cut -d ':' -f 4)
#rgpu="${flowcell}.${lane}"


#picard AddOrReplaceReadGroups -I ${temp_bam} \
#	   -O ${dedup_bam} \
#	   --RGLB ${sname} \
#	   --RGPL ILLUMINA \
#	   --RGPU ${rgpu} \
#	   --RGSM ${sname} \
#	   --RGID ${sname}

#alternative setup without a sample list
#FILES_TO_PROCESS=(${bam_dir}/*.bam)
#INPUT_FILE=${FILES_TO_PROCESS[$SLURM_ARRAY_TASK_ID]}
#BASENAME=$(basename "$INPUT_FILE" .q30.bam)

#other way to define array: array=0-$(($(wc -l < "$sample_list") - 1)) careful with 0-index and 1-index
