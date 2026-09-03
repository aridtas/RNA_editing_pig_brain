#!/bin/bash

#SBATCH --time=45:00
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH --output=slurm.output_random_flanking_%A.txt
#SBATCH --error=slurm.error_random_flanking_%A.txt
#SBATCH --job-name=rflank

# --- Configuration ---
output_dir="/lustre/nobackup/WUR/ABGC/arida001/results/flanking_seqs"
genome_dir="/lustre/nobackup/WUR/ABGC/arida001/mapping/sscrofa_genome"
genome_fasta="${genome_dir}/older_reference_104.fa"
genome_index="${genome_dir}/older_reference_104.fa.fai"

upstream_flank=5
downstream_flank=5
target_count=29840

# --- Intermediate Files ---
candidates_bed="${output_dir}/candidates.bed"
validated_sites="${output_dir}/random_A_sites.bed"
final_fasta="${output_dir}/random_control_flanking.fasta"


bedtools random -l 1 -n 200000 -g "${genome_index}" > "${candidates_bed}"


bedtools getfasta -fi "${genome_fasta}" -bed "${candidates_bed}" -tab \
| awk -v limit="$target_count" '
    BEGIN {
        count=0; 
        OFS="\t";  
    } 
    {
        if ($2 == "A" || $2 == "a") {
            
            # Split coordinate string "chr1:100-101" back into columns
            split($1, coords, "[:-]"); 
            chrom = coords[1];
            start = coords[2];
            end = coords[3];
            
            print chrom, start, end, "random_A_"count, "0", "+"
            
            count++
            if (count == limit) exit
        }
    }' > "${validated_sites}"

# Check counts
actual_count=$(wc -l < "${validated_sites}")

echo "   Found ${actual_count} sites matching 'A'."


bedtools flank \
    -i "${validated_sites}" \
    -g "${genome_index}" \
    -l $upstream_flank \
    -r $downstream_flank \
    -s \
| bedtools getfasta \
    -fi "${genome_fasta}" \
    -bed stdin \
    -fo "${final_fasta}" \
    -s

echo "Done! Control sequences saved to: ${final_fasta}"