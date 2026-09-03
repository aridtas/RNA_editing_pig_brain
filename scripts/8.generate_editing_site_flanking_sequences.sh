#!/bin/bash

#SBATCH --time=15:00
#SBATCH -N 1
#SBATCH -c 1
#SBATCH --mem=1G
#SBATCH --output=slurm.output_flanking_%A.txt
#SBATCH --error=slurm.error_flanking_%A.txt
#SBATCH --job-name=flanks

sites_dir=/lustre/nobackup/WUR/ABGC/arida001/results/editing_sites
output_dir=/lustre/nobackup/WUR/ABGC/arida001/results/flanking_seqs
genome_dir=/lustre/nobackup/WUR/ABGC/arida001/mapping/sscrofa_genome

vcf_parameters='g_2_e_005_r_5.vcf'
tissues=('cortex' 'hippocampus' 'hypothalamus')

for tissue in ${tissues[@]}; do

    input=${sites_dir}/${tissue}_${vcf_parameters}

    grep -v "^#" ${input} | awk 'OFS="\t" {
        # $4 = REF, $5 = ALT
        if ($4=="A" && $5=="G") {
            strand = "+";
        } else if ($4=="T" && $5=="C") {
            strand = "-"; }
        
        print $1, $2-1, $2, $3, 0, strand
    }' >> ${output_dir}/stranded_editing_sites.bed
done


bedtools flank -i ${output_dir}/stranded_editing_sites.bed \
-g ${genome_dir}/older_reference_104.fa.fai -l 5 -r 0 -s \
| bedtools getfasta -fi ${genome_dir}/older_reference_104.fa \
-bed stdin -fo ${output_dir}/upstream_editing.fasta -s


bedtools flank -i ${output_dir}/stranded_editing_sites.bed \
-g ${genome_dir}/older_reference_104.fa.fai -l 0 -r 5 -s \
| bedtools getfasta -fi ${genome_dir}/older_reference_104.fa \
-bed stdin -fo ${output_dir}/downstream_editing.fasta -s


