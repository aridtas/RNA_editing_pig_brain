#!/bin/bash

OUTPUT="bam_read_counts.tsv"

# Create header
echo -e "Sample_Name\tRead_Count" > "$OUTPUT"

for i in *.bam; do
    # Calculate sum of reads from mapped column ($3)
    count=$(samtools idxstats "${i}" | awk '{sum+=$3} END {print sum}')
    
    # Output sample name and count separated by a tab
    echo -e "${i%.dedup.md.bam}\t${count}" >> "$OUTPUT"
done
