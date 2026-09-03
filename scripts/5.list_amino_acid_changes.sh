#!/bin/bash

FILE=$1
OUTPUT="aa_changes_hypothalamus.tsv"

grep -v "^#" "$FILE" | awk '
BEGIN {
    FS="\t"; OFS="\t";

    # Severity Hierarchy (Rank 1 = Most Severe)
    rank["stop_gained"] = 1;
    rank["stop_lost"] = 2;
    rank["start_lost"] = 3;
    rank["missense_variant"] = 4;
    rank["splice_region_variant"] = 5;
    rank["synonymous_variant"] = 6;
}
$9 != "-" {  # CDS only

    variant_id = $1;   # First column
    loc = $2;          # Genomic location
    cons = $7;         # Consequence(s)
    aa_change = $11;   # Amino acid change

    # --- Determine most severe consequence for this line ---
    split(cons, arr, ",");
    best_rank = 100;
    best_cons = "unknown";

    for (i in arr) {
        c = arr[i];
        r = (c in rank) ? rank[c] : 100;
        if (r < best_rank) {
            best_rank = r;
            best_cons = c;
        }
    }

    # --- Store/Update best per site ---
    # We use the location (loc) as a key to ensure we only report the most 
    # severe consequence if a variant appears on multiple lines.
    if (!(loc in stored_rank) || best_rank < stored_rank[loc]) {
        stored_rank[loc] = best_rank;
        stored_id[loc] = variant_id;
        
        if (best_cons == "missense_variant" && aa_change != "-") {
            stored_aa[loc] = aa_change;
        } else {
            stored_aa[loc] = "NA";
        }
    }
}
END {
    # --- Print individual associations instead of counts ---
    for (l in stored_id) {
        if (stored_aa[l] != "NA") {
            print stored_id[l], stored_aa[l];
        }
    }
}' >> "$OUTPUT"

echo "Done! AA changes saved to $OUTPUT"