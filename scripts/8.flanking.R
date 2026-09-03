library(ggplot2)
library(ggseqlogo)
library(Biostrings)



fasta_up <- readDNAStringSet("upstream_editing.fasta")
fasta_down <- readDNAStringSet("downstream_editing.fasta")
seqs_up <- as.character(fasta_up)
seqs_down <- as.character(fasta_down)


p_up <- ggseqlogo(seqs_up, method = "prob") +
  ggtitle("Upstream sequence context of editing sites") +
scale_x_continuous(breaks = 1:5, labels = c("-5", "-4", "-3", "-2", "-1")) +
  xlab("Distance from Editing Site (bp)")

p_down <- ggseqlogo(seqs_down, method = "prob") +
  ggtitle("Downstream sequence context of editing sites")


control_fasta <- readDNAStringSet("random_control_flanking.fasta")
control_seqs <- as.character(control_fasta)
control_upstream   <- control_seqs[seq_along(control_seqs) %% 2 == 1]
control_downstream <- control_seqs[seq_along(control_seqs) %% 2 == 0]

control_p_up <- ggseqlogo(control_upstream, method = "prob") +
  ggtitle("Upstream sequence context of control sites") +
  scale_x_continuous(breaks = 1:5, labels = c("-5", "-4", "-3", "-2", "-1")) +
  xlab("Distance from A Site (bp)")

control_p_down <- ggseqlogo(control_downstream, method = "prob") +
  ggtitle("Downstream sequence context of control sites")



up_seqs   <- as.character(readDNAStringSet("upstream_editing.fasta"))
down_seqs <- as.character(readDNAStringSet("downstream_editing.fasta"))


full_seqs <- paste0(up_seqs,'A', down_seqs)
control_full_seqs <- paste0(control_upstream, "A", control_downstream)


ggseqlogo(control_full_seqs, method = "prob") +
scale_x_continuous(
      breaks = 1:11,
      labels = c(
        "-5", "-4", "-3", "-2", "-1",
        "0",
        "+1", "+2", "+3", "+4", "+5"
      )
    ) +
    
xlab("Distance from Control Adenine (bp)") +
ylab("Frequency") +
ggtitle("Sequence Context Flanking Control Adenines") +
annotate("text", x = 6, y = 1.05, label = "A", vjust = 0)
  

  