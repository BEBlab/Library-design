
library(seqinr)
library(writexl)

set.seed(42)

# Wild-type DNA sequence: change accordingly
wt <- s2c("ATGGCCTCAAACGATTATACCCAACAAGCAACCCAAAGCTATGGGGCCTACCCCACCCAGCCCGGGCAGGGCTATTCCCAGCAGAGCAGTCAGCCCTACGGACAGCAGAGTTACAGTGGTTATAGCCAGTCC")
wt_codons <- splitseq(wt, frame = 0, word = 3)
num_codons <- length(wt_codons)

# Synonymous codon table
syn_dict <- list(
  A = c("GCT", "GCC", "GCA", "GCG"),
  R = c("CGT", "CGC", "CGA", "CGG", "AGA", "AGG"),
  N = c("AAT", "AAC"),
  D = c("GAT", "GAC"),
  C = c("TGT", "TGC"),
  Q = c("CAA", "CAG"),
  E = c("GAA", "GAG"),
  G = c("GGT", "GGC", "GGA", "GGG"),
  H = c("CAT", "CAC"),
  I = c("ATT", "ATC", "ATA"),
  L = c("TTA", "TTG", "CTT", "CTC", "CTA", "CTG"),
  K = c("AAA", "AAG"),
  M = c("ATG"),
  F = c("TTT", "TTC"),
  P = c("CCT", "CCC", "CCA", "CCG"),
  S = c("TCT", "TCC", "TCA", "TCG", "AGT", "AGC"),
  T = c("ACT", "ACC", "ACA", "ACG"),
  W = c("TGG"),
  Y = c("TAT", "TAC"),
  V = c("GTT", "GTC", "GTA", "GTG"),
  STOP = c("TAA", "TAG", "TGA")
)

# Homology regions: change accordingly
constant_region_before <- "CCTCTATACTTTAACGTCAAGGAG"
constant_region_after <- "ACGGACACTTCAGGCTATGGCCAG"

# Prepare output list
results <- list()

# Loop over NNK positions: add the number of codons that you want to mutagenise
for (nnk_pos in 1:44) {
  codons <- wt_codons
  codons[nnk_pos] <- "NNK"  # Insert degenerate NNK
  
  # Pick a random position ≠ NNK pos
  possible_syn_pos <- setdiff(1:num_codons, nnk_pos)
  syn_pos <- sample(possible_syn_pos, 1)
  
  original_codon <- wt_codons[syn_pos]
  aa <- translate(s2c(original_codon))
  
  # Check for synonymous options
  syn_codons <- syn_dict[[aa]]
  syn_codons <- setdiff(syn_codons, original_codon)
  
  # Apply if possible
  if (length(syn_codons) > 0) {
    new_codon <- sample(syn_codons, 1)
    codons[syn_pos] <- new_codon
  } else {
    syn_pos <- 0
    new_codon <- ""
    original_codon <- ""
  }
  
  # Assemble final sequence
  final_seq <- paste0(constant_region_before, paste(codons, collapse = ""), constant_region_after)
  
  results[[nnk_pos]] <- data.frame(
    Full_Sequence = final_seq,
    NNK_Position = nnk_pos,
    Syn_Position = syn_pos,
    Syn_From = original_codon,
    Syn_To = new_codon,
    stringsAsFactors = FALSE
  )
}

# Export to Excel
final_df <- do.call(rbind, results)
#Change filename accordingly
write_xlsx(final_df, "FUS1_NNK_with_random_synonymous.xlsx")

