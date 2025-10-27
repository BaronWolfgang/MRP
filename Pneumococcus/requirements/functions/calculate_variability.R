calculate_variability <- function(aligned_aa_seqs) {
  alignment_matrix <- as.matrix(aligned_aa_seqs)
  
  # Replace terminal gaps with NA per sequence (row)
  replace_terminal_gaps_with_NA <- function(row) {
    non_gap_idx <- which(row != "-")
    if (length(non_gap_idx) == 0) {
      # sequence is entirely gaps
      row[] <- NA_character_
      return(row)
    }
    start <- min(non_gap_idx)
    end   <- max(non_gap_idx)
    if (start > 1) row[1:(start - 1)] <- NA_character_
    if (end < length(row)) row[(end + 1):length(row)] <- NA_character_
    row
  }
  
  for (i in seq_len(nrow(alignment_matrix))) {
    alignment_matrix[i, ] <- replace_terminal_gaps_with_NA(alignment_matrix[i, ])
  }
  
  # Determine consensus sequence manually
  consensus_seq <- unlist(
    apply(alignment_matrix, 2, function(col) {
      residues <- col[!is.na(col)]
      if (length(residues) == 0) {
        NA_character_
      } else {
        names(sort(table(residues), decreasing = TRUE))[1]
      }
    }),
    use.names = FALSE
  )
  
  coverage <- colSums(!is.na(alignment_matrix))
  
  # couning srps
  srp_counts <- sapply(seq_len(ncol(alignment_matrix)), function(i) {
    sum(
      alignment_matrix[, i] != consensus_seq[i] &
          alignment_matrix[, i] != "-" & !is.na(alignment_matrix[, i]) &
          consensus_seq[i] != "-" & !is.na(consensus_seq[i])
      )
  })
  
  srp_frequencies <- srp_counts / coverage
  
  # counting gaps
  gap_counts <- sapply(seq_len(ncol(alignment_matrix)), function(i) {
    sum(alignment_matrix[, i] == "-", na.rm = TRUE)
  })
  
  #NOTE: adjusting for gaps occuring in all 100% of sequences when mapping back to PDB that contains histags
  gap_counts <- ifelse(gap_counts == coverage, 0, gap_counts)
    
  gap_freq <- gap_counts / coverage
  
  # Determine Indel type relative to consensus sequence
  indel_type <- case_when(
    gap_counts > 0 & consensus_seq == "-" ~ "insertion",
    gap_counts > 0 & consensus_seq != "-" ~ "deletion",
    TRUE ~ NA_character_
  )
  
  # Residue table per column
  residue_tables <- apply(alignment_matrix, 2, function(col) {
    sorted_table <- sort(table(col, useNA = "no"), decreasing = TRUE)
    paste0(names(sorted_table), ":", as.integer(sorted_table), collapse = "; ")
  })
  
  sequence_length <- length(consensus_seq)
  n_seq <- length(aligned_aa_seqs)
  
  data_variability <- data.frame(
    position = seq_len(sequence_length),
    consensus_res = consensus_seq,
    residue_tables = residue_tables,
    srp_counts = srp_counts,
    srp_freq = srp_frequencies,
    gap_counts = gap_counts,
    gap_freq = gap_freq,
    indel_type = indel_type,
    n_seq = n_seq,
    coverage = coverage,
    seq_length = sequence_length
  ) %>%
    dplyr::mutate(
      indel_counts = dplyr::case_when(
        indel_type == "deletion" ~ gap_counts, 
        indel_type == "insertion" ~ coverage - gap_counts, 
        TRUE ~ 0
      ),
      indel_freq = indel_counts / coverage,
      var_counts = srp_counts + indel_counts,
      var_freq = var_counts / coverage,
    )
  
  return(data_variability)
}
