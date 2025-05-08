calculate_variability <- function(aligned_aa_seqs) {
  alignment_matrix <- as.matrix(aligned_aa_seqs)
  
  # Determine consensus sequence manually
  consensus_seq <- apply(alignment_matrix, 2, function(col) {
    residues <- col
    unique_residues <- unique(col)
    most_frequent <- names(sort(table(residues), decreasing = TRUE))[1]
  })
  
  # Variability ignoring gaps
  variability_counts <- sapply(seq_len(ncol(alignment_matrix)), function(i) {
    sum(alignment_matrix[, i] != consensus_seq[i] &
          alignment_matrix[, i] != "-" &
          consensus_seq[i] != "-")
  })
  
  variability_frequencies <- variability_counts / length(aligned_aa_seqs)
  
  # counting gaps
  gap_counts <- sapply(seq_len(ncol(alignment_matrix)), function(i) {
    sum(alignment_matrix[, i] == "-")
  })
  
  gap_freq <- gap_counts / length(aligned_aa_seqs)
  
  # Deermine Indel type relative to consensus sequence
  indel_type <- ifelse(consensus_seq == "-", "insertion", "deletion")
  
  # Residue table per column
  residue_tables <- apply(alignment_matrix, 2, function(col) {
    sorted_table <- sort(table(col), decreasing = TRUE)
    paste0(names(sorted_table), ":", as.integer(sorted_table), collapse = "; ")
  })
  
  sequence_length <- length(consensus_seq)
  n_seq <- length(aligned_aa_seqs)
  
  data_variability <- data.frame(
    position = seq_len(sequence_length),
    consensus_res = consensus_seq,
    residue_tables = residue_tables,
    var_counts = variability_counts,
    var_freq = variability_frequencies,
    gap_counts = gap_counts,
    gap_freq = gap_freq,
    indel_type = indel_type,
    n_seq = n_seq,
    seq_length = sequence_length
  ) %>%
    dplyr::mutate(
      indel_counts = dplyr::case_when(
        indel_type == "deletion" ~ gap_counts, 
        indel_type == "insertion" ~ n_seq - gap_counts, 
        TRUE ~ 0
      ),
      indel_freq = indel_counts / n_seq
    )
  
  return(data_variability)
}
