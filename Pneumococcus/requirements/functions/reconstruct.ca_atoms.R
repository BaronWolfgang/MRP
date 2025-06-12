reconstruct.ca_atoms <- function(pdb_object, verbose = NULL) {
  
  ca_atoms <- pdb_object$atom %>%
    filter(type == "ATOM",
           elety == "CA",
           alt %in% c(NA, "A"),
           resid %in% toupper(AMINO_ACID_CODE[1:20])) %>%
    add_count(chain, name = "chain_length") %>%
    arrange(desc(chain_length), chain) %>%
    filter(chain_length == max(chain_length)) %>%
    filter(chain == min(chain))  # get main chain
  
  pdb_ca_seq <- pdb_object$atom %>%
    filter(type == "ATOM",
           elety == "CA",
           alt %in% c(NA, "A"),
           resid %in% toupper(AMINO_ACID_CODE)[1:20]) %>%
    add_count(chain, name = "chain_length") %>%
    arrange(desc(chain_length), chain) %>%  # Sort by length (descending), then alphabetically
    filter(chain_length == max(chain_length)) %>%
    filter(chain == min(chain)) %>%
    pull(resid) %>%
    aa321() %>%
    str_c(collapse = "") %>%
    AAString()
  
  if (!is.null(pdb_object$seqres)) {
    pdb_seqres <- pdb_object$seqres %>% enframe(name = "chain", value = "resid") %>%
      add_count(chain, name = "chain_length") %>%
      arrange(desc(chain_length), chain) %>%  # Sort by length (descending), then alphabetically
      filter(chain_length == max(chain_length)) %>%
      filter(chain == min(chain)) %>%
      pull(resid) %>%
      aa321() %>%
      str_c(collapse = "") %>%
      AAString()
  } else {
    pdb_seqres <- pdb_ca_seq
  }
  
  alignment <- pairwiseAlignment(pattern = pdb_seqres, subject = pdb_ca_seq)
  
  # Extract aligned sequences
  aligned_seqres <- as.character(alignedPattern(alignment))
  aligned_ca_seq <- as.character(alignedSubject(alignment))
  
  alignment_df <- tibble(
    seqres_resid = strsplit(aligned_seqres, "")[[1]],
    ca_seq_resid = strsplit(aligned_ca_seq, "")[[1]],
  ) %>%
    mutate(
      seqresno = row_number(),
      ca_index = cumsum(ca_seq_resid != "-"),
      resno = NA_integer_,  # Start with all NAs
      chain = unique(ca_atoms$chain)  # Safe as a scalar
    ) %>%
    mutate(
      resno = replace(
        resno, 
        ca_seq_resid != "-", 
        pull(ca_atoms, resno)[ca_index[ca_seq_resid != "-"]]
      ),
      # Optional: convert codes
      seqres_resid = aa123(seqres_resid),
      ca_seq_resid = aa123(ca_seq_resid)
    ) #%>%
    #select(-ca_index, -ca_seq_resid)
  
  if(!is.null(verbose)) {
    print(alignment)
    print(alignment_df)
    print(ca_atoms)
  }
  
  # Step 4: Create rows for missing residues
  ca_atoms_reconstructed <- left_join(
    alignment_df,
    ca_atoms,
    by = c("chain","resno")
  )
  
  return(ca_atoms_reconstructed)
}