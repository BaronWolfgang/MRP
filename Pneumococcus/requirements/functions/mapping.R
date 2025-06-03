###
# Mapping functions to map positions in a given sequence by finding the exact match between the the structure-based and sequence-based predictions, back to residue numbers in a pdb object.
###

## function map.seqposno.to.resno
map.seqposno.to.resno <- function(sequence_AA1, pdb_object, verbose=FALSE) {
  longest_chain <- pdb_object$atom %>%
    filter(type == "ATOM",
           elety == "CA",
           alt %in% c(NA, "A"),
           resid %in% toupper(AMINO_ACID_CODE)) %>%
    add_count(chain, name = "chain_length") %>%
    arrange(desc(chain_length), chain) %>%  # Sort by length (descending), then alphabetically
    filter(chain_length == max(chain_length)) %>%  
    filter(chain == min(chain)) %>% 
    pull(chain) %>%
    unique()
    
  ca_atoms <- pdb_object$atom %>%
    filter(
      type == "ATOM",
      elety == "CA",
      chain == longest_chain,
      alt %in% c(NA, "A"),
      resid %in% toupper(AMINO_ACID_CODE)
    )
  
  subject <- aa321(ca_atoms$resid)
  names(subject) <- ca_atoms$resno
  
  subject_resno <- as.numeric(names(subject))
  subject_seq <- str_c(
    unname(subject),
    collapse = ""
  )
  
  pattern_seq <- sequence_AA1
  
  #custom substitution matrix to find the longest continuous substring (LCS) without gaps or mismatches
  customsubmatrix <- matrix(-Inf, 
                            nrow = 20, 
                            ncol = 20, 
                            dimnames = list(AA_ALPHABET[1:20],AA_ALPHABET[1:20]))
  diag(customsubmatrix) <- 1
  
  #align to find the LCS
  alignment <- pairwiseAlignment(pattern_seq, subject_seq, type="local", gapOpening = Inf,
                                 substitutionMatrix = customsubmatrix)
  
  #extract the LCS
  aligned_pattern <- as.character(alignment)
  
  #locate the LCS in both strings
  loc_subject <- str_locate(subject_seq,pattern = aligned_pattern)
  loc_pattern <- str_locate(pattern_seq,pattern = aligned_pattern)
  
  #Find the difference in the position of the LCS
  shift <- loc_subject[1] - loc_pattern[1] 
  
  #Look at the first resno number in the pdb file, add the calculated shift
  required_adjustment <- subject_resno[1] - 1 + shift
  
  #verbose
  if(verbose) {
    print(aligned_pattern)
    print(loc_subject)
    print(alignedSubject(alignment)) 
    print(loc_pattern)
    print(alignedPattern(alignment)) 
    print(shift)
  }
  
  return(c(required_adjustment))
}

## map.seqresno.to.resno
library("Biostrings")

map.seqresno.to.resno <- function(pdb_object, verbose=FALSE) {
  subject_seq <- pdb_object$atom %>%
    filter(type == "ATOM",
           elety == "CA",
           alt %in% c(NA, "A"),
           resid %in% toupper(AMINO_ACID_CODE)) %>%
    select(resno, resid) %>%
    distinct(.keep_all = TRUE) %>%
    mutate(resid = aa321(resid)) %>%
    pull(resid) %>%
    str_c(collapse = "")
  
  required_adjustment <- map.seqposno.to.resno(subject_seq, pdb_object, verbose=verbose)
  
  return(c(required_adjustment))
}

# create resmap
create.resmap <- function(pdb_atom_df, epitope_prediction_df,required_adjustment, join_by = c("resno")) {
  resmap_df <- pdb_atom_df %>% 
    filter(!is.na(chain),
           type == "ATOM",
           elety == "CA",
           alt %in% c(NA, "A"),
           resid %in% toupper(AMINO_ACID_CODE)) %>%
    #select(resid,resno,chain,x,y,z,o,b) %>%
    full_join(epitope_prediction_df %>%
                mutate(resno = Position + required_adjustment),
              by = join_by
              )
  return(resmap_df)
}

########
#For mapping variability to our epitope prediction data, we need a slightly different method as the consensus sequence of the variability data isn't guaranteed to match the prediction data
########
map.seq.to.pdb <- function(sequence_AA1, pdb_object, alignment_method = "overlap") {
  longest_chain <- pdb_object$atom %>%
    filter(type == "ATOM",
           elety == "CA",
           alt %in% c(NA, "A"),
           resid %in% toupper(AMINO_ACID_CODE)) %>%
    add_count(chain, name = "chain_length") %>%
    arrange(desc(chain_length), chain) %>%
    filter(chain_length == max(chain_length)) %>%
    filter(chain == min(chain)) %>%
    pull(chain) %>%
    unique()
  
  ca_atoms <- pdb_object$atom %>%
    filter(
      type == "ATOM",
      elety == "CA",
      chain == longest_chain,
      alt %in% c(NA, "A"),
      resid %in% toupper(AMINO_ACID_CODE)
    )
  
  subject <- aa321(ca_atoms$resid)
  names(subject) <- ca_atoms$resno
  
  subject_resno <- as.numeric(names(subject))
  subject_seq <- str_c(
    unname(subject),
    collapse = ""
  )
  
  pattern_seq <- sequence_AA1
  
  alignment <- pairwiseAlignment(pattern_seq, subject_seq,
                                 gapOpening = 5,
                                 gapExtension = 0,
                                 type = alignment_method)  # Or "global-local" depending on use case
  
  return(alignment)
}

create.resmap.variability <- function(pdb_atoms_df, variability_df, alignment) {
  aligned_pattern <- as.character(alignedPattern(alignment))
  aligned_subject <- as.character(alignedSubject(alignment))
  
  pattern_chars <- strsplit(aligned_pattern, "")[[1]]
  subject_chars <- strsplit(aligned_subject, "")[[1]]
  
  stopifnot(length(pattern_chars) == length(subject_chars))
  
  resnos_aligned <- vector("numeric", length(subject_chars))
  pdb_index <- 1
  
  for (i in seq_along(subject_chars)) {
    if (subject_chars[i] == "-") {
      resnos_aligned[i] <- NA
    } else {
      resnos_aligned[i] <- ca_atoms$resno[pdb_index]
      pdb_index <- pdb_index + 1
    }
  }
  
  subject_chars[subject_chars == "-"] <- NA
  
  pattern_start <- start(ranges(alignment@pattern))
  subject_start <- start(ranges(alignment@subject))
  
  alignment_df <- data.frame(
    resno = resnos_aligned,
    consensus_res = pattern_chars,
    resid = subject_chars,
    stringsAsFactors = FALSE
  ) %>%
    mutate(
      position_new = pattern_start + row_number() - subject_start,
      position = pattern_start + cumsum(consensus_res != "-") - subject_start
      )
  
  resmap <- full_join(
    alignment_df,
    variability_df, 
    by = c("position")
    #by = c("position")
    )
  
  return(resmap)
  
}