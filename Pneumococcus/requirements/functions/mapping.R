###
# Mapping functions to map positions in a given sequence, back to residue numbers in a pdb object.
###

## function map.seqposno.to.resno
map.seqposno.to.resno <- function(sequence_AA1, pdb_object, verbose=FALSE) {
  longest_chain <- pdb_object$atom %>%
    filter(type == "ATOM",
           elety == "CA",
           resid %in% toupper(AMINO_ACID_CODE)) %>%
    add_count(chain, name = "chain_length") %>%
    arrange(desc(chain_length), chain) %>%  # Sort by length (descending), then alphabetically
    filter(chain_length == max(chain_length)) %>%  
    filter(chain == min(chain)) %>% 
    pull(chain) %>%
    unique()
    
  subject <- pdbseq(pdb_object, inds = atom.select(pdb_object, type ="ATOM", elety="CA", chain=longest_chain))
  
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
    filter(type == "ATOM",
           elety == "CA",
           resid %in% toupper(AMINO_ACID_CODE)) %>%
    select(resid,resno,chain) %>%
    full_join(epitope_prediction_df %>%
                mutate(resno = Position + required_adjustment),
              by = join_by
              )
  return(resmap_df)
}
