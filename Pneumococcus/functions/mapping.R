###
# Mapping functions to map positions in a given sequence, back to residue numbers in a pdb object.
###

## function map.seqposno.to.resno
map.seqposno.to.resno <- function(sequence_AA1, pdb_object, verbose=FALSE, type="overlap") {
  pattern <- pdbseq(pdb_object)
  unique_pattern <- pattern[!duplicated(names(pattern))]
  pattern_resno <- as.numeric(names(unique_pattern))
  pattern_seq <- str_c(
    unname(unique_pattern),
    collapse = ""
  )
  
  subject_seq <- sequence_AA1
  
  alignment <- pairwiseAlignment(pattern_seq, subject_seq, type=type)
  
  # Extract the aligned sequences
  aligned_pattern <- alignedPattern(alignment)
  residues_up_to_gap <- str_extract(as.character(aligned_pattern), "^[^-]*")
  
  shift <- str_locate(subject_seq, residues_up_to_gap)[1]
  
  if(is.na(shift)) {
    warning("Missing sequence residue compared to pdb structure")
    aligned_subject <- alignedSubject(alignment)
    residues_up_to_gap <- str_extract(as.character(aligned_subject), "^[^-]*")
    shift <- str_locate(subject_seq, residues_up_to_gap)[1]
  }
  
  required_adjustment <- pattern_resno[1] - shift
  
  #verbose
  if(verbose) {
    print(deparse(substitute(pdb_object)))
    print(aligned_pattern)
    print(alignedSubject(alignment)) 
    print(residues_up_to_gap)
    print(shift)
  }
  
  return(c(required_adjustment))
}

## map.seqresno.to.resno
library("Biostrings")

map.seqresno.to.resno <- function(pdb_object, verbose=FALSE, type="overlap") {
  subject_seq <- pdb_object$atom %>%
    filter(resid %in% toupper(AMINO_ACID_CODE)) %>%
    select(resno, resid) %>%
    distinct(.keep_all = TRUE) %>%
    mutate(resid = aa321(resid)) %>%
    pull(resid) %>%
    str_c(collapse = "")
  
  required_adjustment <- map.seqposno.to.resno(subject_seq, pdb_object, verbose=verbose, type=type)
  
  return(c(required_adjustment))
}

# create resmap
create.resmap <- function(pdb_atom_df, epitope_prediction_df,required_adjustment, join_by = c("resno")) {
  resmap_df <- pdb_atom_df %>% 
    select(resid,resno,chain) %>%
    unique() %>%
    filter(resid %in% toupper(AMINO_ACID_CODE),
           !is.na(chain),
           chain != FALSE) %>%
    right_join(epitope_prediction_df %>%
                mutate(resno = Position + required_adjustment),
              by = join_by
    )
  return(resmap_df)
}
