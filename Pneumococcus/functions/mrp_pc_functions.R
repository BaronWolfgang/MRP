# Mapping functions
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
  #first_10_residues_pattern <- substr(as.character(aligned_pattern), 1, 10)
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
  if(verbose == TRUE) {
    print(deparse(substitute(pdb_object)))
    print(aligned_pattern)
    print(alignedSubject(alignment)) 
    print(residues_up_to_gap)
    print(shift)
  }
  
  return(c(required_adjustment))
}

## map.seqresno.to.resno
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
## bepipred
create.resmap.bepipred <- function(pdb_atom_df, bepipred_epitope_prediction_df,required_adjustment) {
  resmap_df <- pdb_atom_df %>% 
    select(resid,resno,chain) %>%
    unique() %>%
    filter(resid %in% toupper(AMINO_ACID_CODE),
           !is.na(chain),
           chain != FALSE) %>%
    right_join(bepipred_epitope_prediction_df %>%
                #select(Position,AminoAcid,EpitopeProbability) %>%
                mutate(resno = Position + required_adjustment),
              by = c("resno")
    )
  return(resmap_df)
}

# prediction viewer
prediction.viewer <- function(pdb_entry, directory, chain, predictions_df=NULL) {
  file_path <- paste0(directory,pdb_entry,".pdb")
  sele <- paste0(":", chain)
    
  view <- NGLVieweR(file_path) %>%
    stageParameters(backgroundColor = "white", zoomSpeed = 1) %>%
    addRepresentation("cartoon",
                      param = list(colorScheme = "residueindex", 
                                   sele = sele)
    ) 
  
  if(!is.null(predictions_df)) {
    for (row in 1:nrow(predictions_df)) {
      #print(predictions_df)
      #print(predictions_df$resno_string[row])
      #print(predictions_df$color[row])
      view <- view %>%
        addRepresentation(
          "surface",  # Add surface representation once
          param = list(
            sele = predictions_df$resno_string[row],
            colorValue = predictions_df$color[row]
          )
        )
    } 
    
  }
  
  return(view)
}

