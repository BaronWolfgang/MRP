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

# prediction viewer
prediction.viewer <- function(pdb_entry, directory, chain, sele_color_df=NULL) {
  #view_mode = c("colorValue","residueindex")
  
  file_path <- paste0(directory,pdb_entry,".pdb")
  sele <- paste0(":", chain)
    
  view <- NGLVieweR(file_path) %>%
    stageParameters(backgroundColor = "white", zoomSpeed = 1) %>%
    addRepresentation("cartoon",
                      param = list(colorScheme = "residueindex", 
                                   sele = sele)
    ) 
  
  if(!is.null(sele_color_df)) {
    # view_mode == "colorValue" expects a dataframe with columns $resno_string, containing the selected residues to be colored, and $color containing the desired color 
    
    for (row in 1:nrow(sele_color_df)) {
      if (sele_color_df$color[row] == "residueindex") {
        view <- view %>%
          addRepresentation(
            "surface",  # Add surface representation once
            param = list(
              sele = sele_color_df$resno_string[row],
              colorScheme = sele_color_df$color[row]
            )
          )
      }
      else {
        view <- view %>%
          addRepresentation(
            "surface",  # Add surface representation once
            param = list(
              sele = sele_color_df$resno_string[row],
              colorValue = sele_color_df$color[row]
            )
          )
      }
    } 
    
  }

  return(view)
}

#AlphaFold3 functions
gb2fasta <- function(gb_file_path, output_dir = NULL) {
  gb_lines <- readLines(gb_file_path)
  short_name <- sub(
    "\\.gb$", "", sub(
      "^\\d+\\.\\s*", "", basename(gb_file_path)
    )
  )
  
  origin_indices <- grep("^ORIGIN", gb_lines)
  end_indices <- grep("^//", gb_lines)
  header_prefix <- short_name
  
  fasta_content <- list()
  
  for(i in 1:length(origin_indices)) {
    start_index <- origin_indices[i] + 2
    end_index <- end_indices[i] - 1
    
    protein_sequence_raw <- gb_lines[start_index:end_index] 
    
    protein_sequence <- protein_sequence_raw %>%
      gsub("[0-9]+", "", .) %>%
      gsub("\\s+", "", .) %>%
      paste(collapse = "")
    
    header_line <- gb_lines[grep("^TITLE", gb_lines)[i]]
    header <- sub("^TITLE\\s+", "", header_line)
    
    fasta_content[[i]] <- paste0(">",header_prefix,"_", header, "\n", protein_sequence)
  }
  
  writeLines(unlist(fasta_content), paste0(output_dir,gb_file,"_msa_sequences.fasta"))
}
