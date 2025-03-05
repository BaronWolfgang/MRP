###
#convert .gb to .fasta to allow for Multiple Sequence Alignment
###

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
  
  writeLines(unlist(fasta_content), paste0(output_dir,short_name,"_msa_sequences.fasta"))
}

