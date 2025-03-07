### 
#functions used to generate JSON files that can be used as input for the AlphaFold3 server 
###

parse_fasta <- function(fasta_file) {
  if (!file.exists(fasta_file)) {
    stop("Error: FASTA file does not exist at the specified path: ", fasta_file)
  }
  
  fasta_lines <- readLines(fasta_file)  
  sequences <- list()
  
  current_name <- NULL
  current_sequence <- ""
  
  for (line in fasta_lines) {
    if (startsWith(line, ">")) {  
      if (!is.null(current_name)) {
        
        sequences <- append(sequences, list(list(name = current_name, sequence = current_sequence)))
      }
      current_name <- substr(line, 2, nchar(line))  
      current_sequence <- ""  
    } else {
      current_sequence <- paste0(current_sequence, line)  
    }
  }
  
  
  if (!is.null(current_name)) {
    sequences <- append(sequences, list(list(name = current_name, sequence = current_sequence)))
  }
  
  return(sequences)
}

# Function to generate the correct JSON structure for AlphaFold3s
generate_json <- function(fasta_file, suffix = NULL) {
  sequences <- parse_fasta(fasta_file)
  
  json_data <- lapply(sequences, function(seq_data) {
    
    clean_name <- sub("([^(\\s]+)(\\s|\\().*", "\\1", seq_data$name)
    clean_name <- gsub("\\s+", "", clean_name)
    
    if(!is.null(suffix)) {
      clean_name <- paste(clean_name, suffix)
    }
    
    list(
      name = as.character(clean_name),  
      modelSeeds = list(),  
      sequences = list(
        list(
          proteinChain = list(
            sequence = as.character(seq_data$sequence),  
            count = 1 
          )
        )
      )
    )
  })
  
  return(json_data)
}


# Function to write the JSON data to a file
write_json_file <- function(json_data, output_file) {
  output_dir <- dirname(output_file)
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  write_json(json_data, output_file, pretty = TRUE, auto_unbox = TRUE) #auto_unbox is very important as otherwise the Alphafoldserver will not accept the .json
}