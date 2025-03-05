###
# Sometimes using read.pdb() on a file downloaded using get.pdb() has a missing pdbseq entry. 
# This function tries to check the pdb file and recover the pdbseq
###

read.pdb.custom <- function(pdb_entry, pdb_path = NULL) {
  if (!is.null(pdb_path)) {
    pdb_file <- paste0(pdb_path, pdb_entry, ".pdb")
  } else {
    pdb_file <- paste0(pdb_entry, ".pdb") 
  }
  
  pdb_object <- NULL  # Default return value
  
  if (file.exists(pdb_file)) {
    pdb_object <- tryCatch({
      read.pdb(file = pdb_file, rm.alt = FALSE)
    }, error = function(e) {
      message("Error reading file: ", pdb_file, " - Trying to load PDB entry, ", pdb_entry, ", from web.")
      return(NULL)  
    })
  } 
  
  if (is.null(pdb_object)) {  # If local file fails, attempt to fetch from the web
    pdb_object <- tryCatch({
      read.pdb(pdb_entry, rm.alt = FALSE)
    }, error = function(e) {
      message("Error reading PDB entry: ", pdb_entry, "- Failed to load from web, returning NULL")
      return(NULL)
    })
  }
  
  if (!is.null(pdb_object) && is.null(pdbseq(pdb_object))) {
    print(paste("Missing pdbseq:", pdb_entry, "- Retrying to load from web using PDB entry:", pdb_entry))
    pdb_object <- tryCatch({
      read.pdb(pdb_entry, rm.alt = FALSE)
    }, error = function(e) {
      message("Retry failed for PDB entry: ", pdb_entry, "- Returning NULL")
      return(NULL)
    })
  }
  
  return(pdb_object)
}