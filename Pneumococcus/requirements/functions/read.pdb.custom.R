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
  
  # Try reading local PDB file
  if (file.exists(pdb_file)) {
    pdb_object <- tryCatch({
      read.pdb(file = pdb_file, rm.alt = FALSE)
    }, error = function(e) {
      message("Error reading local PDB file: ", pdb_file)
      return(NULL)  
    })
  }
  
  # Try reading from web if local fails
  if (is.null(pdb_object)) {
    pdb_object <- tryCatch({
      read.pdb(pdb_entry, rm.alt = FALSE)
    }, error = function(e) {
      message("Error reading PDB entry from web: ", pdb_entry)
      return(NULL)
    })
  }
  
  # Retry web load if pdbseq is missing
  if (!is.null(pdb_object) && is.null(pdbseq(pdb_object))) {
    message("Missing pdbseq for ", pdb_entry, " - Retrying from web")
    pdb_object <- tryCatch({
      read.pdb(pdb_entry, rm.alt = FALSE)
    }, error = function(e) {
      message("Retry failed for PDB entry: ", pdb_entry)
      return(NULL)
    })
  }
  
  # Try reading CIF if all PDB attempts fail
  if (is.null(pdb_object)) {
    message("Attempting to fetch and parse mmCIF for ", pdb_entry)
    pdb_object <- tryCatch({
      read.cif(pdb_entry)  # You need this function implemented already
    }, error = function(e) {
      message("Failed to read CIF for ", pdb_entry)
      return(NULL)
    })
  }
  
  return(pdb_object)
}
