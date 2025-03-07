#####
#code adapted from https://github.com/Grantlab/bio3d/blob/master/ver_devel/bio3d/R/read.cif.R 
#adapted to handle Alphafold3 .cif files and convert them to bio3d pdb objects
####

#parce_cif
parse_cif <- function(filename, maxlines = -1, multi = FALSE) {
  if (!file.exists(filename)) {
    return(list(error = "Error when reading file"))
  }
  
  # Initialize variables
  natoms <- 0
  models <- 0
  xyz <- numeric()
  
  # Define lists for storing atom data
  atom_data <- list(
    type = character(), eleno = integer(), elety = character(), alt = character(),
    resid = character(), chain = character(), resno = integer(), insert = character(),
    x = numeric(), y = numeric(), z = numeric(), o = numeric(), b = numeric(),
    segid = integer(), elesy = character(), charge = character(), model = integer()
  )
  
  # Read file
  lines <- readLines(filename, warn = FALSE)
  if (maxlines > 0) {
    lines <- head(lines, maxlines)
  }
  
  section_type <- ""
  columns_map <- list()
  prev_model <- -1
  
  for (line in lines) {
    line <- trimws(line)
    
    if (startsWith(line, "#")) {
      section_type <- ""
      columns_map <- list()
      next
    }
    
    if (startsWith(line, "_atom_site.")) {
      section_type <- "_atom_site"
      columns_map[[line]] <- length(columns_map) + 1
      next
    }
    
    if (section_type == "_atom_site" && (startsWith(line, "ATOM") || startsWith(line, "HETATM"))) {
      fields <- unlist(strsplit(line, "\\s+"))
      
      if (!"_atom_site.pdbx_PDB_model_num" %in% names(columns_map)) next
      curr_model <- as.integer(fields[columns_map[["_atom_site.pdbx_PDB_model_num"]]])
      
      if (curr_model != prev_model) {
        models <- models + 1
        prev_model <- curr_model
      }
      
      if (!multi && models > 1) {
        models <- 1
        break
      }
      
      required_columns <- c("_atom_site.Cartn_x", "_atom_site.Cartn_y", "_atom_site.Cartn_z")
      if (!all(required_columns %in% names(columns_map))) next
      
      tmpx <- as.numeric(fields[columns_map[["_atom_site.Cartn_x"]]])
      tmpy <- as.numeric(fields[columns_map[["_atom_site.Cartn_y"]]])
      tmpz <- as.numeric(fields[columns_map[["_atom_site.Cartn_z"]]])
      
      xyz <- c(xyz, tmpx, tmpy, tmpz)
      
      if (models == 1) {
        natoms <- natoms + 1
        atom_data$type <- c(atom_data$type, fields[columns_map[["_atom_site.group_PDB"]]])
        atom_data$eleno <- c(atom_data$eleno, as.integer(fields[columns_map[["_atom_site.id"]]]))
        atom_data$elety <- c(atom_data$elety, fields[columns_map[["_atom_site.label_atom_id"]]])
        atom_data$alt <- c(atom_data$alt, fields[columns_map[["_atom_site.label_alt_id"]]])
        atom_data$resid <- c(atom_data$resid, fields[columns_map[["_atom_site.label_comp_id"]]])
        atom_data$chain <- c(atom_data$chain, fields[columns_map[["_atom_site.auth_asym_id"]]])
        atom_data$resno <- c(atom_data$resno, as.integer(fields[columns_map[["_atom_site.auth_seq_id"]]]))
        atom_data$insert <- c(atom_data$insert, fields[columns_map[["_atom_site.pdbx_PDB_ins_code"]]])
        atom_data$x <- c(atom_data$x, tmpx)
        atom_data$y <- c(atom_data$y, tmpy)
        atom_data$z <- c(atom_data$z, tmpz)
        atom_data$o <- c(atom_data$o, as.numeric(fields[columns_map[["_atom_site.occupancy"]]]))
        atom_data$b <- c(atom_data$b, as.numeric(fields[columns_map[["_atom_site.B_iso_or_equiv"]]]))
        atom_data$segid <- c(atom_data$segid, as.integer(fields[columns_map[["_atom_site.label_entity_id"]]]))
        atom_data$elesy <- c(atom_data$elesy, fields[columns_map[["_atom_site.type_symbol"]]])
        atom_data$charge <- c(atom_data$charge, fields[columns_map[["_atom_site.pdbx_formal_charge"]]])
        atom_data$model <- c(atom_data$model, as.integer(fields[columns_map[["_atom_site.pdbx_PDB_model_num"]]]))
      }
    }
    if (length(atom_data$charge) == 0) {
      atom_data$charge <- rep(NA, length(atom_data$type))
    }
    atom <- data.frame(atom_data)
  }
  
  return(list(atom = atom, xyz = xyz, models = models))
}

#read.cif.custom
read.cif.custom <- function(file, maxlines = -1, multi=FALSE,
                            rm.insert=FALSE, rm.alt=TRUE, verbose=TRUE) {
  cl <- match.call()
  
  if(missing(file)) {
    stop("read.cif: please specify a CIF 'file' for reading")
  }
  
  if(!is.logical(multi)) {
    stop("read.cif: 'multi' must be logical TRUE/FALSE")
  }
  
  ## parse the CIF file with read_cif_data
  pdb <- parse_cif(file = file, maxlines = maxlines, multi = multi)
  
  pdb$atom <- pdb$atom %>%
    mutate(across(where(is.character), ~na_if(.x, "?")),
           across(where(is.character), ~na_if(.x, ".")),
           across(where(is.character), ~na_if(.x, "")))
  
  if(!is.null(pdb$error))
    stop(paste("Error in reading CIF file", file))
  else
    class(pdb) <- c("pdb") ## this should be cif perhaps?
  
  if(pdb$models > 1)
    pdb$xyz <- matrix(pdb$xyz, nrow=pdb$models, byrow=TRUE)
  
  pdb$xyz <- as.xyz(pdb$xyz)
  
  ## Remove 'Alt records'
  if (rm.alt) {
    if ( sum( !is.na(pdb$atom$alt) ) > 0 ) {
      first.alt <- sort( unique(na.omit(pdb$atom$alt)) )[1]
      cat(paste("   PDB has ALT records, taking",first.alt,"only, rm.alt=TRUE\n"))
      alt.inds <- which( (pdb$atom$alt != first.alt) ) # take first alt only
      if(length(alt.inds)>0) {
        pdb$atom <- pdb$atom[-alt.inds,]
        pdb$xyz <- trim.xyz(pdb$xyz, col.inds=-atom2xyz(alt.inds))
      }
    }
  }
  
  ## Remove 'Insert records'
  if (rm.insert) {
    if ( sum( !is.na(pdb$atom$insert) ) > 0 ) {
      cat("   PDB has INSERT records, removing, rm.insert=TRUE\n")
      insert.inds <- which(!is.na(pdb$atom$insert)) # rm insert positions
      pdb$atom <- pdb$atom[-insert.inds,]
      pdb$xyz <- trim.xyz(pdb$xyz, col.inds=-atom2xyz(insert.inds))
    }
  }
  
  if(any(duplicated(pdb$atom$eleno)))
    warning("duplicated element numbers ('eleno') detected")
  
  ## construct c-alpha attribute
  ca.inds <-  atom.select.pdb(pdb, string="calpha", verbose=FALSE)
  pdb$calpha <- seq(1, nrow(pdb$atom)) %in% ca.inds$atom
  
  ## set call
  pdb$call <- cl
  
  return(pdb)
}