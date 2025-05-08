pdbs_path = "data/pdb/"
if (!dir.exists(pdbs_path)) dir.create(pdbs_path, recursive = TRUE)

DNA_alnfiles_path <- "data/bioedit/fas/DNA"
if (!dir.exists(DNA_alnfiles_path)) dir.create(DNA_alnfiles_path, recursive = TRUE)

AA_alnfiles_path <- "data/bioedit/fas/AA"
if (!dir.exists(AA_alnfiles_path)) dir.create(AA_alnfiles_path, recursive = TRUE)