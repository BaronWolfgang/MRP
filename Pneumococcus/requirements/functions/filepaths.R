pdbs_path = "data/pdb/"
if (!dir.exists(pdbs_path)) dir.create(pdbs_path, recursive = TRUE)

AF3_pdbs_path = "output/AlphaFold3/generated_pdb/vws/"
if (!dir.exists(AF3_pdbs_path)) dir.create(AF3_pdbs_path, recursive = TRUE)

DNA_alnfiles_path <- "data/bioedit/fas/DNA"
if (!dir.exists(DNA_alnfiles_path)) dir.create(DNA_alnfiles_path, recursive = TRUE)

AA_alnfiles_path <- "data/bioedit/fas/AA"
if (!dir.exists(AA_alnfiles_path)) dir.create(AA_alnfiles_path, recursive = TRUE)

AA_filtered_alnfiles_path <- "data/bioedit/fas/AA/filtered"
if (!dir.exists(AA_alnfiles_path)) dir.create(AA_filtered_alnfiles_path, recursive = TRUE)

pdb_files_list <- list.files(pdbs_path, pattern = "\\.pdb$") %>% sub("\\.pdb$", "", .)

AF3_pdb_files_list <- list.files(AF3_pdbs_path, pattern = "\\.pdb$") %>% sub("\\.pdb$", "", .)