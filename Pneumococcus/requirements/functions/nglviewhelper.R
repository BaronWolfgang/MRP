### 
#functions that allow for easy chain selection and custom coloring 
###

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