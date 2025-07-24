### 
#functions that allow for easy chain selection and custom coloring 
###

# prediction viewer
prediction.viewer <- function(file_path, chain, annotation_df=NULL) {
  #view_mode = c("colorValue","residueindex")
  
  sele <- paste0(":", chain)
  
  view <- NGLVieweR(file_path) %>%
    stageParameters(backgroundColor = "white", zoomSpeed = 1) %>%
    addRepresentation("cartoon",
                      param = list(colorScheme = "residueindex", 
                                   sele = sele)
    ) 
  
  if(!is.null(annotation_df)) {
    # view_mode == "colorValue" expects a dataframe with columns $resno_string, containing the selected residues to be colored, and $color containing the desired color 
    
    for (row in 1:nrow(annotation_df)) {
      if (annotation_df$color[row] == "residueindex") {
        view <- view %>%
          addRepresentation(
            "surface",  # Add surface representation once
            param = list(
              sele = annotation_df$resno_string[row],
              colorScheme = annotation_df$color[row]
            )
          )
      }
      else {
        view <- view %>%
          addRepresentation(
            "surface",  # Add surface representation once
            param = list(
              sele = annotation_df$resno_string[row],
              colorValue = annotation_df$color[row]
            )
          )
      }
    } 
    
  }
  
  return(view)
}