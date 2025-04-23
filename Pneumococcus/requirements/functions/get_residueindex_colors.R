rainbow_ngl <- c("red", "orange", "yellow", "green", "blue")

get_residueindex_colors <- function(residues, rev = FALSE) {
  n <- length(residues)
  palette_colors <- colorRampPalette(rainbow_ngl)(n)
  if (rev) palette_colors <- rev(palette_colors)
  data.frame(resno = residues, color = palette_colors)
}

#v2
get_residueindex_colors_v2 <- function(residues, rev = FALSE) {
  #residues <- sort(unique(residues))
  min_res <- min(residues, na.rm = TRUE)
  max_res <- max(residues, na.rm = TRUE)
  
  # Define evenly spaced breaks across the residue range
  breaks <- seq(min_res, max_res, length.out = 5)
  
  # Define rainbow palette
  palette <- rainbow_ngl
  if (rev) palette <- rev(palette)
  
  # Create color interpolation function
  color_fun <- colorRamp2(breaks, palette, space= "sRGB")
  
  raw_colors <- color_fun(residues)
  clean_colors <- substr(raw_colors, 1, 7)  # Keep only #RRGGBB
  
  data.frame(resno = residues, color = clean_colors)
}