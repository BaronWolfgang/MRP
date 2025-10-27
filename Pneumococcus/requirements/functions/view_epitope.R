#!/usr/bin/env Rscript

# ============================================
# Generate a self-contained interactive HTML UI
# ============================================
# Usage:
#   Rscript generate_structure_widget.R protein_specific.csv structure.pdb output_widget.html
# Example:
#   Rscript generate_structure_widget.R data/my_protein.csv data/pdb/my_protein.pdb output/my_protein.html

# ---- 1. Parse command-line arguments ----
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: Rscript generate_structure_widget.R <protein_specific_csv> <pdb_file> <output_html>")
}

protein_specific_csv <- args[1]
pdb_file <- args[2]
output_html <- args[3]

if (!file.exists(protein_specific_csv)) stop("CSV not found: ", protein_specific_csv)
if (!file.exists(pdb_file)) stop("PDB not found: ", pdb_file)

# ---- 2. Load required packages ----
required_packages <- c(
  "bio3d", "NGLVieweR", "htmlwidgets", "htmltools", "shiny", "dplyr", "ggplot2", "readr", "scales", "stringr", "tools", "grDevices"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  library(pkg, character.only = TRUE)
}

# ---- 3. Source color utility ----
color_function_path <- "requirements/functions/get_residueindex_colors.R"
if (!file.exists(color_function_path)) {
  stop("Required function file not found: ", color_function_path)
}

message("Sourcing color function from: ", color_function_path)
source(color_function_path)

# ---- 4. Load data ----
message("Reading CSV: ", protein_specific_csv)
data <- readr::read_csv(protein_specific_csv, show_col_types = FALSE)

if (!all(c("resno", "chain") %in% colnames(data))) {
  stop("CSV must contain columns: resno, chain")
}

# ---- 5. Define Shiny UI ----
ui <- shiny::fluidPage(
  shiny::titlePanel("Protein Structure Viewer"),
  shiny::sidebarLayout(
    shiny::sidebarPanel(
      h4("Display options"),
      selectInput(
        "struc_representation", "Structure representation",
        c("cartoon", "ribbon", "surface", "ball+stick"), selected = "cartoon"
      ),
      selectInput(
        "struc_color_scheme", "Base color scheme",
        c("uniform", "residueindex", "bfactor"), selected = "uniform"
      ),
      selectInput(
        "viewmode", "Residue score column",
        choices = names(data)[sapply(data, is.numeric)]
      ),
      checkboxInput("normalize", "Normalize values", value = FALSE),
      selectInput("viewmethod", "Coloring method", c("score", "topn", "threshold")),
      conditionalPanel(
        condition = "input.viewmethod == 'topn'",
        sliderInput("topn_value", "Top n%", min = 0, max = 100, value = 25, step = 5)
      ),
      conditionalPanel(
        condition = "input.viewmethod == 'threshold'",
        uiOutput("threshold_slider_ui")
      ),
      actionButton("apply", "Apply Coloring")
    ),
    shiny::mainPanel(
      plotOutput("score_plot", height = "250px"),
      NGLVieweROutput("structure", height = "600px")
    )
  )
)

# ---- 6. Server logic ----
server <- function(input, output, session) {
  
  data_reactive <- reactive({
    d <- data
    selected_column <- input$viewmode
    d$view <- d[[selected_column]]
    if (isTRUE(input$normalize)) {
      d$view <- scales::rescale(d$view, to = c(0, 1))
    }
    d
  })
  
  output$threshold_slider_ui <- renderUI({
    req(input$viewmode)
    d <- data_reactive()
    values <- d$view
    sliderInput("threshold_value", "Threshold",
                min = min(values, na.rm = TRUE),
                max = max(values, na.rm = TRUE),
                value = mean(values, na.rm = TRUE), step = 0.01)
  })
  
  output$score_plot <- renderPlot({
    req(input$viewmode)
    d <- data_reactive()
    color_df <- get_residueindex_colors_v2(d$resno)
    d <- dplyr::left_join(d, color_df, by = "resno")
    ggplot2::ggplot(d, aes(x = resno, y = view)) +
      geom_col(aes(fill = color), width = 1) +
      scale_fill_identity() +
      labs(x = "Residue Number", y = input$viewmode,
           title = paste("Per-residue", input$viewmode)) +
      theme_bw()
  })
  
  output$structure <- renderNGLVieweR({
    NGLVieweR(pdb_file) %>%
      addRepresentation(
        type = input$struc_representation,
        param = list(
          name = "base",
          colorScheme = input$struc_color_scheme
        )
      ) %>%
      stageParameters(backgroundColor = "white") %>%
      setQuality("medium") %>%
      setFocus(0)
  })
  
  observeEvent(input$apply, {
    d <- data_reactive()
    proxy <- NGLVieweR_proxy("structure")
    
    for (i in 1:21) {
      proxy %>% removeSelection(paste0("sel", i))
    }
    proxy %>% removeSelection("threshold") %>% removeSelection("topn")
    
    if (input$viewmethod == "score") {
      color_df <- get_residueindex_colors_v2(d$resno)
      d <- dplyr::left_join(d, color_df, by = "resno")
      color_groups <- d %>%
        group_by(color) %>%
        summarise(sele = paste0(resno, collapse = " "), .groups = "drop")
      
      for (i in seq_len(nrow(color_groups))) {
        proxy <- proxy %>%
          addSelection("surface",
                       param = list(
                         name = paste0("sel", i),
                         sele = paste0("/0 AND :", unique(na.omit(d$chain)),
                                       " AND (", color_groups$sele[i], ")"),
                         colorValue = color_groups$color[i]
                       ))
      }
    }
    
    if (input$viewmethod == "threshold") {
      req(input$threshold_value)
      top_residues <- d %>% filter(view >= input$threshold_value)
      sele_string <- paste0(top_residues$resno, collapse = " ")
      proxy <- proxy %>%
        addSelection("surface",
                     param = list(
                       name = "threshold",
                       sele = paste0("/0 AND :", unique(na.omit(d$chain)),
                                     " AND (", sele_string, ")"),
                       colorScheme = "residueindex"
                     ))
    }
    
    if (input$viewmethod == "topn") {
      req(input$topn_value)
      cutoff <- quantile(d$view, probs = 1 - input$topn_value / 100, na.rm = TRUE)
      top_residues <- d %>% filter(view >= cutoff)
      sele_string <- paste0(top_residues$resno, collapse = " ")
      proxy <- proxy %>%
        addSelection("surface",
                     param = list(
                       name = "topn",
                       sele = paste0("/0 AND :", unique(na.omit(d$chain)),
                                     " AND (", sele_string, ")"),
                       colorScheme = "residueindex"
                     ))
    }
  })
}

# ---- 7. Save as self-contained HTML ----
message("Rendering self-contained Shiny HTML app → ", output_html)
app <- shinyApp(ui, server)
htmlwidgets::saveWidget(shiny::browsable(app), file = output_html, selfcontained = TRUE)
message("Done: ", output_html)
