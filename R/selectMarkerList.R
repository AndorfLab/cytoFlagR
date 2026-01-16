# Load required packages 
suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(bslib))
suppressPackageStartupMessages(library(shinyjs))
suppressPackageStartupMessages(library(shinyWidgets))

# CSS string to style the app
css_function <- "
* {
    font-family: Verdana, sans-serif !important;
}
.btn {
    background-color:#539DDD !important; /* Change button colour */
    border-color: #539DDD !important;
    color: #000000 !important; 
}
.btn:hover {
    background-color: #2E6CA6 !important;
    border-color: #2E6CA6 !important;
}"


### Shiny module to let the user select a list of markers from a data frame column
# df: A data frame containing marker names in a single column or among multiple columns
# column_name: Optional string. If df has multiple columns, choose which column to use 
# Returns a character vector of the selected markers chosen by the user
selectMarkerList <- function(df, 
                             column_name = NULL) {
  
  # Validate that input is a data frame
  if (!is.data.frame(df)) {
    stop("Input must be a data frame")
  }
  
  # Require explicit column_name when df has more than one column
  if (ncol(df) != 1 && is.null(column_name)) {
    stop("Please specify column name to use")
  }
  
  # If no column_name provided and df has one column, use the first column
  if (is.null(column_name)) {
    values <- df[[1]]
  }
  
  else {
    
    # Ensure the specified column exists in df
    if (!column_name %in% colnames(df)) {
      stop("Specified column name does not exist in dataframe, please use the column name from your dataframe")
    }
    
    # Extract the values from the specified column
    values <- df[[column_name]]
  }
  
  # Get unique values and sort for display
  unique_values <- sort(unique(values))
  
  # Define the UI
  ui <- fluidPage(
    shinyjs::useShinyjs(),
    shiny::tags$head(tags$style(HTML(css_function))),
    theme = bslib::bs_theme(bg = "#222528", fg = "white", primary = "#2E6CA6"),
    shiny::titlePanel(h4("Choose markers for assessment",
                         style="font-size:25px;color:#539DDD;
                  font-family: monospace;font-weight: bold;"),
                      windowTitle = "cytoFlagR"),
    br(),
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        fluidRow(
          shiny::column(6, shiny::actionButton("selectAll", "Select All")),
          shiny::column(6, shiny::actionButton("deselectAll", "Deselect All"))
        ),
        hr(),
        shinyWidgets::prettyCheckboxGroup(inputId = "selected_values", 
                                          label = "Choose markers from list below:", 
                                          choices = unique_values, selected = unique_values)
      ),
      shiny::mainPanel(
        h4("Selected Markers"),
        verbatimTextOutput("selected"),
        shiny::actionButton("done", "Done")
      )
    )
  )
  
  # Define the server
  server <- function(input, 
                     output, 
                     session) {
    
    # Select all markers
    shiny::observeEvent(input$selectAll, {
      shinyWidgets::updatePrettyCheckboxGroup(session, "selected_values", 
                                              choices = unique_values, 
                                              selected = unique_values)
    })
    
    # Deselect all markers
    shiny::observeEvent(input$deselectAll, {
      shinyWidgets::updatePrettyCheckboxGroup(session, "selected_values", 
                                              choices = unique_values, 
                                              selected = character(0))
    })
    
    # Display selected markers in the main panel
    output$selected <- shiny::renderPrint({
      input$selected_values
    })
    
    # Finish and return the selected markers to the R session
    shiny::observeEvent(input$done, {
      shiny::stopApp(input$selected_values)
    })
  }
  
  # Launch the Shiny app
  selected_values <- runApp(shinyApp(ui, server))
  
  return(selected_values)
}
