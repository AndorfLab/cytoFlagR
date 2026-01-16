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

### Shiny module to let the user select a single reference marker from a data frame column
# df: A data.frame containing marker names
# column_name: String specifying which column of df contains the marker values
# Returns a selected reference marker
referenceMarkerSelection <- function(df, 
                                     column_name = NULL) {
  
  # Validate that input is a data frame
  if (!is.data.frame(df)) {
    stop("Input must be a data frame")
  }
  
  # Validate that the specified column exists in the data frame
  if (!column_name %in% colnames(df)) {
    stop("Specified column name does not exist in dataframe, 
           please use the column name from your dataframe")
  }
  
  # Extract the values from the specified column
  values <- df[[column_name]]
  
  # Get unique values and sort for display
  unique_values <- sort(unique(values))
  
  # Define the UI
  ui <- fluidPage(
    shinyjs::useShinyjs(),
    shiny::tags$head(tags$style(HTML(css_function))),
    theme = bslib::bs_theme(bg = "#222528", fg = "white", primary = "#2E6CA6"),
    shiny::titlePanel(h4("Reference marker selection for biaxial density plots",
                         style="font-size:20px;color:#539DDD;
                  font-family: monospace;font-weight: bold;"),
                      windowTitle = "cytoFlagR"),
    br(),
    br(),
    shinyWidgets::prettyRadioButtons(inputId = "selected_value", 
                                     label = "Choose a refence marker", 
                                     choices = unique_values, selected = unique_values[1], 
                                     status = "primary", animation = "smooth"),
    br(),
    h5("Selected Reference Marker: ",style="font-size:20px;color:#539DDD;"),
    verbatimTextOutput("selected"),
    br(),
    shiny::actionButton("done", "Save")
  )
  
  # Define the server
  server <- function(input, 
                     output, 
                     session) {
    
    # Display the currently selected reference marker
    output$selected <- shiny::renderPrint({
      input$selected_value
    })
    
    # Stop the app and return the selected marker when the user clicks Save
    shiny::observeEvent(input$done, {
      shiny::stopApp(input$selected_value)
    })
  }
  
  # Launch the Shiny app
  selected_value <- runApp(shinyApp(ui, server))
  
  return(selected_value)
}