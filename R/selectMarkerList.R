
suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(bslib))
suppressPackageStartupMessages(library(shinyjs))
suppressPackageStartupMessages(library(shinyWidgets))

css_function<-"
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

selectMarkerList<-function(df, column_name = NULL) {
  ## input validation
  if(!is.data.frame(df)) {
    stop("Input must be a data frame")
  }
  if(ncol(df) != 1 && is.null(column_name)) {
    stop("Please specify column name to use")
  }
  if(is.null(column_name)) {
    values<-df[[1]]
  }
  else {
    if(!column_name %in% colnames(df)) {
      stop("Specified column name does not exist in dataframe, please use the column name from your dataframe")
    }
    values<-df[[column_name]]
  }
  
  unique_values<-sort(unique(values))
  
  ## ui function
  ui<-fluidPage(
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
  ## server function
  server<-function(input, output, session) {
    
    shiny::observeEvent(input$selectAll, {
      shinyWidgets::updatePrettyCheckboxGroup(session, "selected_values", 
                                              choices = unique_values, 
                                              selected = unique_values)
    })
    
    shiny::observeEvent(input$deselectAll, {
      shinyWidgets::updatePrettyCheckboxGroup(session, "selected_values", 
                                              choices = unique_values, 
                                              selected = character(0))
    })
    
    output$selected<-shiny::renderPrint({
      input$selected_values
    })
    
    shiny::observeEvent(input$done, {
      shiny::stopApp(input$selected_values)
    })
  }
  selected_values<-runApp(shinyApp(ui, server))
  return(selected_values)
}
