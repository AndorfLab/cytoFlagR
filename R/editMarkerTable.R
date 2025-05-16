
suppressPackageStartupMessages(library(shiny))
suppressPackageStartupMessages(library(DT))
suppressPackageStartupMessages(library(bslib))
suppressPackageStartupMessages(library(readr))
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

editMarkerTable<-function(df) {
  ui<-fluidPage(
    shinyjs::useShinyjs(),
    shiny::tags$head(tags$style(HTML(css_function))),
    theme = bslib::bs_theme(bg = "#222528", fg = "white", primary = "#2E6CA6"),
    shiny::titlePanel(tags$h1("cytoFlagR automated marker threshold table",
                         style="font-size:25px;color:#539DDD;
                  font-family: monospace;font-weight: bold;"),
                      windowTitle = "cytoFlagR"),
    br(),
    tags$div(style="margin-top: -15px; margin-bottom: 20px;",
             tags$h6("To retain the current marker cutoff values, 
                     scroll to the bottom and click Download new cutoffs!", 
                     style="font-size:18px;color:#E1E1E1;
                  font-family: monospace;")),
    br(),
    DTOutput("marker_threshold_table"),
    verbatimTextOutput("updated_value"),
    shiny::downloadButton("save","Download new cutoffs")
  )
  server<-function(input, output, session) {
    
    values <- reactiveVal(df)
    
    output$marker_threshold_table<-DT::renderDT({
      DT::datatable(values(), editable = TRUE, options = list(pageLength = nrow(df)), 
                    class = 'cell-border stripe') %>% 
        DT::formatStyle(columns = 1:ncol(df), target = 'row',
                        backgroundColor = "#646668", 
                        color = '#FDFDFD')
    })
    shiny::observeEvent(input$marker_threshold_table_cell_edit, {
      info<-input$marker_threshold_table_cell_edit
      str(info)
      ### get the edit cutoff values
      updated_df<-values()
      updated_df[info$row, info$col]<-info$value
      values(updated_df) ### reactive update to cutoff values
      
      ### display values that were changed
      output$updated_value<-shiny::renderPrint({
        cat("Updated ", colnames(updated_df)[info$col],
            " in row ", info$row, "to", info$value)
      })
    })
    ### save updated marker threshold dataframe
    output$save<-shiny::downloadHandler(
      # filename = function() {
      #   
      # },
      file<-"updated_automated_marker_threshold_values.csv",
      content = function(file) {
        updated_auto_cutoffs<-values()
        ### save it to working global environment
        assign("updated_auto_cutoffs", updated_auto_cutoffs, envir = .GlobalEnv)
        
        ### save it to output folder
        readr::write_delim(updated_auto_cutoffs, file = file)
        shiny::showNotification("Marker thresholds updated, download and quit!", type = "message")
      }
    )
  }
  return(shinyApp(ui, server))
}
