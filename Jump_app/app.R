#-----------------------------------------------------------------------------------------
#
# Load and analyse squat or countermovement jump data
# Robert Schuster (ACU SPRINT)
# September 2022
#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
#-----------------------------------------------------------------------------------------

# TO DO ----------------------------------------------------------------------------------
# - Add 'Both' option to jumpType and automatic jump type detection
# ----------------------------------------------------------------------------------------

library(shiny)
library(DT)
source("Jump_functions.R")

# UI -------------------------------------------------------------------------------------
ui <- fluidPage(
  
  # Application title
  titlePanel(img(src = "ACU_logo.png", height = 70, width = 200)),
  
  # Sidebar
  sidebarLayout(
    sidebarPanel(
      h3("SPRINT Jump Analyser"),
      
      selectInput('jumpType', 'Select the type of jump to analyse', c('Countermovement jump','Squat jump'), selectize = F),
      fileInput("file", 
                "Select the files you want to analyse",
                multiple = T,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv")),
      tags$hr(), # horizontal line
      # length of period to determine CM threshold
      sliderInput("thl", "Quiet standing period length [s]:",
                  min = 0.1, max = 1,
                  value = 0.5, step = 0.1),
      tags$hr(), # horizontal line
      # Save filename
      textInput("fileName", "Enter the file name you want to save the results to:"),
      
      # Download button
      p('(You can select the rows you want to download by clicking them)'),
      downloadButton("downloadData", "Download selected rows", icon = icon("download"), 
                     style="color: #333; background-color: #FFF; border-color: #333")
    ),
    
    # Performance metrics table and graphs of each rep
    mainPanel(
      dataTableOutput('results'),
      uiOutput('repTabs')
    )
  )
)

# Server logic ---------------------------------------------------------------------------
server <- function(input, output) {
  # remove default input file size restriction (increase to 30MB)
  options(shiny.maxRequestSize = 30*1024^2)
  # Load files into workspace
  getData <- reactive({
    if (!is.null(input$file)) {
      numfiles = nrow(input$file)
      perfMet = list()
      ds <- list()
      for (i in 1:numfiles) {
        data <- importTrial(input$file[[i,'datapath']],input$file[[i,'name']],input$jumpType) # import and prepare trial
        data <- nReps(data) # check for multiple reps and cut trial accordingly
        data <- flight(data, input$thl) # determine start and end of flight and start of movement
        data <- qualityCheck(data)
        data <- perfMetrics(data) # calculate performance metrics
        
        perfMet[[i]] = data$pm
        ds[[data$fn]] <- data
      }
      ds$pm <- do.call(rbind, perfMet)
      return(ds)
    }
  })
  
  # Table of performance metrics
  output$results <- renderDataTable(
    getData()$pm,
    rownames = T,
    options = list(select = T,
                   searching = F,
                   scrollX = T,
                   scrollY = T)
  )
  
  # Plot each rep in a separate tab and print associated quality warnings
  repTabs <- function(fn) {
    reps <- names(getData()[[fn]])[grep(fn,names(getData()[[fn]]))]
    rTabs <- lapply(1:length(reps), function(r) {
      rn <- reps[r]
      tabPanel(paste('Rep',r), renderPlot({
        t <- getData()[[fn]][[rn]]$Time
        f <- getData()[[fn]][[rn]]$Total
        bm <- mean(f[1:(getData()[[fn]]$freq*1)]) # alternatively: getData()[[fn]]$bodymass
        sf <- getData()[[fn]]$flight[r,1]
        ef <- getData()[[fn]]$flight[r,2]
        sj <- getData()[[fn]]$sj[r]
        # plot
        plot(x = t, y = f, type = "l", lwd = 2, 
             xlab = "Time [s]", ylab = "Force [N]")
        abline(h = bm, col = "green", lwd = 2, lty = "dashed")
        abline(v = t[sf], col = "red", lwd = 2, lty = "dotted")
        abline(v = t[ef], col = "red", lwd = 2, lty = "dotted")
        abline(v = t[sj], col = "blue", lwd = 2, lty = "dotdash")
      }),
      # warning messages
      if (any(grepl(rn,names(getData()[[fn]]$warn)))) {
        msg <- paste(unlist(getData()[[fn]]$warn[[rn]]), collapse = '<br/>')
        HTML(paste("<b>Warning:</b>", msg, sep = '<br/>'))
      })
    })
    do.call(tabsetPanel, rTabs)
  }
  
  output$repTabs <- renderUI({
    if (!is.null(input$file)) {
      fns <- input$file$name
      fTabs <- lapply(1:length(fns), function(f) {
        if (is.null(fns)) {
          return(NULL)
        } else {
          fn <- basename(fns[f])
          tabPanel(fn, repTabs(fn))
        }
      })
      do.call(tabsetPanel, fTabs)
    }
  })
  
  # Save csv of performance metrics
  # https://stackoverflow.com/questions/70039664/reactive-element-as-filename-in-downloadhandler
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0(input$fileName, ".csv")
    },
    content = function(file) {
      write.csv(getData()$pm[input$results_rows_selected,], file)
    }
  )
}

# Run the application --------------------------------------------------------------------
shinyApp(ui = ui, server = server)
