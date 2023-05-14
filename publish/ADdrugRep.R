# Load R packages
library(shiny)
library(shinythemes)


# Define UI
ui <- fluidPage(theme = shinytheme("cerulean"),
                # Page header
                #headerPanel("Drug Repurposing for Alzheimer's Disease based on DeepWalk"),
                navbarPage(
                  # theme = "cerulean",  # <--- To use a theme, uncomment this
                  "Drug Repurposing for Alzheimer's Disease based on DeepWalk",
                  tabPanel("Search for drug score",
                           sidebarPanel(
                             tags$h3("Input:"),
                             textInput("txt1", "Drug Name or ID:", ""),
                             textInput("txt2", "Surname:", ""),
                             
                           ), # sidebarPanel
                           mainPanel(
                             h1("Header 1"),
                             
                             h4("Output 1"),
                             verbatimTextOutput("txtout"),
                             
                           ) # mainPanel
                           
                  ), # Navbar 1, tabPanel
                  tabPanel("Drug target interaction", "This panel is intentionally left blank"),
                  tabPanel("AD associated genes", "This panel is intentionally left blank")
                  
                ) # navbarPage
) # fluidPage


# Define server function  
server <- function(input, output) {
  
  output$txtout <- renderText({
    paste( input$txt1, input$txt2, sep = " " )
  })
} # server


# Create Shiny object
shinyApp(ui = ui, server = server)