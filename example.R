library(shiny)
setwd("/users/PCON0041/shunian/shunian/AD_drugRep_project_202209/Code/Rshiny/")
ui <- fluidPage(
  h1("Hello, world!"),
  img(src = "flowers-276014_960_720.jpg",
      align = "left",
      alt = "Test Image")
)
server <- function(input, output, session) {
}
shinyApp(ui, server)