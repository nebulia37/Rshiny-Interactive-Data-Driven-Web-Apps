data(mtcars)

choices = data.frame(
  var = names(mtcars),
  num = 1:length(names(mtcars))
)
# List of choices for selectInput
mylist <- as.list(choices$num)
# Name it
names(mylist) <- choices$var
# ui
ui <- fluidPage(
  # ... other UI elements ...
  fluidRow(
    column(cyjShinyOutput('cyjShiny'), width=10),
    column(htmlOutput("legend"),width=2),
    tags$style(HTML("
    .cytoscape-legend {
      padding: 10px;
      border: 1px solid #ccc;
    }
    .cytoscape-legend-item {
      display: flex;
      align-items: center;
      margin-bottom: 5px;
    }
    .cytoscape-legend-item .color-box {
      width: 20px;
      height: 20px;
      margin-right: 5px;
    }
  "))
  )
)

create_legend <- function() {
  # Add legend items based on your network's node and edge styles
  # Adjust the color and label values to match your style
  legend_items <- list(
    list(color = "#EA3423", label = "AD node"),
    list(color = "#F2A93B", label = "Gene node"),
    list(color = "#B5D7E4", label = "Drug node")
  )
  
  legend_html <- lapply(legend_items, function(item) {
    tags$div(
      class = "cytoscape-legend-item",
      tags$div(class = "color-box", style = paste0("background-color: ", item$color, ";")),
      tags$div(item$label)
    )
  })
  
  return(legend_html)
}

# server
server <- function(input, output) {
  output$legend <- renderUI({
    create_legend()
  })
  
  load("network_forEachDrug.RData") # list of edges and data.frame of nodes
  
  
  networkInput <- reactive({ 
    single_drug_name = "efavirenz"
    
    tbl.edges <- network_forEachDrug[[1]][single_drug_name][[1]]
    colnames(tbl.edges) <- c("source", "target","interaction")
    indiv_net_node <- unique(c(tbl.edges$source,tbl.edges$target))
    tbl.nodes <- network_forEachDrug[[2]] %>% filter(id %in% indiv_net_node)
    #tbl.nodes <- network_forEachDrug[[2]][single_drug_name]
    graph.json <- toJSON(dataFramesToJSON(tbl.edges, tbl.nodes), auto_unbox=TRUE)
    
    graph.json
  })
  
  output$cyjShiny <- renderCyjShiny({
    StyleFile='ADrep_style.js'
    cyj <- cyjShiny(graph = networkInput(), layoutName = "cola", styleFile = StyleFile)
    cyj
  })
}
# app
shinyApp(ui = ui, server = server)