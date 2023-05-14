# Load R packages
library(shiny)
library(shinythemes)
library(cyjShiny)
library(DT)
library(dplyr)
library(graph)
library(BiocGenerics)
library(parallel)
library(BiocManager)
options(repos = BiocManager::repositories())


all_druglist = read.delim(file = "new.ranked_score_sum.txt",stringsAsFactors = F,sep = "\t")
choices_drug = data.frame(
  drug = all_druglist$Drug,
  var = all_druglist$Drug
)
# List of choices for selectInput
druglist <- as.list(choices_drug$drug)
# Name it
names(druglist) <- choices_drug$var

drug_target_data = read.delim(file = "new.filtered_Final_drug_target_info.txt",stringsAsFactors = F, sep = "\t")
all_genelist = unique(drug_target_data$Gene_name)
choices_gene =  data.frame(
  gene = all_genelist,
  var = all_genelist
)

genelist <- as.list(choices_gene$gene)
names(genelist) <- choices_gene$var

# Define UI
ui <- fluidPage(
          theme = shinytheme("united"),
          titlePanel("Drug repurposing for Alzheimer's Disease based on DeepWalk"),
          navbarPage("",
                     tabPanel(icon("home"),
                              fluidRow(
                                column(
                                  br(),
                                  p("Through this application, it is intended to represent the result of our drug repurposing for Alzheimer's Disease(AD) project. 
                      We constructed a graph  that contains", tags$span(tags$b(" protein-protein interactions"), style = "color: red;"), ",",  tags$span(tags$b("drug-target interactions"), style = "color: red;"),"and",  tags$span(tags$b("AD-associated genes"), style = "color: red;"),". 
                      Then, we apply the", tags$span(tags$b("DeepWalk"), style = "color: red;"), "method on the graph to learn node representations. 
                      The model then scores and ranks drug candidates based on ", tags$span(tags$b("their significance of paths to AD node"), style = "color: red;"), "within the network through the learned embeddings. Users can search for the rank information of a certain drug in 'Search drug rank' section",
                                    style="text-align:justify;color:black;background-color:lavender;padding:15px;border-radius:10px"),
                                  br(),
                                  p("Besides the result of drug repurposing, we also represent the ", tags$span(tags$b("reliable drug target interactions"), style = "color: red;"), " assembled from four databases and ", tags$span(tags$b("AD associated gene list"), style = "color: red;"), " we obtained from multiple sources, the detail of which are in the ", tags$span(tags$b("'Drug target interactions'"), style = "color: red;"), "and",tags$span(tags$b(" 'AD associated genes' "), style = "color: red;")," section respectively.",style="text-align:justify;color:black;background-color:papayawhip;padding:15px;border-radius:10px"),
                                  width=7),
                                column(tags$img(src="flowchart.jpeg",width="600px",height="300px"),width=5)),
                              hr(),
                              tags$style(".fa-database {color:#E87722}"),
                              h3(p(em("Top50 ranked drugs repurposed for AD "),icon("database",lib = "font-awesome"),style="color:black;text-align:center")),
                              fluidRow(column(DT::dataTableOutput("DrugData"),
                                              width = 12)),
                              hr(),
                              p(em("Developed by"),br("Shunian Xiang"),style="text-align:center; font-family: times")), # tabPanel 1 for home
                     tabPanel("Search drug rank", 
                    sidebarPanel(
                      tags$h3("Input drug name for search:"),
                      selectInput("drug_name_id", label = "Drug Name",choices = druglist),
                     hr()),
                    mainPanel(
                      h2("Get the rank of a drug by searching"),
                      h4("The higher the rank, the more chance that the drug may have the potential to treat AD"),
                      #fluidRow(column(DT::dataTableOutput("single_drug_table"),
                      #                width = 12)),
                      # output rank and score of the drug, and the network of the drug to AD
                      uiOutput("ui_panel2"),
                      hr()
                    )), # tabPanel 2 for drug repurposing result search 
              tabPanel("Drug target interactions",
                    sidebarPanel(
                      tags$h3("Input drug name/gene symbol for search:"),
                      radioButtons("search_type", label ="",choiceNames = list("Search drug","Search gene"),choiceValues = list(
                        "drug", "gene"),selected = "drug"),
                      #textInput("search_name_id", "Drug Name or Gene Symbol:",value = "tacrine"),
                      selectInput("search_name_id", "Drug Name or Gene Symbol:",choices  = druglist)),
                       # sidebarPanel
                    mainPanel(
                      h2("Reliable drug target interactions for 1591 FDA approved drugs"),
                      verbatimTextOutput("txtout_3"),
                      #fluidRow(column(DT::dataTableOutput("single_drug_table"),
                      #                width = 12)),
                      #tableOutput("outTable"),
                      fluidRow(
                  column(
                    br(),
                      uiOutput("descript_drugTarget",style="text-align:justify;color:black;background-color:lavender;padding:15px;border-radius:10px"),
                    br(),
                    p("Users can obtain the drug target interactions for certain drug or gene by searching the drug name or gene symbol.",style="text-align:justify;color:black;background-color:papayawhip;padding:15px;border-radius:10px"),
                      width=12),
                      uiOutput("ui_panel3"),
                      hr()
                    ))), # tabPanel 3 for drug target search
              tabPanel("AD associated genes",
                  sidebarPanel(
                    tags$h3("Input evidence sources:"),
                    selectInput("ADgene_source",p("Please select the evidence source:",style="color:black; text-align:center"),choices=c("All"=1,"Experimentally validated"=2,"Large scale GWAS"=3,"Database"=4)),
                    h4("References"),
                    p("[1] Fang J, Zhang P. et al. Endophenotype-based in silico network medicine discovery combined with insurance record data mining identifies sildenafil as a candidate drug for Alzheimer's disease. Nat Aging. 2021 Dec;1(12):1175-1188."),
                    p("[2] Lambert, J.-C. et al. Meta-analysis of 74,046 individuals identifies 11 new susceptibility loci for Alzheimer’s disease. Nat. genetics 45, 1452 (2013)."),
                    p("[3] Marioni RE, Harris SE, Zhang Q, et al. GWAS on family history of Alzheimer’s disease. Transl Psychiatry. 2018;8:99."),
                    p("[4] Jansen IE, Savage JE, Watanabe K, et al. Genome-wide meta-analysis identifies new loci and functional pathways influencing Alzheimer’s disease risk. Nat Genet. 2019;51:404–413."),
                    p("[5] Kunkle BW, Grenier-Boley B, Sims R, et al. Genetic meta-analysis of diagnosed Alzheimer's disease identifies new risk loci and implicates Abeta, tau, immunity and lipid processing. Nat Genet. 2019;51:414–430."),
                    p("[6] De Rojas I, Moreno-Grau S, Tesi N, et al. Common variants in Alzheimer’s disease and risk stratification by polygenic risk scores. Nat Commun. 2021;12:3417."),
                    p("[7] Wightman, D.P., Jansen, I.E., Savage, J.E. et al. A genome-wide association study with 1,126,563 individuals identifies new risk loci for Alzheimer’s disease. Nat Genet, 1276–1282 (2021)."),
                    p("[8] Bellenguez C, Küçükali F, Jansen IE, et al. New insights into the genetic etiology of Alzheimer’s disease and related dementias. Nat Genet. 2022;54:412–436.")
                    #selectInput("PruebaAnalitica",p("Please select the test you want to try:",style="color:black; text-align:center"),choices=c("Shapiro-Wilk"=1,"Anderson-Darling"=2,"Cramér-von Mises"=3,"Kolmogorov-Smirnov"=4,"Jarque-Bera"=5))
                  ),
                  mainPanel(h2("Comprehensive AD associated gene list"),fluidRow(
                  column(
                    width = 12,
                    br(),
                      uiOutput("descript_ADgenes",style="text-align:justify;color:black;background-color:lavender;padding:15px;border-radius:10px"),
                    br())),
                  DT::dataTableOutput("AD_geneTable"))), # tabPanel 4 for AD genes
              tabPanel("Download",
                       h3("Full Data Download",style="text-align:justify;color:dark grey;background-color:lavender;padding:15px;border-radius:10px"),
                       h4("Download drugs ranked by our drug repurposing method for Alzheimer's Disease"),
                       downloadButton("downloadData1", "Click to Save"),
                       h4("Download drug target interactions"),
                       downloadButton("downloadData2", "Click to Save"),
                       h4("Download AD associated genes"),
                       downloadButton("downloadData3", "Click to Save"),
                       h4("Download protein protein interactions"),
                       downloadButton("downloadData4", "Click to Save")))
)
# Define server function  
create_legend <- function() {
  # Add legend items based on your network's node and edge styles
  # Adjust the color and label values to match your style
  legend_items <- list(
    list(color = "#EA3423", label = "AD node"),
    list(color = "#B5D7E4", label = "Gene node"),
    list(color = "#F2A93B", label = "Drug node")
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

server <- function(input, output, session) {
    #####-----------------------------
    #for tabPanel 1 (homepage) output table
    ####------------------------------
    result_data <-read.delim("./new.ranked_score_sum.txt",sep = "\t",header = TRUE)
    top50_drugs = result_data[1:50,]
    output$DrugData <- DT::renderDataTable(
        server = FALSE,
        DT::datatable({
            top50_drugs
        },
        extensions = 'Buttons',
        options = list(lengthMenu=list(c(5,15,20),c('5','15','20')),
                     pageLength=10,
                     initComplete = JS(
                     "function(settings, json) {",
                     "$(this.api().table().header()).css({'background-color': 'moccasin', 'color': '1c1b1b'});",
                     "}"),
                     dom = 'Bfrtip',
                     buttons = c('copy', 'csv', 'excel'),
                   columnDefs=list(list(className='dt-center',targets="_all"))
        ),
        #filter = "top",
        selection = 'multiple',
        style = 'bootstrap',
        class = 'cell-border stripe',
        rownames = FALSE,
        colnames = c("Rank","NodeID","Drug","Score")
    ))
  
   #####-----------------------------
    #for tabPanel2 (search drug) output table and plot
    ####------------------------------
    datasetInput_2 <- reactive({ 
        single_drug_data_2 = result_data %>% filter(Drug==tolower(input$drug_name_id))  
        single_drug_data_2
    })

    output$single_drug_table =  DT::renderDataTable(
      datasetInput_2(),
      options = list(paging = FALSE, searching = FALSE, dom = 't', initComplete = JS(
        "function(settings, json) {",
        "$(this.api().table().header()).css({'background-color': 'moccasin', 'color': '1c1b1b'});",
        "}")),
      #style = 'bootstrap',
      
      class = 'cell-border stripe',
      rownames = FALSE,
      colnames = c("Rank","NodeID","Drug","Score")
    )

    output$ui_panel2 <- renderUI({
        
        if( nrow(datasetInput_2()) == 0)
            print("The input is not in the drug list")
        else{
            if( nrow(datasetInput_2()) != 0){
              fluidRow(
                DT::dataTableOutput("single_drug_table"),
                hr(),
                h4("Below is the visulization of all the paths linking the drug node and AD node with a path length<=3:"),
                # output the network bewtween the drug and AD node
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
                # output the table of drug rank info
                
              )
        }
        }
    })
  
    
    
    ## load RData contains list of dataframes of edges and nodes info for each drug from the drug to AD 
    load("network_forEachDrug.RData") # list of edges and data.frame of nodes
    
    
    networkInput <- reactive({ 
      single_drug_name = tolower(input$drug_name_id)
     
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
    
    output$legend <- renderUI({
      create_legend()
    })
    #####-----------------------------
    #for tabPanel3 (search drug target) output table and plot
    ####------------------------------
    drug_target_data <- read.delim("./new.filtered_Final_drug_target_info.txt", sep = "\t", header=TRUE)
    
    url1_drug <- a("CHEMBL", href="https://www.ebi.ac.uk/chembl/")
    url2_drug <- a("BindingDB", href="https://www.bindingdb.org/rwd/bind/index.jsp")
    url3_drug<- a("the Therapeutic Target Database",href="https://db.idrblab.net/ttd/")
    url4_drug<- a("IUPHAR/BPS Guide to PHARMACOLOGY",href="https://www.guidetopharmacology.org/")
    output$descript_drugTarget <- renderUI({
      tagList("We assembled drug-target interactions and bioactivity data from four commonly used databases:",url1_drug, "(v31),",url2_drug, " (downloaded in November 2022), ", url3_drug," (downloaded in November 2022) and", url4_drug, " (downloaded in November 2022). Of the interactions present in these databases, we retained only those that satisfied the following criteria: (1) binding affinities, including Ki, Kd, IC50, or EC50 ≤10 μM (10,000 nM). (2) gene targets and their respective proteins must have a unique UniProt accession number. (3) protein targets be marked as ‘reviewed’ in the UniProt database22. (4) protein targets found in Homo Sapiens. As a result of these criteria, ",tags$span(tags$b("10,701"), style = "color: red;")," drug-gene interactions between ",tags$span(tags$b("1,591"), style = "color: red;")," US FDA-approved drugs and ",tags$span(tags$b("1,254"), style = "color: red;")," unique genes were obtained and integrated with our PPI network")
    })
    
    observeEvent(input$search_type, {
      if (input$search_type == "drug") {
        updateSelectInput(session, "search_name_id", choices = druglist)
      } else if (input$search_type == "gene") {
        updateSelectInput(session, "search_name_id", choices = genelist)
      }
    })

    datasetInput_3 <- reactive({
      if(input$search_type=="drug"){
        single_drug_data_3 = drug_target_data %>% filter(Drug_name==input$search_name_id)
      }else if(input$search_type=="gene"){
        single_drug_data_3 = drug_target_data %>% filter(Gene_name==input$search_name_id)
      }
      single_drug_data_3
    })
    
    output$ui_panel3 <- renderUI({
      if(nrow(datasetInput_3()) == 0){
        print("The input is not in our list")
      }else{
        DT::dataTableOutput("drug_targetTable")
      }
    })
    
    output$drug_targetTable = DT::renderDT(
      server = FALSE,
      datasetInput_3(),
        extensions = 'Buttons',
        options = list(lengthMenu=list(c(5,15,20),c('5','15','20')),
                       pageLength=10,
                       initComplete = JS(
                         "function(settings, json) {",
                         "$(this.api().table().header()).css({'background-color': 'moccasin', 'color': '1c1b1b'});",
                         "}"),
                       dom = 'Bfrtip',
                       buttons = c('copy', 'csv', 'excel'),
                       columnDefs=list(list(className='dt-center',targets="_all"))
        ),
        selection = 'multiple',
        style = 'bootstrap',
        class = 'cell-border stripe',
        rownames = FALSE,
        colnames = c("Drug_name", "DrugBank_ID", "CHEMBL_ID", "Gene_name","UniProt","Source","EntrezID")
        )
  #####-----------------------------
  #for tabPanel4 (present AD genes) output table
  ####------------------------------
  url1 <- a("the OMIM database", href="https://www.omim.org/")
  url2 <- a("the Comparative Toxicogenomics Database (CTD)", href="http://ctdbase.org/")
  url3<- a("HuGE Navigator",href="https://phgkb.cdc.gov/PHGKB/hNHome.action")
  url4<- a("DisGeNET",href="https://www.disgenet.org/")
  url5<- a("ClinVar",href="https://www.ncbi.nlm.nih.gov/clinvar/")
  url6<- a("Open Targets",href="https://www.opentargets.org/")
  
  output$descript_ADgenes <- renderUI({
      tagList("We assembled multiple data sources to get AD-associated genes.",tags$span(tags$b("First"), style = "color: red;"),", we downloaded ",tags$span(tags$b("54"), style = "color: red;")," genes related to amyloid and ",tags$span(tags$b("27"), style = "color: red;")," genes related to tauopathy from the Supplementary Data of [1], where the author collected experimentally validated genes that satisfied at least one of the following criteria: i) gene validation in large-scale amyloid or tauopathy GWAS studies; ii) in vivo experimental model evidence that knockdown or overexpression of the gene leads to AD-like amyloid or tau pathology. ",tags$span(tags$b("Second"), style = "color: red;"),", we obtained ",tags$span(tags$b("188"), style = "color: red;")," unique late-onset AD common risk genes identified by eight large-scale genetic studies [2-8].",tags$span(tags$b("Third"), style = "color: red;"),", we identified a set of ",tags$span(tags$b("93"), style = "color: red;")," AD-associated genes  curated in at least 2 of 6 following disease gene databases: ",url1,", ",url2,", ",url3,", ",url4,", ",url5,", ",url6,". After moving duplicates, we obtained ",tags$span(tags$b("233"), style = "color: red;")," AD-associated genes in total")
  })
  
  all_AD_gene_withSource <- read.delim("./all_AD_genes_withSource.txt",sep="\t",header=TRUE,row.names = NULL)
  experimental_AD_gene_withSource <- read.delim("./Experimentally_validated_genes.txt",header = TRUE,sep = "\t",row.names = NULL)
  AD_risk_gene_withSource <- read.delim("./GWAS_AD_risk_gene_withSource.txt",header = TRUE,sep="\t",row.names = NULL)
  DB_AD_genes_withSource<- read.delim("./final_AD_genes_DB_withSource.txt",header = TRUE,sep="\t",row.names = NULL)
  all_PPIs<-read.delim("./PPI_net.txt",header=TRUE,sep="\t",row.names=NULL)

  output$AD_geneTable = DT::renderDataTable(
      server = FALSE,
      if(input$ADgene_source==1){
        DT::datatable({
          all_AD_gene_withSource
        },
        extensions = 'Buttons',
        options = list(lengthMenu=list(c(5,15,20),c('5','15','20')),
                       pageLength=10,
                       initComplete = JS(
                         "function(settings, json) {",
                         "$(this.api().table().header()).css({'background-color': 'moccasin', 'color': '1c1b1b'});",
                         "}"),
                       dom = 'Bfrtip',
                       buttons = c('copy', 'csv', 'excel'),
                       columnDefs=list(list(className='dt-center',targets="_all"))
        ),
        selection = 'multiple',
        style = 'bootstrap',
        class = 'cell-border stripe',
        rownames = FALSE,
        colnames = c("Gene","Source"))
      } else if(input$ADgene_source==2){
        DT::datatable({
          experimental_AD_gene_withSource
        },
        extensions = 'Buttons',
        options = list(lengthMenu=list(c(5,15,20),c('5','15','20')),
                       pageLength=10,
                       initComplete = JS(
                         "function(settings, json) {",
                         "$(this.api().table().header()).css({'background-color': 'moccasin', 'color': '1c1b1b'});",
                         "}"),
                       dom = 'Bfrtip',
                       buttons = c('copy', 'csv', 'excel'),
                       columnDefs=list(list(className='dt-center',targets="_all"))
        ),

        style = 'bootstrap',
        class = 'cell-border stripe',
        rownames = FALSE,
        colnames = c("Gene","GeneEntrezID","Protein","Alias","EvidenceDescription","Reference"))
      }
     else if(input$ADgene_source==3){
        DT::datatable({
          AD_risk_gene_withSource
        },
        extensions = 'Buttons',
        options = list(lengthMenu=list(c(5,15,20),c('5','15','20')),
                       pageLength=10,
                       initComplete = JS(
                         "function(settings, json) {",
                         "$(this.api().table().header()).css({'background-color': 'moccasin', 'color': '1c1b1b'});",
                         "}"),
                       dom = 'Bfrtip',
                       buttons = c('copy', 'csv', 'excel'),
                       columnDefs=list(list(className='dt-center',targets="_all"))
        ),
        style = 'bootstrap',
        class = 'cell-border stripe',
        rownames = FALSE,
        colnames = c("Gene","Source"))
     } else if(input$ADgene_source==4){
       DT::datatable({
         DB_AD_genes_withSource
       },
       extensions = 'Buttons',
       options = list(lengthMenu=list(c(5,15,20),c('5','15','20')),
                      pageLength=10,
                      initComplete = JS(
                        "function(settings, json) {",
                        "$(this.api().table().header()).css({'background-color': 'moccasin', 'color': '1c1b1b'});",
                        "}"),
                      dom = 'Bfrtip',
                      buttons = c('copy', 'csv', 'excel'),
                      columnDefs=list(list(className='dt-center',targets="_all"))
       ),
       style = 'bootstrap',
       class = 'cell-border stripe',
       rownames = FALSE,
       colnames = c("Gene","Source"))
     }
  )
  #####-----------------------------
  #for tabPanel5 (Download)
  ####------------------------------
  # Downloadable csv of selected dataset ----
  output$downloadData1 <- downloadHandler(
    filename = function() {
      paste("ranked_drugs_", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(result_data, file, row.names = FALSE)
    }
  )
  
  output$downloadData2 <- downloadHandler(
    filename = function() {
      paste("drug_target_interactions_", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(drug_target_data, file, row.names = FALSE)
    }
  )
  
  output$downloadData3 <- downloadHandler(
    filename = function() {
      paste("AD_associated_genes_", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(all_AD_gene_withSource, file, row.names = FALSE)
    }
  )

  output$downloadData4 <- downloadHandler(
    filename = function() {
      paste("Protein_protein_interactions_", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(all_PPIs, file, row.names = FALSE)
    }
  )
} # server

# Create Shiny object
shinyApp(ui = ui, server = server)