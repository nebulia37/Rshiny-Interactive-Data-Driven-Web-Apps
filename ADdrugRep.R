# Load R packages
library(shiny)
library(shinythemes)
library(DT)
library(dplyr)

setwd("/fs/ess/PCON0041/shunian/AD_drugRep_project_202209/Code/Rshiny/")
# Define UI
ui <- fluidPage(
          theme = shinytheme("united"),
          titlePanel("Drug repurposing for Alzheimer's Disease based on DeepWalk"),
          navbarPage("",
              tabPanel(icon("home"),
                  fluidRow(column(tags$img(src="Antioquia.png",width="200px",height="260px"),width=2),
                  column(
                    br(),
                    p("Through this application, it is intended to represent the result of our drug repurposing for Alzheimei's Disease(AD) project. 
                      We constructed a graph  that contains protein-protein interactions, drug-target interactions and AD-associated genes. 
                      Then, we apply the DeepWalk method on the graph to learn node representations. 
                      The model then scores and ranks drug candidates based on their proximity to AD-related genes within the network through the learned embeddings. Users can search for the rank information of a certain drug in 'Search drug rank' section",
                      style="text-align:justify;color:black;background-color:lavender;padding:15px;border-radius:10px"),
                    br(),
                    p("Besides the result of drug repurposing, we also represent the reliable drug target interaction assembled from four databases and AD associated gene list we obtained from multiple sources, the detail of which are in the 'Drug target interactions' and 'AD associated genes' section respectively.",style="text-align:justify;color:black;background-color:papayawhip;padding:15px;border-radius:10px"),
                      width=8),
                  column(
                    br(),
                    tags$img(src="Gobernacion.png",width="200px",height="130px"),
                    br(),
                    br(),
                    p("For more information please check the",em("Anuario Estadístico de Antioquia's"),"page clicking",
                    br(),
                    a(href="http://www.antioquiadatos.gov.co/index.php/anuario-estadistico-de-antioquia-2016", "Here",target="_blank"),style="text-align:center;color:black"),
                    width=2)),
                    hr(),
                    tags$style(".fa-database {color:#E87722}"),
                    h3(p(em("Top50 ranked drugs repurposed for AD "),icon("database",lib = "font-awesome"),style="color:black;text-align:center")),
                    fluidRow(column(DT::dataTableOutput("DrugData"),
                                     width = 12)),
                                    hr(),
                                    p(em("Developed by"),br("Shunian Xiang"),style="text-align:center; font-family: times")), # tabPanel 1 for home
              tabPanel("Search drug score", 
                    sidebarPanel(
                      tags$h3("Input drug name/drugbank ID for search:"),
                      radioButtons("drugname_type", "",c("drug name","drugbank ID"),selected = "drug name"),
                      textInput("drug_name_id", "Drug Name or ID:",value="galantamine"),
                      actionButton("submitbutton_1", "Submit", class = "btn btn-primary")), # sidebarPanel
                    mainPanel(
                      h2("Get the rank of a drug by searching"),
                      h4("The higher the rank, the more chance that the drug may have the potential to treat AD"),
                      verbatimTextOutput("txtout_2"),
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
                      textInput("search_name_id", "Drug Name or Gene Symbol:",value = "tacrine"),
                      actionButton("submitbutton_2", "Submit", class = "btn btn-primary")), # sidebarPanel
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
                       downloadButton("downloadData3", "Click to Save")))
)
# Define server function  
server <- function(input, output) {
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
        if(input$drugname_type=="drug name"){
          single_drug_data = result_data %>% filter(Drug==input$drug_name_id)
        }else if(input$drugname_type=="drugbank ID"){
          single_drug_data = result_data %>% filter(NodeID==input$drug_name_id)
        }
        single_drug_data
    })
  
    output$ui_panel2 <- renderUI({
        if(input$submitbutton_1>0 && nrow(datasetInput_2()) == 0)
            return("The input is not in the drug list")
    
        DT::dataTableOutput("single_drug_table")
        hr()
        # output the network bewtween the drug and AD node
        cyjShinyOutput('cyjShiny')
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
    # output network 
    StyleFile = "stylefile.js"
    loadStyleFile(StyleFile)
    ## load RData contains list of dataframes of edges and nodes info for each drug from the drug to AD 
    load("network_forEachDrug.RData")
    
    
    networkInput <- reactive({ 
      if(input$drugname_type=="drug name"){
        single_drug_name = input$drug_name_id
      }else if(input$drugname_type=="drugbank ID"){
        single_drug_name = result_data$drug_name[match(input$drug_name_id,result_data$drugbank_id)]
      }
      
      tbl.nodes <- network_forEachDrug[single_drug_name][0]
      tbl.edges <- network_forEachDrug[single_drug_name][1]
      graph.json <- toJSON(dataFramesToJSON(tbl.edges, tbl.nodes), auto_unbox=TRUE)
      
      graph.json
    })
    
    output$cyjShiny <- renderCyjShiny({
      cyjShiny(graph=networkInput(), layoutName="cola")
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
      tagList("We assembled drug-target interactions and bioactivity data from four commonly used databases:",url1_drug, "(v31),",url2_drug, " (downloaded in November 2022), ", url3_drug," (downloaded in November 2022) and", url4_drug, " (downloaded in November 2022). Of the interactions present in these databases, we retained only those that satisfied the following criteria: (1) binding affinities, including Ki, Kd, IC50, or EC50 ≤10 μM (10,000 nM). (2) gene targets and their respective proteins must have a unique UniProt accession number. (3) protein targets be marked as ‘reviewed’ in the UniProt database22. (4) protein targets found in Homo Sapiens. As a result of these criteria, 10,701 drug-gene interactions between 1,591 US FDA-approved drugs and 1,254 unique genes were obtained and integrated with our PPI network")
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
      if(input$submitbutton_2>0 && nrow(datasetInput_3()) == 0){
        print("The input is not in our list")
      }else{
        DT::dataTableOutput("drug_targetTable")
      }
    })
    
    output$drug_targetTable = DT::renderDataTable(
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
      tagList("We assembled multiple data sources to get AD-associated genes. First, we downloaded 54 genes related to amyloid and 27 genes related to tauopathy from the Supplementary Data of [1], where the author collected experimentally validated genes that satisfied at least one of the following criteria: i) gene validation in large-scale amyloid or tauopathy GWAS studies; ii) in vivo experimental model evidence that knockdown or overexpression of the gene leads to AD-like amyloid or tau pathology. Second, we obtained 188 unique late-onset AD common risk genes identified by eight large-scale genetic studies [2-8]. Third, we identified a set of 93 AD-associated genes  curated in at least 2 of 6 following disease gene databases: ",url1,", ",url2,", ",url3,", ",url4,", ",url5,", ",url6,". After moving duplicates, we obtained 233 AD-associated genes in total")
  })
  
  all_AD_gene_withSource <- read.delim("./all_AD_genes_withSource.txt",sep="\t",header=TRUE,row.names = NULL)
  experimental_AD_gene_withSource <- read.delim("./Experimentally_validated_genes.txt",header = TRUE,sep = "\t",row.names = NULL)
  AD_risk_gene_withSource <- read.delim("./GWAS_AD_risk_gene_withSource.txt",header = TRUE,sep="\t",row.names = NULL)
  DB_AD_genes_withSource<- read.delim("./final_AD_genes_DB_withSource.txt",header = TRUE,sep="\t",row.names = NULL)
  
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
} # server

# Create Shiny object
shinyApp(ui = ui, server = server)