###--------create network for each drug---------
edge_info <- read.csv(file = "/users/PCON0041/shunian/shunian/AD_drugRep_project_202209/Data/Intermediate_data/Graphs_and_embeddings/new.all_drugName_nodeID.edge.txt",header = T)
index_info <- read.csv(file = "/users/PCON0041/shunian/shunian/AD_drugRep_project_202209/Data/Intermediate_data/Graphs_and_embeddings/new.all_drugName_nodeID.index.txt",header = T)
node_id_name <- read.delim(file = "/users/PCON0041/shunian/shunian/AD_drugRep_project_202209/Data/Intermediate_data/Graphs_and_embeddings/New_graph/new.Final_NodeID_Name.txt",sep="\t",header=T)

index_info$drugName = node_id_name$Name[match(index_info$drugNodeID,node_id_name$NodeID)]
edge_info$source_name = node_id_name$Name[match(edge_info$source,node_id_name$NodeID)]
edge_info$target_name = node_id_name$Name[match(edge_info$target,node_id_name$NodeID)]
edge_info$interaction = rep("interact",dim(edge_info)[1])

network_forEachDrug = list()
edges = list()
network_forEachDrug[[1]] = list()
network_forEachDrug[[2]]=data.frame()
## network_forEachDrug[[1]] contains edges info
for(i in 1:dim(index_info)[1]){
  if(is.na(index_info$endLoc[i])){
    edges[[i]] = NA
    index_info$endLoc[i]=index_info$endLoc[i-1]
  }else{
    if(i == 1){
      edges[[i]] = edge_info[1:index_info$endLoc[i],3:5]
    }else{
      edges[[i]] = edge_info[(index_info$endLoc[i-1]+1):index_info$endLoc[i],3:5]
      #print(edges[[i]])
    }
  }
}
names(edges) <- index_info$drugName
network_forEachDrug[[1]] = edges

## network_forEachDrug[[2]] contains nodes info
nodes <- node_id_name%>% mutate(
  Type = case_when(
    NodeID==114483835 ~ "AD",
    NodeID>114483835 ~ "drug",
    NodeID<114483835 ~ "gene"
  )) %>% select(Name,Type)
names(nodes) = c("id","type")
network_forEachDrug[[2]] = nodes
save(network_forEachDrug,file="network_forEachDrug.RData")
###-------create json file for style in cytoscape
library("jsonlite")
# {"selector":"node", "css": {
#"text-valign":"center",
#"text-halign":"center",
#"background-color": "lightgreen",
#"border-color": "black",
#"shape": "ellipse",
#"width": "mapData(count, 0, 300, 30, 80)",
#"height": "mapData(count, 0, 300, 30, 80)",
#"content": "data(id)",
#"border-width": "1px"
#}},
ls1=list(selector="node[type='drug']:",
        css=list("text-valign"="center",
                      "text-halign"="center",
                      "background-color"="lightblue",
                      "border-color"="black",
                      shape= "round-rectangle",
                      width = "40px",
                      height = "40px",
                      content="data(id)",
                      "font-size" = "14px"))
ls2=list(selector="node[type='AD']:",
         css=list("text-lign"="center",
                  "text-halign"="center",
                  "background-color"="red",
                  "border-color"="black",
                  shape= "hexagon",
                  width = "40px",
                  height = "40px",
                  content = "data(id)",
                  "font-size" = "14px"))
ls3=list(selector="node[type='gene']:",
         css=list("text-valign"="center",
                  "text-halign"="center",
                  "background-color"="lightpurple",
                  "border-color"="black",
                  shape= "circle",
                  content = "data(id)",
                  width = "40px",
                  height = "40px",
                  "font-size" = "14px"))

ls4=list(selector="edge",
         css=list("line-color" = "maroon","curve-style" = "bezier"))

ADrep_style = list(ls1,ls2,ls3,ls4)
ADrep_style_json = jsonlite::toJSON(ADrep_style, auto_unbox = TRUE) %>% jsonlite::prettify()
jsonlite::write_json(ADrep_style_json, "ADrep_style.js", pretty = T)
write_json(ADrep_style_json, file="ADrep_style.json")

