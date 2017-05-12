#### Visualize Adverse Outcome Pathway (AOP) WIKI Data ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##  Set working directory and import key event relationships
library(igraph)
workingDir <-"C://Users/NPollesc/Desktop/GitHub/AOPwiki/"
dataFile <- "data/all-KERs.txt"
source("C://Users/NPollesc/Desktop/GitHub/AOPwiki/R_files/AOP_net_functions.R")
allKERs<-read.table(paste(workingDir, dataFile, sep=""), sep="\t", stringsAsFactors=FALSE, header=TRUE)

##  Identify all unique KEs
allKEs<-c(allKERs[,2],allKERs[,4])
uniqueKEs<-unique(allKEs)
keID<-data.frame(ID=1:length(uniqueKEs),KE=uniqueKEs)

KERs<-cbind(match(allKERs[,2], keID[,2]),match(allKERs[,4], keID[,2]))
AOPgraph<-graph_from_edgelist(KERs, directed=T)

## Create names for vertices in AOPwiki graph
V(AOPgraph)$KE_name<-as.character(keID$KE)
V(AOPgraph)$name<-keID$KE

#Plot the AOP wiki
par(bg="black")
set.seed(1)
plot(AOPgraph, vertex.size=5, edge.arrow.size=.1, edge.width=2)

#### CONNECTED COMPONENTS ANALYSIS ####
V(AOPgraph)$color<-unlist(color.comps(AOPgraph,"strong")$vcol)
E(AOPgraph)$color<-unlist(color.comps(AOPgraph,"strong")$ecol)


## Size change to highlight strongly connected components 
# TASK: Can consider incorporating this as an option in color.comps()
V(AOPgraph)$size<-1
E(AOPgraph)$width<-1
V(AOPgraph)$size[which(V(AOPgraph)$color!="gray")]<-2
E(AOPgraph)$width[which(E(AOPgraph)$color!="gray")]<-2

set.seed(1)
plot(AOPgraph, vertex.size=V(AOPgraph)$size, edge.arrow.size=.1, vertex.label=NA)#, vertex.color="orange",edge.color="gray")

# #### Graph condensation ####
# 
# condense.map(AOPgraph)
# condense.graph(AOPgraph,condense.map(AOPgraph))
# plot(AOPgraph.cond,vertex.size=1, edge.arrow.size=.1, vertex.label=NA)
# 
# is.dag(AOPgraph)
# is.dag(AOPgraph.cond)
# 
# #Attempt a topological sorting of the entire AOPwiki (Hint, it is unreadable :) )
# plot(AOPgraph.cond,vertex.size=2, edge.arrow.size=.1,layout=topo.layout(AOPgraph.cond))

#### CENTRALITY MEASURES FOR THE AOPWIKI ####
#  Which Key Event is most represented? 
names(which(table(KERs)==max(table(KERs))))
V(AOPgraph)$KE_name[72]

# Which key event has the most incident nodes?
sort(degree(AOPgraph))
V(AOPgraph)$KE_name[345]

# Which key event is involved in the most shortest paths between other key events?
sort(betweenness(AOPgraph))
V(AOPgraph)$KE_name[345]

# Which key event is the closest to the rest?
sort(closeness(AOPgraph))
V(AOPgraph)$KE_name[711]


