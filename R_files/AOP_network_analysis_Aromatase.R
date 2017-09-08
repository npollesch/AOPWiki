#### Load proper libraries for AOP Networks and Graph Analyses ####

library("aop")
library(graph)
library(igraph)
library(readr)
library(ggplot2)
library(rgl)
library(plotrix)

workingDir <- "C://Users/NPollesc/Desktop/GitHub/AOPwiki/" #EPA Dir
# workingDir<- "C://Users/Nathan Pollesch/Documents/GitHub/AOPWiki/" #Personal Dir
setwd(workingDir)
ARMimport<- "data/Aromatase network.txt"
THYimport<- "data/Thyroid network.txt"

source(paste(workingDir,"R_files/AOP_net_functions.R",sep="")) #imports custom functions
#KERdata<-read.table(paste(workingDir, KERimport, sep=""), sep="\t", stringsAsFactors=FALSE, header=TRUE)
#KEdata<-read.table(paste(workingDir, KEimport, sep=""), sep="\t", stringsAsFactors=FALSE, header=TRUE)
#KEPdata<-read.table(paste(workingDir, KEplus, sep=""), sep="\t", stringsAsFactors=FALSE, header=TRUE)

#### AOP: Import and format network ####
# Name the AOP network
AOP.name<- "Aromatase"
# Specify import files (loads from working directory)
importnetwork <- "data/aromatase_AOP.cyjs"
importnetattr <- "data/aromatase_AOP_data_shortnames.csv"

#Import network as cytoscape object and turn in aop object
#stea_aop <- convert_cytoscape_to_aop("//aa.ad.epa.gov/ORD/DUL/USERS/NPollesc/Net MyDocuments/projects/Pellston_2017/Ed_Anze/Programs/steatosis/Steatosis_AOP/steatosis.cyjs")
AOP_net<-convert_cytoscape_to_aop(paste(workingDir,importnetwork,sep=""))

#Import additional data exported from cytoscape, including node names and event types
#stea_info <- read_csv("//aa.ad.epa.gov/ORD/DUL/USERS/NPollesc/Net MyDocuments/projects/Pellston_2017/Ed_Anze/Programs/steatosis/Steatosis_AOP/stea_nodes.csv")
AOP_data<- read_csv(paste(workingDir,importnetattr,sep=""))
#convert aop object to graphNEL object
AOP_graph <- convert_aop_to_graph(AOP_net)
#convert graphNEL object to igraph object for analysis using igraph package
AOP<-igraph.from.graphNEL(AOP_graph, name = TRUE, weight = TRUE,
                         unlist.attrs = TRUE)

## ASSIGN NAMES, ATTRIBUTES, and NODE COLORS BY KED

# Assigning names can be tricky, since it depends on how cytoscape exports the network .json file

# V(AOP)names were imported as 'SUID', not 'name', this tells how to match the entries in table_node to the entries in sg graph object by SUID

match(as.character(AOP_data$SUID),V(AOP)$name)

# Assign names to igraph object (note that somehow the node names in V(sg)$name are in descending order, so 'rev' is needed to map actual names properly)
V(AOP)$name<-AOP_data$name[match(V(AOP)$name,as.character(AOP_data$SUID))]
# Assign key event descriptor types (MIE,AO,KE,etc...) And color MIE and AO
V(AOP)$ked<-AOP_data$`Event Type`[match(V(AOP)$name,as.character(AOP_data$name))]

#V(AOP)$KEID<-AOP_data$KEID[match(V(AOP)$name,as.character(AOP_data$name))]
# Assign level of biological organization data from KEPdata sheet in AOP_wiki_analysis file
#V(AOP)$lobo<-KEPdata$Level.of.Organisation[match(V(AOP)$KE_ID,KEPdata$Event)]
#V(AOP)$KE_ID
#V(AOP)$KE_ID[which(is.na(V(AOP)$lobo))]
#match(V(AOP)$KE_ID,KEPdata$Event)

# # Assign KEs with 0 in-degree to be MIE and 0-out-degree to be AO
# indeg<-degree(AOP,mode="in")
# outdeg<-degree(AOP,mode="out")
# indeg
# outdeg
# for (i in 1:length(V(AOP)))
# {if(indeg[i]==0)
#   V(AOP)[i]$KED<-"MIE"
# else V(AOP)[i]$KED<-"KE"
# }
# for (i in 1:length(V(AOP)))
# {if(outdeg[i]==0)
#   V(AOP)[i]$KED<-"AO"
# }
# V(AOP)$name[which(V(AOP)$ked!=V(AOP)$KED)]
# 
# V(AOP)$ked[which(V(AOP)$ked!=V(AOP)$KED)]<-V(AOP)$KED[which(V(AOP)$ked!=V(AOP)$KED)]

E(AOP)$color<-"gray"
V(AOP)$color<-"Yellow"
V(AOP)$color[which(V(AOP)$ked=="MIE")]<-"Green"
V(AOP)$color[which(V(AOP)$ked=="AO")]<-"Red"
# # If one wishes to reshape the vertices, these are the igraph shapes
# names(igraph:::.igraph.shapes)


#### PLOTS NETWORKS ####
lyt=layout.davidson.harel(AOP)

lyt
# creates plot of igraph object with KED coloring
set.seed(3)
plot(AOP, edge.arrow.size=.35, edge.curved=.1, vertex.size=8, main=paste(AOP.name,"AOP Network"),layout=layout.davidson.harel(AOP), vertex.label=paste(strtrim(V(AOP)$name,20),"...",sep=""), vertex.label.color="black",vertex.label.cex=.9, edge.curved=T )
legend('topleft',c("MIE","KE", "AO"), pch=22,
       col="#777777",  pt.bg=c("Green","Yellow","Red"), pt.cex=2, cex=.8, bty="n", ncol=1)

jpeg(filename = paste("AOP_",AOP.name,"_plot.jpeg"),
     width = 800, height = 800, units = "px", pointsize = 12,
     quality = 100)
set.seed(3)
plot(AOP, edge.arrow.size=.35, edge.curved=.1, vertex.size=8, main=paste(AOP.name,"AOP Network"),layout=layout.davidson.harel(AOP), vertex.label.cex=.9, vertex.label.color="black", vertex.label=paste(strtrim(V(AOP)$name,20),"...",sep=""), edge.curved=T )
legend('bottomright',c("MIE","KE", "AO"), pch=22,
       col="#777777",  pt.bg=c("Green","Yellow","Red"), pt.cex=2, cex=.8, bty="n", ncol=1)
dev.off()


#### AOP: Network Form Including Cycles, DAG, and condensation graphs ####
# is.dag() reports TRUE for DAG, FALSE if cycles exist
is.dag(AOP)

# Colors cycles for cyclic graphs (Error in function for when graph is acyclic. Need to fix.)
# dev.off()
# V(AOP)$color<-cycle.color(AOP)[[1]]
# E(AOP)$color<-cycle.color(AOP)[[2]]
# 
# #Find colors used for connected components, so that they can be added to the legend
# which(V(AOP)$color!="Red"&V(AOP)$color!="Yellow"&V(AOP)$color!="Green")
# V(AOP)$color[which(V(AOP)$color!="Red"&V(AOP)$color!="Yellow"&V(AOP)$color!="Green")]
# V(AOP)$color[11]
# 
# set.seed(5)
# plot(AOP,edge.arrow.size=.25, vertex.label=NA,edge.curved=.1, vertex.size=10, main=paste(AOP.name,"AOP Network - Connected Components"),layout=layout.davidson.harel(AOP))
# legend('topleft',c("MIE","KE", "AO","Connected Component"), pch=22,
#        col="#777777",  pt.bg=c("Green","Yellow","Red",V(AOP)$color[11]), pt.cex=2, cex=.8, bty="n", ncol=1)
# 
# jpeg(filename = paste("AOP_",AOP.name,"_CC_plot.jpeg"),
#      width = 800, height = 800, units = "px", pointsize = 12,
#      quality = 75)
# set.seed(3)
# plot(AOP,edge.arrow.size=.25, edge.curved=.1, vertex.size=10, main=paste(AOP.name,"AOP Network - Connected Components"),layout=layout.davidson.harel(AOP))
# legend('bottomleft',c("MIE","KE", "AO","Connected Component"), pch=22,
#        col="#777777",  pt.bg=c("Green","Yellow","Red",V(AOP)$color[11]), pt.cex=2, cex=.8, bty="n", ncol=1)
# dev.off()
# #(B.1) if condensation is desired, create map and condense graph
# 
# condense.map(AOP)
# condense.graph(AOP,condense.map(AOP))
# 
# #(B.2) verify that new graph is DAG and plot
# is.dag(AOP.cond)
# 
# set.seed(36)
# plot(AOP.cond, vertex.label=NA, main=paste(AOP.name,"AOP Network - Condensed"), edge.arrow.size=.25, edge.curved=.1, vertex.size=10,layout=layout.davidson.harel(AOP.cond))
# legend('topleft',c("MIE","KE", "AO","Contracted Node"), pch=22,
#        col="#777777",  pt.bg=c("Green","Yellow","Red","Purple"), pt.cex=2, cex=.8, bty="n", ncol=1)
# 
# 
# jpeg(filename = paste("AOP_",AOP.name,"_cond_plot.jpeg"),
#      width = 800, height = 800, units = "px", pointsize = 12,
#      quality = 75)
# set.seed(3)
# plot(AOP.cond, main=paste(AOP.name,"AOP Network - Condensed"), edge.arrow.size=.25, edge.curved=.1, vertex.size=10,layout=layout.davidson.harel(AOP.cond))
# legend('bottomleft',c("MIE","KE", "AO","Contracted Node"), pch=22,
#        col="#777777",  pt.bg=c("Green","Yellow","Red","Purple"), pt.cex=2, cex=.8, bty="n", ncol=1)
# dev.off()

#### AOP: Topological Sorting of Network ####
#every dag has a topological sorting, so first test if the graph is a DAG
#is.dag(AOP.cond)
is.dag(AOP)

#find topological layout
#topo.layout(AOP.cond)
topo.layout(AOP)

# #plots topologically sorted network using defined layout, plot options are specified for visualization given topological sorting layout
# par(mar=c(5,1,5,5))
# plot(AOP.cond, layout=topo.layout(AOP.cond), edge.arrow.size=.25, edge.curved=1, vertex.size=5, edge.width=1.5, edge.color="Gray", vertex.label.degree=0, vertex.label=NA,vertex.label.dist=.75, main=paste(AOP.name," AOP Network - Topological Sorting"))
# legend('bottomright',c("MIE","KE", "AO","Contracted Node"), pch=22,
#        col="#777777",  pt.bg=c("Green","Yellow","Red","Purple"), pt.cex=2, cex=.8, bty="n", ncol=1)
# 
# par(mar=c(5,1,5,5))
# jpeg(filename = paste("AOP_",AOP.name,"_topo_plot.jpeg"),
#      width = 1000, height = 800, units = "px", pointsize = 12,
#      quality = 75)
# plot(AOP.cond, layout=topo.layout(AOP.cond), edge.arrow.size=.25, edge.curved=1, vertex.size=5, edge.width=1.5, edge.color="Gray", vertex.label.degree=0, vertex.label.dist=.75, main=paste(AOP.name," AOP Network - Topological Sorting"))
# legend('bottomright',c("MIE","KE", "AO","Contracted Node"), pch=22,
#        col="#777777",  pt.bg=c("Green","Yellow","Red","Purple"), pt.cex=2, cex=.8, bty="n", ncol=1)
# dev.off()


#If AOP not AOP.cond
dev.off()
par(mar=c(5,1,5,5))
plot(AOP, layout=topo.layout(AOP), vertex.label.color="black",edge.arrow.size=.25, edge.curved=1, vertex.size=5, edge.width=1.5, edge.color="Gray", vertex.label.degree=0, vertex.label.dist=1, main=paste(AOP.name," AOP Network - Topological Ordering"),vertex.label=paste(strtrim(V(AOP)$name,35),"...",sep=""))
legend('bottomright',c("MIE","KE", "AO"), pch=22,
       col="#777777",  pt.bg=c("Green","Yellow","Red"), pt.cex=2, cex=.8, bty="n", ncol=1)

#save as JPEG
par(mar=c(5,1,5,5))
jpeg(filename = paste("AOP_",AOP.name,"_topo_plot.jpeg"),
     width = 1000, height = 800, units = "px", pointsize = 12,
     quality = 75)
plot(AOP, layout=topo.layout(AOP), vertex.label.color="black",edge.arrow.size=.25, edge.curved=1, vertex.size=5, edge.width=1.5, edge.color="Gray", vertex.label.degree=0, vertex.label.dist=1, main=paste(AOP.name," AOP Network - Topological Sorting"), vertex.label=paste(strtrim(V(AOP)$name,35),"...",sep="") )
legend('bottomright',c("MIE","KE", "AO"), pch=22,
       col="#777777",  pt.bg=c("Green","Yellow","Red"), pt.cex=2, cex=.8, bty="n", ncol=1)
dev.off()

# Topolayout with PO coloration (Must run MIE to AO path occurrence analysis in later section)

# Create a list of ked designations to plot near 

# V(AOP)$po_col
dev.off()
par(mar=c(5,1,5,5))
plot(AOP, layout=topo.layout(AOP), vertex.label.color="black", vertex.color=V(AOP)$po_col, edge.arrow.size=.25, edge.curved=1, vertex.size=5, edge.width=1.5, edge.color="Gray", vertex.label.degree=0, main=paste(AOP.name," AOP Network - Topological Ordering, MIE to AO Path Occurrence Coloration"),vertex.label=paste(V(AOP)$ked,": ",paste(strtrim(V(AOP)$name,30),"...",sep=""),sep=""), vertex.label.dist=1.1)

#save as JPEG
par(mar=c(5,1,5,5))
jpeg(filename = paste("AOP_",AOP.name,"_topo_plot_PO.jpeg"),
     width = 1000, height = 800, units = "px", pointsize = 12,
     quality = 75)
par(mar=c(5,1,5,5))
plot(AOP, layout=topo.layout(AOP), vertex.label.color="black", vertex.color=V(AOP)$po_col, edge.arrow.size=.25, edge.curved=1, vertex.size=5, edge.width=1.5, edge.color="Gray", vertex.label.degree=0, main=paste(AOP.name," AOP Network - Topological Ordering, MIE to AO Path Occurrence Coloration"),vertex.label=paste(V(AOP)$ked,": ",paste(strtrim(V(AOP)$name,30),"...",sep=""),sep=""), vertex.label.dist=1.1)
dev.off()


##Notes on topological sorting of AOP networks:
  #This sorting can be used to prioritize measurement and have a quick visual respresentation of potential up-stream effects.  I.e. if you are towards the beginning, only few KE can affect you, but towards the end, many more could (not that they necessarily do)
  #This sorting also lines up logically with MIE -> KE -> AO

#### AOP: Path count function and visualizations ####

#!TASK: NEED TO FIX THIS FUNCTION - PERHAPS SOME PATHS DON"T EXIST BETWEEN MIE AND ALL AOS? ####
# V(AOP.cond)$aoppn<-aop.paths(AOP.cond, nor=TRUE)
# V(AOP.cond)$aopp<-aop.paths(AOP.cond, nor=FALSE)

V(AOP)$aoppn<-aop.paths(AOP, nor=TRUE)
V(AOP)$aopp<-aop.paths(AOP, nor=FALSE)

## Node coloring option based on no of path occurances from MIE to AO
# heatgrad=rev(heat.colors(n=max(aop.paths(AOP.cond, norm=FALSE))))
# heatgrad<-append(heatgrad, "#FFFFFF", after=0)
# 
# wbpal=colorRampPalette(c("white","blue"))(n=max(aop.paths(AOP.cond, norm=FALSE))+1)


heatgrad=rev(heat.colors(n=max(aop.paths(AOP, norm=FALSE))))
heatgrad<-append(heatgrad, "#FFFFFF", after=0) #in the instance that a KE is not included in paths between MIE and AO...

wbpal=colorRampPalette(c("white","blue"))(n=max(aop.paths(AOP, norm=FALSE))+1)



## shows the color scheme with associated values
# pie(rep(1,max(aop.paths(AOP.cond,norm=FALSE))+1),col=rev(heat.colors(n=max(aop.paths(AOP.cond, norm=FALSE))+1)))

## assign color attribute of vertices
#V(AOP.cond)$po_col<-wbpal[aop.paths(AOP.cond, norm=FALSE)+1]# applies color scheme to nodes


V(AOP)$po_col<-wbpal[aop.paths(AOP, norm=FALSE)+1]# applies color scheme to nodes

# ### Topologically sorted plot with path occurance coloration
# #dev.off()
# plot(AOP.cond, layout=topo.layout(AOP.cond), edge.arrow.size=.25, edge.curved=1, vertex.size=6, edge.width=1.5, edge.color="Gray", vertex.label.degree=0, vertex.label.dist=.75, main="Steatosis AOP - Path Occurrence between MIE and AO")
# color.legend(-1.,.5,-1.1,1, grad="y",legend=c(paste("Low (",min(V(AOP.cond)$aopp),")"),paste("High (",max(V(AOP.cond)$aopp),")")), rect.col=heatgrad[c(trunc(seq(0,max(V(AOP.cond)$aopp),length.out=10)),0)])

# AOP Networks with path occurrence colorings 

# for AOP.cond
# dev.off()
# set.seed(3)
# plot(AOP.cond, main=paste(AOP.name,"AOP Network - Path Occurrences"), edge.arrow.size=.5, edge.curved=.1, vertex.size=10,layout=layout.davidson.harel(AOP.cond),vertex.color=V(AOP)$po_col,vertex.label=paste(strtrim(V(AOP)$name,20),"...",sep=""))
# color.legend(1.,.5,1.1,1, grad="y",legend=c(paste("Low (",min(V(AOP.cond)$aopp),")"),paste("High (",max(V(AOP.cond)$aopp),")")), rect.col=wbbpal[c(trunc(seq(0,max(V(AOP.cond)$aopp),length.out=10)),0)])
# 
# jpeg(filename = paste("AOP_",AOP.name,"_PO_plot.jpeg"),
#      width = 900, height = 800, units = "px", pointsize = 12,
#      quality = 75)
# set.seed(3)
# plot(AOP.cond, main=paste(AOP.name,"AOP Network - Path Occurrences"), edge.arrow.size=.5, edge.curved=.1, vertex.size=10,layout=layout.davidson.harel(AOP.cond),vertex.color=V(AOP)$po_col,vertex.label=paste(strtrim(V(AOP)$name,20),"...",sep=""))
# color.legend(1.,.5,1.1,1, grad="y",legend=c(paste("Low (",min(V(AOP.cond)$aopp),")"),paste("High (",max(V(AOP.cond)$aopp),")")), rect.col=wbpal[c(trunc(seq(0,max(V(AOP.cond)$aopp),length.out=10)),0)])
# dev.off()



dev.off()
set.seed(3)
plot(AOP, main=paste(AOP.name,"AOP Network - Path Occurrences"), vertex.label.cex=.9, edge.arrow.size=.5, edge.curved=.1, vertex.size=10,layout=layout.davidson.harel(AOP),vertex.color=V(AOP)$po_col,vertex.label=paste(strtrim(V(AOP)$name,20),"...",sep=""), vertex.label.color="Black")
color.legend(1.,-.50,1.1,-1, grad="y",legend=c(paste("Low (",min(V(AOP)$aopp),")"),paste("High (",max(V(AOP)$aopp),")")), rect.col=wbpal[c(trunc(seq(0,max(V(AOP)$aopp),length.out=10)),0)])

jpeg(filename = paste("AOP_",AOP.name,"_PO_plot.jpeg"),
     width = 800, height = 800, units = "px", pointsize = 12,
     quality = 100)
set.seed(3)
plot(AOP, main=paste(AOP.name,"AOP Network - Path Occurrences"), vertex.label.cex=.9, edge.arrow.size=.5, edge.curved=.1, vertex.size=10,layout=layout.davidson.harel(AOP),vertex.color=V(AOP)$po_col, vertex.label.color="Black")
color.legend(1.,-.50,1.1,-1, grad="y",legend=c(paste("Low (",min(V(AOP)$aopp),")"),paste("High (",max(V(AOP)$aopp),")")), rect.col=wbpal[c(trunc(seq(0,max(V(AOP)$aopp),length.out=10)),0)])
dev.off()

### BARPLOT OF aop.paths() data
# par(mfrow=c(1,1),las=2, mar=c(5,15,4,3)) # make label text perpendicular to axis and increase margin
# #V(sg)$color<-heatgrad[node_counts$count+1]# applies color scheme to nodes 
# aop.paths.bp<-barplot(aop.paths(AOP.cond), xlim=c(0,max(aop.paths(AOP.cond))+1), main="Steatosis AOP - Number of occurrences in paths between MIEs and AOs", horiz=TRUE, names.arg=paste(strtrim(V(AOP)$name,35),"...",sep=""), col=V(AOP.cond)$po_col, axes=FALSE, cex.names=.75, cex.main=.8, adj=0)
# 
# #bpall<-barplot(node_counts$count, xlim=c(0,max(counts)+1), main="Number of occurances in paths between MIEs and AOs", horiz=TRUE, names.arg=names(counts), col=V(sg)$color, axes=FALSE, cex.names=.75, cex.main=.8, adj=0, beside=FALSE)
# axis(1,at=seq(0,45,5))
# text(1, aop.paths.bp, aop.paths(AOP.cond,norm=FALSE), cex=1,pos=4) #adds values to bars
# text(.25, aop.paths.bp, V(AOP.cond)$ked ,cex=.5,pos=4 ) #adds OE value text next to bars

## BARPLOT FOR AOP not AOP.cond
dev.off()
par(mfrow=c(1,1),las=2, mar=c(5,15,4,3)) # make label text perpendicular to axis and increase margin
#V(sg)$color<-heatgrad[node_counts$count+1]# applies color scheme to nodes 
aop.paths.bp<-barplot(aop.paths(AOP), xlim=c(0,max(aop.paths(AOP))+1), main=paste(AOP.name,"AOP - Number of occurrences in paths between MIEs and AOs"), horiz=TRUE, names.arg=paste(strtrim(V(AOP)$name,35),"...",sep=""), col=V(AOP)$po_col, axes=FALSE, cex.names=.75, cex.main=.8, adj=0)

#bpall<-barplot(node_counts$count, xlim=c(0,max(counts)+1), main="Number of occurances in paths between MIEs and AOs", horiz=TRUE, names.arg=names(counts), col=V(sg)$color, axes=FALSE, cex.names=.75, cex.main=.8, adj=0, beside=FALSE)
axis(1,at=seq(0,60,5))
text(1, aop.paths.bp, aop.paths(AOP,norm=FALSE), cex=.9,pos=4) #adds values to bars
text(-.1, aop.paths.bp, V(AOP)$ked ,cex=.5,pos=4 ) #adds OE value text next to bars

jpeg(filename = paste("AOP_",AOP.name,"_PO_barplot.jpeg"),
     width = 900, height = 800, units = "px", pointsize = 12,
     quality = 100)
par(mfrow=c(1,1),las=2, mar=c(5,15,4,3)) # make label text perpendicular to axis and increase margin
#V(sg)$color<-heatgrad[node_counts$count+1]# applies color scheme to nodes 
aop.paths.bp<-barplot(aop.paths(AOP), xlim=c(0,max(aop.paths(AOP))+1), main=paste(AOP.name,"AOP - Number of occurrences in paths between MIEs and AOs"), horiz=TRUE, names.arg=paste(strtrim(V(AOP)$name,35),"...",sep=""), col=V(AOP)$po_col, axes=FALSE, cex.names=.75, cex.main=.8, adj=0)

#bpall<-barplot(node_counts$count, xlim=c(0,max(counts)+1), main="Number of occurances in paths between MIEs and AOs", horiz=TRUE, names.arg=names(counts), col=V(sg)$color, axes=FALSE, cex.names=.75, cex.main=.8, adj=0, beside=FALSE)
axis(1,at=seq(0,60,5))
text(1, aop.paths.bp, aop.paths(AOP,norm=FALSE), cex=.9,pos=4) #adds values to bars
text(-.1, aop.paths.bp, V(AOP)$ked ,cex=.5,pos=4 ) #adds OE value text next to bars

dev.off()

#### Robustness Measures for MIE and AO ####

#### ~~ Edge connectivity ####
# This function produces a triplet of (MIE, AO, EdgeConnectivity) for each
# comibination of MIE to AO. names=TRUE returns the MIE to AO pairs as names
# otherwise a matrix is return with vertex indices for graph supplied
ecl<-aop.edge.connectivity(AOP,names=F)
ecln<-aop.edge.connectivity(AOP,names=T)

ecln
#Sorted by highest edge connectivity, i.e. most robustly connected MIE to AO in network
ecl[order(ecl[,3],decreasing=T),]

#Creates a new matrix that includes names
 # ecln[,1]<-V(AOP)$name[ecl[,1]]
 # ecln[,2]<-V(AOP)$name[ecl[,2]]
 # ecln[,3]<-ecl[,3]

#### ~~ Simple Path Node Redundancy and Simple Path Visualization ####
# This is a custom analysis that calculates the number of nodes involved 
# in simple paths between nodes that can be removed while still maintaining
# a path between the source and target node. It is calculated as the
# difference in magnitude between the union of all nodes in simple paths
# and the nodes in the(a) shortest path.
simple.path.redundancy(AOP,51,47)
simple.path.redundancy(AOP,51,42)

# This is a custom visualization that can be used to plot all simple paths
# between MIE and AOs.

# USE locations from simple.path.coloring to specify a color
# This usage allows for adding in multiple path colors
# Note: loc=T color assignments will overwrite previous assignments

E(AOP)$asp_clr<-"gray"
E(AOP)$asp_clr[simple.path.coloring(AOP,51,47)]<-"purple"

set.seed(3)
plot(AOP,edge.color=E(AOP)$asp_clr,edge.arrow.size=.2,vertex.size=5,vertex.label=V(AOP))

# USE non-location option of simple.path.coloring to color graph
# Note: This loc=F option only allows for a single path to be highlighted
# Note (cont): highlight color is specified as clr="color"
E(AOP)$asp_clr<-simple.path.coloring(AOP,11,14,loc=F,clr="blue")

set.seed(3)
plot(AOP,edge.color=E(AOP)$asp_clr,edge.arrow.size=.2,vertex.size=5,vertex.label=V(AOP))




#### LOBO Layout Visualization ####
lobo_list=c("Molecular","Cellular","Tissue","Organ","Individual","Population","") #creates an ordering of biological organization
V(AOP)$lobo[which(is.na(V(AOP)$lobo))]<-"" #NAs are not recognized by the lobo.layout function, only those created in the ordering lobo_list above.
V(AOP)$lobo_o<-match(V(AOP)$lobo,lobo_list)
tcols=rainbow(length(unique(V(AOP)$lobo))) #creates a color scheme for visualization
V(AOP)$lobo_col<-tcols[V(AOP)$lobo_o]

V(AOP)$lobo_col
lobo.layout(AOP)
plot(AOP, layout=lobo.layout(AOP),vertex.size=5,  edge.curved=.3, edge.color="gray", edge.arrow.size=.1, vertex.label=V(AOP)$KE_ID, vertex.color=V(AOP)$lobo_col)



