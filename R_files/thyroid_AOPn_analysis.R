#### SECTION 1: LOAD PACKAGES, SET DIRECTORIES, AND LOAD CUSTOM AOP FUNCTIONS ####

## Load packages
library("aop") # Used for converting cytoscape files to graph objects
library(graph) 
library(igraph)
library(readr)
library(plotrix)

workingDir <-"C://Users/NPollesc/Desktop/GitHub/AOPwiki/" #EPA Dir

## Set working directory
#getwd() # View current working directory
#workingDir <-"" #specify the working directory.

setwd(workingDir)

## Load custom AOP function file
source(paste(workingDir,"R_files/AOPn_supplementary_functions.R",sep="")) #imports custom functions
source(paste(workingDir,"AOP_net_functions.R",sep="")) #imports custom functions, AOP_net_functions.R needs to be in working directory


#### SECTION 2: IMPORT AND FORMAT AOP NETWORK FILE ####
# Name the AOP network (you specify the name, it is used in titles for plots and output files)
AOP.name<-"Thyroid"
# Specify import files (loads from working directory)
importnetwork<-"data/thyroid_AOP.cyjs"
importnetattr<-"data/thyroid_AOP_data.csv"

#Import network as cytoscape object and turn in aop object
AOP_net<-convert_cytoscape_to_aop(paste(workingDir,importnetwork,sep=""))
#Import additional data exported from cytoscape, including node names and event types
AOP_data<- read_csv(paste(workingDir,importnetattr,sep=""))
#convert aop object to graphNEL object
AOP_graph <- convert_aop_to_graph(AOP_net)
#convert graphNEL object to igraph object for analysis using igraph package
AOP<-igraph.from.graphNEL(AOP_graph, name = TRUE, weight = TRUE,
                         unlist.attrs = TRUE)

#### Section 3: ASSIGN NAMES, ATTRIBUTES, AND NODE COLORS BY KED ####

## Assigning names can be tricky, since it depends on how cytoscape exports the network .json file

## Assign names to igraph object. The SUID attribute of the .json file is used to match the $name attributed created by default in section 1
V(AOP)$name<-AOP_data$name[match(V(AOP)$name,as.character(AOP_data$SUID))]

## Assign key event designator types (MIE,AO,KE,etc...) and the AOPwiki Key Event ID as attributes
V(AOP)$ked<-AOP_data$`Event Type`[match(V(AOP)$name,as.character(AOP_data$name))]
V(AOP)$KE_ID<-AOP_data$KE_ID[match(V(AOP)$name,as.character(AOP_data$name))]

## Remove comments to view the attributes just assigned
# V(AOP)$ked 
# V(AOP)$KE_ID 

## Assigns edge color as gray and node colors according to key event designator (MIE,KE,AO)
E(AOP)$color<-"gray"
V(AOP)$color<-"Yellow"
V(AOP)$color[which(V(AOP)$ked=="MIE")]<-"Green"
V(AOP)$color[which(V(AOP)$ked=="AO")]<-"Red"

#### Section 4: NETWORK PLOTS ####

## Creates and export plot of AOP network with KED coloring
set.seed(2)
plot(AOP, edge.arrow.size=.25, edge.curved=.1, vertex.size=10, main=paste(AOP.name,"AOP Network"),layout=layout.davidson.harel(AOP), vertex.label.cex=.75, vertex.label=NA, edge.curved=T )
legend('topleft',c("MIE","KE", "AO"), pch=22,
       col="#777777",  pt.bg=c("Green","Yellow","Red"), pt.cex=2, cex=.8, bty="n", ncol=1)
# Set jpeg output attributes
jpeg(filename = paste("AOP_",AOP.name,"_plot.jpeg"), 
     width = 800, height = 800, units = "px", pointsize = 12,
     quality = 75)
set.seed(3)
plot(AOP, edge.arrow.size=.25, edge.curved=.1, vertex.size=10, main=paste(AOP.name,"AOP Network"),layout=layout.davidson.harel(AOP), vertex.label=NA)
legend('bottomleft',c("MIE","KE", "AO"), pch=22,
       col="#777777",  pt.bg=c("Green","Yellow","Red"), pt.cex=2, cex=.8, bty="n", ncol=1)
dev.off()

## Create and export plot of AOP network with topological sorting
is.dag(AOP) # is.dag() reports TRUE for DAG, FALSE if cycles exist

## Plot network with topological sort
dev.off()
par(mar=c(5,1,5,5)) #changes plot margin to accomodate topological sort
plot(AOP, layout=topo.layout(AOP), edge.arrow.size=.25, edge.curved=1, vertex.size=5, edge.width=1.5, edge.color="Gray", vertex.label.degree=0, vertex.label.dist=.75, main=paste(AOP.name," AOP Network - Topological Ordering"), vertex.label=paste(strtrim(V(AOP)$name,35),"...",sep=""))
legend('bottomright',c("MIE","KE", "AO"), pch=22,
       col="#777777",  pt.bg=c("Green","Yellow","Red"), pt.cex=2, cex=.8, bty="n", ncol=1)

## Export plot as JPEG
par(mar=c(5,1,5,5))
jpeg(filename = paste("AOP_",AOP.name,"_topo_plot_NP.jpeg"),
     width = 1000, height = 800, units = "px", pointsize = 12,
     quality = 75)
plot(AOP, layout=topo.layout(AOP), edge.arrow.size=.25, edge.curved=1, vertex.size=5, edge.width=1.5, edge.color="Gray", vertex.label.degree=0, vertex.label.dist=.75, main=paste(AOP.name," AOP Network - Topological Sorting"), vertex.label=paste(strtrim(V(AOP)$name,35),"...",sep="") )
legend('bottomright',c("MIE","KE", "AO"), pch=22,
       col="#777777",  pt.bg=c("Green","Yellow","Red"), pt.cex=2, cex=.8, bty="n", ncol=1)
dev.off()


#### Section 5: AOP PATH COUNT FUNCTION AND VISUALIZATIONS ####

## Applies custom AOP path count function aop.paths() included in AOPn_supplementary_functions.R file
V(AOP)$aopp<-aop.paths(AOP, nor=FALSE) 

## Assign AOP path occurrence values to colors and colors nodes
wbpal=colorRampPalette(c("white","blue"))(n=max(aop.paths(AOP, norm=FALSE))+1)

V(AOP)$po_col<-wbpal[aop.paths(AOP, norm=FALSE)+1] # applies color scheme to nodes

## Plot AOP network with path occurrence coloring
dev.off()
set.seed(3)
plot(AOP, main=paste(AOP.name,"AOP Network - Path Occurrences"), edge.arrow.size=.5, edge.curved=.1, vertex.size=10,layout=layout.davidson.harel(AOP),vertex.color=V(AOP)$po_col,vertex.label=paste(strtrim(V(AOP)$name,20),"...",sep=""), vertex.label.color="Black")
color.legend(1.,.5,1.1,1, grad="y",legend=c(paste("Low (",min(V(AOP)$aopp),")"),paste("High (",max(V(AOP)$aopp),")")), rect.col=wbpal[c(trunc(seq(0,max(V(AOP)$aopp),length.out=10)),0)])

jpeg(filename = paste("AOP_",AOP.name,"_PO_plot.jpeg"),
     width = 900, height = 800, units = "px", pointsize = 12,
     quality = 75)
set.seed(3)
plot(AOP, main=paste(AOP.name,"AOP Network - Path Occurrences"), edge.arrow.size=.5, edge.curved=.1, vertex.size=10,layout=layout.davidson.harel(AOP),vertex.color=V(AOP)$po_col, vertex.label.color="Black")
color.legend(1.,.5,1.1,1, grad="y",legend=c(paste("Low (",min(V(AOP)$aopp),")"),paste("High (",max(V(AOP)$aopp),")")), rect.col=wbpal[c(trunc(seq(0,max(V(AOP)$aopp),length.out=10)),0)])
dev.off()

## Barplot of AOP network with coloring based on path occurrence
dev.off()
par(mfrow=c(1,1),las=2, mar=c(5,15,4,3)) # make label text perpendicular to axis and increase margin

aop.paths.bp<-barplot(aop.paths(AOP), xlim=c(0,max(aop.paths(AOP))+1), main=paste(AOP.name,"AOP - Number of occurrences in paths between MIEs and AOs"), horiz=TRUE, names.arg=paste(strtrim(V(AOP)$name,35),"...",sep=""), col=V(AOP)$po_col, axes=FALSE, cex.names=.75, cex.main=.8, adj=0)

#bpall<-barplot(node_counts$count, xlim=c(0,max(counts)+1), main="Number of occurances in paths between MIEs and AOs", horiz=TRUE, names.arg=names(counts), col=V(sg)$color, axes=FALSE, cex.names=.75, cex.main=.8, adj=0, beside=FALSE)
axis(1,at=seq(0,45,5))
text(1, aop.paths.bp, aop.paths(AOP,norm=FALSE), cex=1,pos=4) #adds values to bars
text(.25, aop.paths.bp, V(AOP)$ked ,cex=.5,pos=4 ) #adds OE value text next to bars

jpeg(filename = paste("AOP_",AOP.name,"_PO_barplot.jpeg"),
     width = 900, height = 800, units = "px", pointsize = 12,
     quality = 100)
par(mfrow=c(1,1),las=2, mar=c(5,15,4,3))
aop.paths.bp<-barplot(aop.paths(AOP), xlim=c(0,max(aop.paths(AOP))+1), main=paste(AOP.name,"AOP - Number of occurrences in paths between MIEs and AOs"), horiz=TRUE, names.arg=paste(strtrim(V(AOP)$name,35),"...",sep=""), col=V(AOP)$po_col, axes=FALSE, cex.names=.75, cex.main=.8, adj=0)
axis(1,at=seq(0,45,5))
text(1, aop.paths.bp, aop.paths(AOP,norm=FALSE), cex=1,pos=4) #adds values to bars
text(.25, aop.paths.bp, V(AOP)$ked ,cex=.5,pos=4 ) #adds OE value text next to bars
dev.off()
