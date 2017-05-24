#### Load proper libraries for AOP Networks and Graph Analyses ####

library("aop")
library(graph)
library(igraph)
library(readr)
library(ggplot2)
library(rgl)
library(plotrix)

workingDir <-"C://Users/NPollesc/Desktop/GitHub/AOPwiki/" #EPA Dir
# workingDir<- "C://Users/Nathan Pollesch/Documents/GitHub/AOPWiki/" #Personal Dir
setwd(workingDir)
#ARMimport<- "data/Aromatase network.txt"
#THYimport<- "data/Thyroid network.txt"

source(paste(workingDir,"R_files/AOP_net_functions.R",sep="")) #imports custom functions
#KERdata<-read.table(paste(workingDir, KERimport, sep=""), sep="\t", stringsAsFactors=FALSE, header=TRUE)
#KEdata<-read.table(paste(workingDir, KEimport, sep=""), sep="\t", stringsAsFactors=FALSE, header=TRUE)
#KEPdata<-read.table(paste(workingDir, KEplus, sep=""), sep="\t", stringsAsFactors=FALSE, header=TRUE)

#### AOP: Import and format network ####
#Import network as cytoscape object and turn in aop object
#stea_aop <- convert_cytoscape_to_aop("//aa.ad.epa.gov/ORD/DUL/USERS/NPollesc/Net MyDocuments/projects/Pellston_2017/Ed_Anze/Programs/steatosis/Steatosis_AOP/steatosis.cyjs")
AOP_net<-convert_cytoscape_to_aop(paste(workingDir,"data/thyroid_network.cyjs",sep=""))

#Import additional data exported from cytoscape, including node names and event types
#stea_info <- read_csv("//aa.ad.epa.gov/ORD/DUL/USERS/NPollesc/Net MyDocuments/projects/Pellston_2017/Ed_Anze/Programs/steatosis/Steatosis_AOP/stea_nodes.csv")
AOP_data<- read_csv(paste(workingDir,"data/thyroid_nodes.csv",sep=""))
#convert aop object to graphNEL object
AOP_graph <- convert_aop_to_graph(AOP_net)
#convert graphNEL object to igraph object for analysis using igraph package
AOP<-igraph.from.graphNEL(AOP_graph, name = TRUE, weight = TRUE,
                         unlist.attrs = TRUE)

## ASSIGN NAMES, ATTRIBUTES, and NODE COLORS BY KED
# Name the AOP network
AOP.name<-"Thyroid"

#V(AOP)names were imported as 'SUID', not 'name', this tells how to match the entries in table_node to the entries in sg graph object by SUID
match(V(AOP)$name,as.character(AOP_data$SUID))
# assign names to igraph object (note that somehow the node names in V(sg)$name are in descending order, so 'rev' is needed to map actual names properly)
V(AOP)$name<-rev(AOP_data$name)

# assign key event descriptor types (MIE,AO,KE,etc...) And color MIE and AO
V(AOP)$ked<-rev(AOP_data$'Event Type')
E(AOP)$color<-"gray"
V(AOP)$color<-"Yellow"
V(AOP)$color[which(V(AOP)$ked=="MIE")]<-"Green"
V(AOP)$color[which(V(AOP)$ked=="AO")]<-"Red"

#### PLOTS NETWORKS ####

# creates plot of igraph object with KED coloring
dev.off()
set.seed(1)
plot(AOP, edge.arrow.size=.25, edge.curved=.1, vertex.size=10, main=paste(AOP.name,"AOP Network"),layout=layout.davidson.harel(AOP), vertex.label.cex=.75, vertex.label=NA )
legend('bottomleft',c("MIE","KE", "AO"), pch=22,
       col="#777777",  pt.bg=c("Green","Yellow","Red"), pt.cex=2, cex=.8, bty="n", ncol=1)

jpeg(filename = paste("AOP_",AOP.name,"_plot.jpeg"),
     width = 800, height = 800, units = "px", pointsize = 12,
     quality = 75)
set.seed(3)
plot(AOP, edge.arrow.size=.25, edge.curved=.1, vertex.size=10, main=paste(AOP.name,"AOP Network"),layout=layout.davidson.harel(AOP) )
legend('bottomleft',c("MIE","KE", "AO"), pch=22,
       col="#777777",  pt.bg=c("Green","Yellow","Red"), pt.cex=2, cex=.8, bty="n", ncol=1)
dev.off()


#### AOP: Network Form Including Cycles, DAG, and condensation graphs ####
# is.dag() reports TRUE for DAG, FALSE if cycles exist
is.dag(AOP)

# Colors cycles for cyclic graphs (Error in function for when graph is acyclic. Need to fix.)
dev.off()
V(AOP)$color<-cycle.color(AOP)[[1]]
E(AOP)$color<-cycle.color(AOP)[[2]]



#Find colors used for connected components, so that they can be added to the legend
which(V(AOP)$color!="Red"&V(AOP)$color!="Yellow"&V(AOP)$color!="Green")
V(AOP)$color[which(V(AOP)$color!="Red"&V(AOP)$color!="Yellow"&V(AOP)$color!="Green")]
V(AOP)$color[11]

set.seed(3)
plot(AOP,edge.arrow.size=.25, edge.curved=.1, vertex.size=10, main=paste(AOP.name,"AOP Network - Connected Components"),layout=layout.davidson.harel(AOP))
legend('bottomleft',c("MIE","KE", "AO","Connected Component"), pch=22,
       col="#777777",  pt.bg=c("Green","Yellow","Red",V(AOP)$color[11]), pt.cex=2, cex=.8, bty="n", ncol=1)

jpeg(filename = paste("AOP_",AOP.name,"_CC_plot.jpeg"),
     width = 800, height = 800, units = "px", pointsize = 12,
     quality = 75)
set.seed(3)
plot(AOP,edge.arrow.size=.25, edge.curved=.1, vertex.size=10, main=paste(AOP.name,"AOP Network - Connected Components"),layout=layout.davidson.harel(AOP))
legend('bottomleft',c("MIE","KE", "AO","Connected Component"), pch=22,
       col="#777777",  pt.bg=c("Green","Yellow","Red",V(AOP)$color[11]), pt.cex=2, cex=.8, bty="n", ncol=1)
dev.off()
#(B.1) if condensation is desired, create map and condense graph

condense.map(AOP)
condense.graph(AOP,condense.map(AOP))

#(B.2) verify that new graph is DAG and plot
is.dag(AOP.cond)

set.seed(7)
plot(AOP.cond, main=paste(AOP.name,"AOP Network - Condensed"), edge.arrow.size=.25, edge.curved=.1, vertex.size=10,layout=layout.davidson.harel(AOP.cond))
legend('bottomleft',c("MIE","KE", "AO","Contracted Node"), pch=22,
       col="#777777",  pt.bg=c("Green","Yellow","Red","Purple"), pt.cex=2, cex=.8, bty="n", ncol=1)


jpeg(filename = paste("AOP_",AOP.name,"_cond_plot.jpeg"),
     width = 800, height = 800, units = "px", pointsize = 12,
     quality = 75)
set.seed(3)
plot(AOP.cond, main=paste(AOP.name,"AOP Network - Condensed"), edge.arrow.size=.25, edge.curved=.1, vertex.size=10,layout=layout.davidson.harel(AOP.cond))
legend('bottomleft',c("MIE","KE", "AO","Contracted Node"), pch=22,
       col="#777777",  pt.bg=c("Green","Yellow","Red","Purple"), pt.cex=2, cex=.8, bty="n", ncol=1)
dev.off()

#### AOP: Topological Sorting of Network ####
#every dag has a topological sorting, so first test if the graph is a DAG
is.dag(AOP.cond)
is.dag(AOP)

#find topological layout
topo.layout(AOP.cond)
topo.layout(AOP)

#plots topologically sorted network using defined layout, plot options are specified for visualization given topological sorting layout
plot(AOP.cond, layout=topo.layout(AOP.cond), edge.arrow.size=.25, edge.curved=1, vertex.size=5, edge.width=1.5, edge.color="Gray", vertex.label.degree=0, vertex.label.dist=.75, main=paste(AOP.name," AOP Network - Topological Sorting"))
legend('bottomright',c("MIE","KE", "AO","Contracted Node"), pch=22,
       col="#777777",  pt.bg=c("Green","Yellow","Red","Purple"), pt.cex=2, cex=.8, bty="n", ncol=1)

jpeg(filename = paste("AOP_",AOP.name,"_topo_plot.jpeg"),
     width = 1000, height = 800, units = "px", pointsize = 12,
     quality = 75)
plot(AOP.cond, layout=topo.layout(AOP.cond), edge.arrow.size=.25, edge.curved=1, vertex.size=5, edge.width=1.5, edge.color="Gray", vertex.label.degree=0, vertex.label.dist=.75, main=paste(AOP.name," AOP Network - Topological Sorting"))
legend('bottomright',c("MIE","KE", "AO","Contracted Node"), pch=22,
       col="#777777",  pt.bg=c("Green","Yellow","Red","Purple"), pt.cex=2, cex=.8, bty="n", ncol=1)
dev.off()


#If AOP not AOP.cond
plot(AOP, layout=topo.layout(AOP), edge.arrow.size=.25, edge.curved=1, vertex.size=5, edge.width=1.5, edge.color="Gray", vertex.label.degree=0, vertex.label.dist=.75, main=paste(AOP.name," AOP Network - Topological Sorting"))
legend('bottomright',c("MIE","KE", "AO"), pch=22,
       col="#777777",  pt.bg=c("Green","Yellow","Red"), pt.cex=2, cex=.8, bty="n", ncol=1)

#save as JPEG
jpeg(filename = paste("AOP_",AOP.name,"_topo_plot.jpeg"),
     width = 1000, height = 800, units = "px", pointsize = 12,
     quality = 75)
plot(AOP, layout=topo.layout(AOP), edge.arrow.size=.25, edge.curved=1, vertex.size=5, edge.width=1.5, edge.color="Gray", vertex.label.degree=0, vertex.label.dist=.75, main=paste(AOP.name," AOP Network - Topological Sorting"))
legend('bottomright',c("MIE","KE", "AO"), pch=22,
       col="#777777",  pt.bg=c("Green","Yellow","Red"), pt.cex=2, cex=.8, bty="n", ncol=1)
dev.off()



##Notes on topological sorting of AOP networks:
  #This sorting can be used to prioritize measurement and have a quick visual respresentation of potential up-stream effects.  I.e. if you are towards the beginning, only few KE can affect you, but towards the end, many more could (not that they necessarily do)
  #This sorting also lines up logically with MIE -> KE -> AO

#### AOP: Path count function and visualizations ####

#!TASK: NEED TO FIX THIS FUNCTION - PERHAPS SOME PATHS DON"T EXIST BETWEEN MIE AND ALL AOS? ####
V(AOP.cond)$aoppn<-aop.paths(AOP.cond, nor=TRUE)
V(AOP.cond)$aopp<-aop.paths(AOP.cond, nor=FALSE)

V(AOP)$aoppn<-aop.paths(AOP, nor=TRUE)
V(AOP)$aopp<-aop.paths(AOP, nor=FALSE)

# Node coloring option based on no of path occurances from MIE to AO
heatgrad=rev(heat.colors(n=max(aop.paths(AOP.cond, norm=FALSE))))
heatgrad<-append(heatgrad, "#FFFFFF", after=0)

heatgrad=rev(heat.colors(n=max(aop.paths(AOP, norm=FALSE))))
heatgrad<-append(heatgrad, "#FFFFFF", after=0) #in the instance that a KE is not included in paths between MIE and AO...

# shows the color scheme with associated values
# pie(rep(1,max(aop.paths(AOP.cond,norm=FALSE))+1),col=rev(heat.colors(n=max(aop.paths(AOP.cond, norm=FALSE))+1)))

#assign color attribute of vertices
V(AOP.cond)$color<-heatgrad[aop.paths(AOP.cond, norm=FALSE)+1]# applies color scheme to nodes


V(AOP)$color<-heatgrad[aop.paths(AOP, norm=FALSE)+1]# applies color scheme to nodes

# ### Topologically sorted plot with path occurance coloration
# #dev.off()
# plot(AOP.cond, layout=topo.layout(AOP.cond), edge.arrow.size=.25, edge.curved=1, vertex.size=6, edge.width=1.5, edge.color="Gray", vertex.label.degree=0, vertex.label.dist=.75, main="Steatosis AOP - Path Occurrence between MIE and AO")
# color.legend(-1.,.5,-1.1,1, grad="y",legend=c(paste("Low (",min(V(AOP.cond)$aopp),")"),paste("High (",max(V(AOP.cond)$aopp),")")), rect.col=heatgrad[c(trunc(seq(0,max(V(AOP.cond)$aopp),length.out=10)),0)])

# AOP Networks with path occurrence colorings 

# for AOP.cond
dev.off()
set.seed(3)
plot(AOP.cond, main=paste(AOP.name,"AOP Network - Path Occurrences"), edge.arrow.size=.5, edge.curved=.1, vertex.size=10,layout=layout.davidson.harel(AOP.cond))
color.legend(1.,.5,1.1,1, grad="y",legend=c(paste("Low (",min(V(AOP.cond)$aopp),")"),paste("High (",max(V(AOP.cond)$aopp),")")), rect.col=heatgrad[c(trunc(seq(0,max(V(AOP.cond)$aopp),length.out=10)),0)])

jpeg(filename = paste("AOP_",AOP.name,"_PO_plot.jpeg"),
     width = 900, height = 800, units = "px", pointsize = 12,
     quality = 75)
set.seed(3)
plot(AOP.cond, main=paste(AOP.name,"AOP Network - Path Occurrences"), edge.arrow.size=.5, edge.curved=.1, vertex.size=10,layout=layout.davidson.harel(AOP.cond))
color.legend(1.,.5,1.1,1, grad="y",legend=c(paste("Low (",min(V(AOP.cond)$aopp),")"),paste("High (",max(V(AOP.cond)$aopp),")")), rect.col=heatgrad[c(trunc(seq(0,max(V(AOP.cond)$aopp),length.out=10)),0)])
dev.off()



dev.off()
set.seed(3)
plot(AOP, main=paste(AOP.name,"AOP Network - Path Occurrences"), edge.arrow.size=.5, edge.curved=.1, vertex.size=10,layout=layout.davidson.harel(AOP))
color.legend(1.,.5,1.1,1, grad="y",legend=c(paste("Low (",min(V(AOP)$aopp),")"),paste("High (",max(V(AOP)$aopp),")")), rect.col=heatgrad[c(trunc(seq(0,max(V(AOP)$aopp),length.out=10)),0)])

jpeg(filename = paste("AOP_",AOP.name,"_PO_plot.jpeg"),
     width = 900, height = 800, units = "px", pointsize = 12,
     quality = 75)
set.seed(3)
plot(AOP, main=paste(AOP.name,"AOP Network - Path Occurrences"), edge.arrow.size=.5, edge.curved=.1, vertex.size=10,layout=layout.davidson.harel(AOP))
color.legend(1.,.5,1.1,1, grad="y",legend=c(paste("Low (",min(V(AOP)$aopp),")"),paste("High (",max(V(AOP)$aopp),")")), rect.col=heatgrad[c(trunc(seq(0,max(V(AOP)$aopp),length.out=10)),0)])
dev.off()

### BARPLOT OF aop.paths() data
par(mfrow=c(1,1),las=2, mar=c(5,15,4,3)) # make label text perpendicular to axis and increase margin
#V(sg)$color<-heatgrad[node_counts$count+1]# applies color scheme to nodes 
aop.paths.bp<-barplot(aop.paths(AOP.cond), xlim=c(0,max(aop.paths(AOP.cond))+1), main="Steatosis AOP - Number of occurrences in paths between MIEs and AOs", horiz=TRUE, names.arg=V(AOP.cond)$name, col=V(AOP.cond)$color, axes=FALSE, cex.names=.75, cex.main=.8, adj=0)

#bpall<-barplot(node_counts$count, xlim=c(0,max(counts)+1), main="Number of occurances in paths between MIEs and AOs", horiz=TRUE, names.arg=names(counts), col=V(sg)$color, axes=FALSE, cex.names=.75, cex.main=.8, adj=0, beside=FALSE)
axis(1,at=seq(0,45,5))
text(1, aop.paths.bp, aop.paths(AOP.cond,norm=FALSE), cex=1,pos=4) #adds values to bars
text(.25, aop.paths.bp, V(AOP.cond)$ked ,cex=.5,pos=4 ) #adds OE value text next to bars

## BARPLOT FOR AOP not AOP.cond
dev.off()
par(mfrow=c(1,1),las=2, mar=c(5,15,4,3)) # make label text perpendicular to axis and increase margin
#V(sg)$color<-heatgrad[node_counts$count+1]# applies color scheme to nodes 
aop.paths.bp<-barplot(aop.paths(AOP), xlim=c(0,max(aop.paths(AOP))+1), main="Epilepsy AOP - Number of occurrences in paths between MIEs and AOs", horiz=TRUE, names.arg=V(AOP)$name, col=V(AOP)$color, axes=FALSE, cex.names=.75, cex.main=.8, adj=0)

#bpall<-barplot(node_counts$count, xlim=c(0,max(counts)+1), main="Number of occurances in paths between MIEs and AOs", horiz=TRUE, names.arg=names(counts), col=V(sg)$color, axes=FALSE, cex.names=.75, cex.main=.8, adj=0, beside=FALSE)
axis(1,at=seq(0,45,5))
text(1, aop.paths.bp, aop.paths(AOP,norm=FALSE), cex=1,pos=4) #adds values to bars
text(.25, aop.paths.bp, V(AOP)$ked ,cex=.5,pos=4 ) #adds OE value text next to bars

jpeg(filename = paste("AOP_",AOP.name,"_PO_barplot.jpeg"),
     width = 900, height = 800, units = "px", pointsize = 12,
     quality = 75)
par(mfrow=c(1,1),las=2, mar=c(5,15,4,3))
aop.paths.bp<-barplot(aop.paths(AOP), xlim=c(0,max(aop.paths(AOP))+1), main=paste(AOP.name,"AOP - Number of occurrences in paths between MIEs and AOs"), horiz=TRUE, names.arg=V(AOP)$name, col=V(AOP)$color, axes=FALSE, cex.names=.75, cex.main=.8, adj=0)
axis(1,at=seq(0,45,5))
text(1, aop.paths.bp, aop.paths(AOP,norm=FALSE), cex=1,pos=4) #adds values to bars
text(.25, aop.paths.bp, V(AOP)$ked ,cex=.5,pos=4 ) #adds OE value text next to bars
dev.off()

#~ Misc, notes, and tasks ####

#Create random DAGS
library(pcalg)
set.seed(101)
myDAG <- randomDAG(n = 12, prob= 0.2, lB = 0.1, uB = 1)

myDAGi<-igraph.from.graphNEL(myDAG)
dev.off()
plot(myDAGi)

#Create random graphs

g <- erdos.renyi.game(30, 2/30, dir=TRUE)
plot(g)

components(g, mode="strong")

V(g)$color<-cycle.color(g)[[1]]

cycle.color(g)

E(g)$color<-cycle.color(g)[[2]]

# betweenness values provide DESCRIBE
V(sg.cond)$btw<-betweenness(sg.cond,normalized=TRUE)

##create adjacency matrix 
A_stea<-get.adjacency(sg)
AM_stea<-as.matrix(A_stea)

####~ Example: DAG, Cycles, and Condensation Maps of Random Graphs ####
#(A.1) generate random igraph object, call it 'g'
g <- erdos.renyi.game(50, 2/50, dir=TRUE)

#(A.2) test for cycles and plot
is.dag(g)
plot(g)

#(A.3) if cycles exist, highlight them using cycle.color() and plot
V(g)$color<-cycle.color(g)[[1]]
E(g)$color<-cycle.color(g)[[2]]
plot(g)

#(B.1) if condensation is desired, create map and condense graph
condense.map(g)
condense.graph(g,condense.map(g))

#(B.2) verify that new graph is DAG and plot
is.dag(g.cond)
plot(g.cond)
