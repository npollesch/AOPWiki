#libraries
library(igraph)
library(prodlim)
library(RColorBrewer)

#Directories
workingDir<-"C:\\Users\\obrienja\\Documents\\GitHub\\AOPWiki\\R_files_Jason\\"

#load source of functions
source(paste(workingDir,"AOP-Net-functions-JOB.R",sep="")) #imports JASON custom functions


### This script uses several objects that are created from AOP-Net-AOPWiki.R script:
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###     "AOPg" igraph object
###     "laops_AOPg" pathList of linear AOPs (can be created based on MIE/AO pairs or origin/terminus pairs)
###     "KERWdata" data frame that contains KER info


### User Defined Adj vs NonAdj

#instances where KER is ALWAYS defined by users as direct (adjacent)
numE<-unique(KERWdata$KER_ID)
allDir<-vector()
for(i in 1:length(numE)){
  E<-KERWdata$DIRECT_INDIRECT[KERWdata$KER_ID==numE[i]]
  allDir<-c(allDir,all(E=="directly leads to"))
}
sum(allDir)
# 780

#instances where KER is ALWAYS defined by users indirect (non adjacent)
numE<-unique(KERWdata$KER_ID)
allDir<-vector()
for(i in 1:length(numE)){
  E<-KERWdata$DIRECT_INDIRECT[KERWdata$KER_ID==numE[i]]
  allDir<-c(allDir,all(E=="indirectly leads to"))
}
sum(allDir)
# 163

# Instances where the same KER is "direct" in one AOP but "indirect" in another, as defined by users
numE<-unique(KERWdata$KER_ID)
allDir<-vector()
for(i in 1:length(numE)){
  E<-KERWdata$DIRECT_INDIRECT[KERWdata$KER_ID==numE[i]]
  allDir<-c(allDir,all(E==E[1]))
}
sum(!allDir)
# 12 instances where the same KER is "direct" in one AOP but "indirect" in another



### Algorithm Based Identification of Adj vs NonAdj

# extract list of paths that contain only adjacent KERs (adjPaths)
adjPaths_AOPg<-remove.nonAdj(laops_AOPg)

# combine adjacent paths into an igraph object
AOPg_adj<-graph_from_pathList(adjPaths_AOPg)

#Map vertex attributes from AOpg to AOPg_adj
mapCoords<-match(V(AOPg_adj)$name,V(AOPg)$name)
V(AOPg_adj)$KE_KED<-V(AOPg)$KE_KED[mapCoords]
V(AOPg_adj)$KE_PD<-V(AOPg)$KE_PD[mapCoords]  #if used
V(AOPg_adj)$col<-V(AOPg)$col[mapCoords]
V(AOPg_adj)$plotX<-V(AOPg)$plotX[mapCoords]
V(AOPg_adj)$plotY<-V(AOPg)$plotY[mapCoords]
V(AOPg_adj)$AOP_ID<-V(AOPg)$AOP_ID[mapCoords]
E(AOPg_adj)$col<-"grey"

# compare edges to original graph to identify the "removed" ones (ie the non adjacent KERs)
non.adj.KERs<-edge_difference(AOPg_adj, AOPg)
nrow(non.adj.KERs)
# 156 KERs identified as "non-adjacent"

# summary of KERs
KER_Summary<-as.data.frame(as_edgelist(AOPg))
colnames(KER_Summary)<-c("KEup","KEdown")
KER_Summary$adj<-"adjacent"
KER_Summary$adj[row.match(non.adj.KERs,KER_Summary[,1:2])]<-"NON ADJACENT"
sum(KER_Summary$adj=="adjacent")
# 762 adjacent edges between MIE and AO paths (902 when based on origin and terminus paths)
sum(KER_Summary$adj=="NON ADJACENT")
# 296 adjacent edges between MIE and AO paths (156 when based on origin and terminus paths)

# add adjacency as an edge attribute to original AOPg
eCoords<-as.character(as.vector(t(KER_Summary[,1:2])))
E(AOPg, P=eCoords)$adjacency<-KER_Summary$adj
E(AOPg)$col<-"grey"
E(AOPg)$col[E(AOPg)$adjacency=="NON ADJACENT"]<-"orange"
E(AOPg)$width<-1
E(AOPg)$width[E(AOPg)$adjacency=="NON ADJACENT"]<-1.5


# Plot
# all edges (non adjacent are thick and orange)
plotLay<-cbind(V(AOPg)$plotX,V(AOPg)$plotY) 
par(mar=c(0,0,0,0))
plot(AOPg, vertex.size=2.5, vertex.color=V(AOPg)$col, vertex.label=NA, edge.width=E(AOPg)$width, edge.arrow.size=0.2,edge.curved=0, edge.color=E(AOPg)$col, layout=plotLay)

# only adjacent edges
plotLay<-cbind(V(AOPg_adj)$plotX,V(AOPg_adj)$plotY)
par(mar=c(0,0,0,0))
plot(AOPg_adj, vertex.size=2.5, vertex.color=V(AOPg_adj)$col, edge.width=1, edge.arrow.size=0.2, vertex.label=NA, layout=plotLay)



### linear AOP analysis of adjacent-only graph

#number of linear AOPs
laops_AOPg_adj<-linear.AOPs(AOPg_adj, use_KE_PD=TRUE)
sum(sapply(laops_AOPg_adj,length))
# 4654 linear AOPs when only adjacent MIE to AO paths are considered (5361 for origin/terminus paths)

# MIE pair with largest number of simple paths
m<-max(sapply(laops_AOPg_adj,length))
m
# 1887
names(laops_AOPg_adj)[sapply(laops_AOPg_adj,length)==m]
# MIE   AO   simplePaths
# 201  563           102



### Component Analysis

# weak components
AOPg_adj_wComps<-components(AOPg_adj, mode="weak")
# 26 weak components when adjacent MIE to AO paths are considered

# strong components- AOPg_adj
AOPg_adj_sComps<-components(AOPg_adj, mode="strong")
ntcomps_adj<-which(AOPg_adj_sComps$csize>1) # non-trivial ccs (i.e. with more than 1 node)
length(ntcomps_adj)
#4 strong comps when adjacent MIE to AO paths are considered

# assign attribute scc to nodes based on their membership 
V(AOPg_adj)$scc<-AOPg_adj_sComps$membership 

# color and width attributes of edges in scc
E(AOPg_adj)$scc_col<-E(AOPg_adj)$col 
E(AOPg_adj)$scc_width<-1
sccPal_adj<-brewer.pal(n=length(ntcomps_adj), name="Dark2")
for(i in 1:length(ntcomps_adj)){
  subG<-induced_subgraph(AOPg_adj, V(AOPg_adj)[V(AOPg_adj)$scc==ntcomps_adj[i]])
  E<-as.character(as.vector(t(as_edgelist(subG))))
  E(AOPg_adj, P=E)$scc_col<-sccPal_adj[i]
  E(AOPg_adj, P=E)$scc_width<-2
}

#plot
plotLay<-cbind(V(AOPg_adj)$plotX,V(AOPg_adj)$plotY) 
par(mar=c(0,0,0,0))
plot(AOPg_adj, vertex.size=2.5, vertex.color=V(AOPg_adj)$col, vertex.label=NA, edge.width=E(AOPg_adj)$scc_width, edge.arrow.size=0.2, edge.color=E(AOPg_adj)$scc_col, layout=plotLay)


#contract strong components -AOPg_adj
AOPg_adj.con<-contract.scc(AOPg_adj)
#plot
plotLay<-cbind(V(AOPg_adj.con)$plotX,V(AOPg_adj.con)$plotY) 
par(mar=c(0,0,0,0))
plot(AOPg_adj.con, vertex.size=V(AOPg_adj.con)$size, vertex.color=V(AOPg_adj.con)$col, vertex.label=NA, edge.width=E(AOPg_adj.con)$width, edge.arrow.size=0.2, edge.curved=0, edge.color=E(AOPg_adj.con)$col, layout=plotLay)

# linear AOPs in contracted network (AOPg_adj.con)
laops_AOPg_adj.con<-linear.AOPs(AOPg_adj.con)
sum(sapply(laops_AOPg_adj.con,length))
# 3087 linear AOPs in contracted network (almost halved!)

# MIE pair with largest number of simple paths
m<-max(sapply(laops_AOPg_adj.con,length))
m
# 62
names(laops_AOPg_adj.con)[sapply(laops_AOPg_adj.con,length)==m]
# [1] "201 563"   MIE=424 AO=563 has 62 linear AOPs between them when SCCs are contracted 



