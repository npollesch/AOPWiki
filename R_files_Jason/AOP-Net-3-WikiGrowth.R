#libraries
library(igraph)
library(prodlim)

#Directories
workingDir<-"C:\\Users\\obrienja\\Documents\\GitHub\\AOPWiki\\R_files_Jason\\"

#load source of functions
source(paste(workingDir,"AOP-Net-functions-JOB.R",sep="")) #imports JASON custom functions


### This script uses several objects that are created from AOP-Net-AOPWiki.R script:
### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
###     "AOPg" igraph object

g<-AOPg

### Growth by AOP ID  (NOTE: TAKES ABOUT 2 HOURS TO COMPLETE)

#list all AOP IDs in graph and list in order
IDs<-unlist(V(g)$AOP_ID)
IDs<-unique(IDs)
IDs<-IDs[order(IDs)]

#plce to store metrics
vC<-vector()        # vectex count
eC<-vector()        # edge count
mieC<-vector()      # MIE count
aoC<-vector()       # AO count
oC<-vector()        # origin count
tC<-vector()        # terminus count
adjC<-vector()      # adjacent KER count
nonC<-vector()      # non-adjacent KER count
wcC<-vector()       # weak component count
scC<-vector()       # strong component count
borrowC<-vector()   # instances a new AOP borrows a KE that is already in the network
laopC<-vector()     # number of linear AOPs
laopC_adj<-vector() # number of linear AOPs, only considering adjacent MIE to AO paths

# vector to store list of vertices that are currently in the growing network
vSet<-vector()

# Build and measure the network one AOP_ID at a time
for(i in 1:length(IDs)){

  #location of vertices with AOP_ID=IDs[i]
  v<-sapply(V(g)$AOP_ID, function(x){any(x%in%IDs[i])})  
  
  # number of times the new AOP uses a KE that is already in the network (borrowed a KE)
  b<-sum(V(g)[v]%in%vSet)
  borrowC<-c(borrowC,b)
  
  #add new vertices to list and create subG
  vSet<-unique(c(vSet, V(g)[v]))
  subG<-induced_subgraph(g, vSet)
  
  #number of vertices
  vC<-c(vC,vcount(subG))
  
  #number of edges
  eC<-c(eC,ecount(subG))
  
  #number of MIEs and AOs
  mieC<-c(mieC, sum(V(subG)$KE_KED=="MIE"))
  aoC<-c(aoC, sum(V(subG)$KE_KED=="AO"))
  
  #identify and count origin and terminus
  subG<-add_KE_PD(subG) 
  oC<-c(oC, sum(V(subG)$KE_PD=="origin"))
  tC<-c(tC, sum(V(subG)$KE_PD=="terminus"))
  
  #linear AOPs
  laops<-linear.AOPs(subG, use_KE_PD=FALSE)
  laopC<-c(laopC, sum(sapply(laops,length)))
  
  #identify and count adjacent and non-adjacent KERs
  subG<-add_KER_adjacency(subG)
  adjC<-c(adjC, sum(E(subG)$adjacency=="adjacent"))
  nonC<-c(nonC, sum(E(subG)$adjacency=="non-adjacent"))
  
  # create subgraph of only adjacent edges
  subG_adj<-subgraph.edges(subG, eids=E(subG)[E(subG)$adjacency=="adjacent"])
  
  #adj only linear AOPs
  laops_adj<-linear.AOPs(subG_adj, use_KE_PD=FALSE)
  laopC_adj<-c(laopC_adj, sum(sapply(laops_adj,length)))
  
  #components
  wcC<-c(wcC,components(subG, mode="weak")$no)
  scC<-c(scC,sum(components(subG, mode="strong")$csize>1))
}

wikiGrowth<-data.frame(
  AOP_ID=IDs,
  Vs=vC,
  Es=eC,
  MIEs=mieC,
  AOs=aoC,
  origins=oC,
  termini=tC,
  borrow=borrowC,
  adjE=adjC,
  non_adjE=nonC,
  weakC=wcC,
  strongC=scC,
  LAOPS=laopC,
  ADJ_LAOPS=laopC_adj)

write.table(wikiGrowth, paste(workingDir,"results\\wikiGrowth_by_AOP.txt", sep=""), sep="\t", row.names=FALSE)



