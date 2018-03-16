#libraries
library(igraph)
library(prodlim)
library(RColorBrewer)

#Directories
workingDir<-"C:\\Users\\obrienja\\Documents\\GitHub\\AOPWiki\\R_files_Jason\\"
dataDir<-"C:\\Users\\obrienja\\Documents\\GitHub\\AOPWiki\\data\\"

#load source of functions
source(paste(workingDir,"AOP-Net-functions-JOB.R",sep="")) #imports JASON custom functions

# Identifies location of data files
KERimport <- "all-KERs.txt"
KEimport <- "all-KEs.txt"
KEked <- "all-KEs-KED.txt" # Additional ontology information file (ex: MIE, KE, AO)
KEplus <- "all-KEs-plus.txt" # Additional ontology information file
KERwoe<- "all-ke-ker-woe.txt" # WOE for KER

# import data
KERdata<-read.table(paste(dataDir, KERimport, sep=""), sep="\t", stringsAsFactors=FALSE, header=TRUE)
KEdata<-read.table(paste(dataDir, KEimport, sep=""), sep="\t", stringsAsFactors=FALSE, header=TRUE)
KEKdata<-read.table(paste(dataDir, KEked, sep=""), sep="\t", stringsAsFactors=FALSE, header=F)
KEPdata<-read.table(paste(dataDir, KEplus, sep=""), sep="\t", stringsAsFactors=FALSE, header=TRUE)
KERWdata<-read.table(paste(dataDir, KERwoe, sep=""), sep="\t", stringsAsFactors=FALSE, header=T)

# Format data for easier handling later
KEPdata[,1]<-as.numeric(substring(KEPdata[,1],5)) #strips the characters Aop: from AOPID column and turns the result numeric
KEPdata[,3]<-as.numeric(substring(KEPdata[,3],7)) #strips the characters Event: from the Event ID column and turns the result numeric
KEKdata[,2]<-as.numeric(substring(KEKdata[,2],7)) #strips the characters Event: from the Event ID column and turns the result numeric
KERWdata[,2]<-as.numeric(substring(KERWdata[,2],7)) #strips the characters Event: from the Event ID column and turns the result numeric
KERWdata[,3]<-as.numeric(substring(KERWdata[,3],7)) #strips the characters Event: from the Event ID column and turns the result numeric
KERWdata[,4]<-as.numeric(substring(KERWdata[,4],14))



### Create AOP wiki igraph object (AOPg)

# Identify all unique KEs
allKEs<-c(KERdata[,2],KERdata[,4])
uniqueKEs<-unique(allKEs)
keID<-data.frame(ID=1:length(uniqueKEs),KE=uniqueKEs)
KERs<-cbind(match(KERdata[,2], keID[,2]),match(KERdata[,4], keID[,2]))

# Create AOP igraph object
AOPg<-graph_from_edgelist(KERs, directed=T)

# Add names for key event nodes
V(AOPg)$KE_name<-as.character(keID$KE)

# Add event IDS
V(AOPg)$KE_EID<-KEdata[match(V(AOPg)$KE_name,KEdata[,3]),1] # adds event ID number
V(AOPg)$name<-V(AOPg)$KE_EID #changes default 'name' object to be the KE event ID number
V(AOPg)$name<-KEdata[match(V(AOPg)$KE_name,KEdata[,3]),1]

# Add AOP IDs. Note: Each KE may belong to more than one AOP, so the AOP_IDs object is a list not just a single AOP ID
for(i in 1:length(V(AOPg))){
  if(length(which(!is.na(match(KEPdata[,3],V(AOPg)$KE_EID[i]))))>0){
    V(AOPg)[i]$AOP_ID<-list(unique(KEPdata[which(!is.na(match(KEPdata[,3],V(AOPg)[i]$KE_EID))),1]))}
  else{V(AOPg)[i]$AOP_ID<-NA}
}

# Add key event designation data ("MIE, AO, or KE")
V(AOPg)$KE_KED<-KEKdata[match(V(AOPg)$KE_EID,KEKdata[,2]),4] # finds KED (Key Event Designator) to add to V(AOPg) data
length(V(AOPg)$KE_KED[which(is.na(V(AOPg)$KE_KED))]) #The number of KEs without KEDs
V(AOPg)$KE_KED[which(is.na(V(AOPg)$KE_KED))]<-"KE" # ALL KEs without KED (NA values from file) are assigned as generic KE

# Add Path Designator (origin or terminus)
AOPg<-add_KE_PD(AOPg)

# add colour attribute (MIE=green, KE=white, AOP=red, non-MIE origin = dark green, non-AO terminus = brown)
V(AOPg)$col<-"white"
V(AOPg)$col[V(AOPg)$KE_KED=="MIE"]<-"green"
V(AOPg)$col[V(AOPg)$KE_KED=="AO"]<-"red"
#V(AOPg)$col[V(AOPg)$KE_PD=="origin"]<-"khaki"  #if KE_PD is used
#V(AOPg)$col[V(AOPg)$KE_PD=="terminus"]<-"pink" #if KE_PD is used 


# map edge attributes
# keID to KE_EID to KE_name map
keID_map<-data.frame(keID=as.integer(V(AOPg)),KE_EID=V(AOPg)$KE_EID, name=V(AOPg)$KE_name)
ker_by_eid<-data.frame(
  ID=KERdata$ID,
  from=as.character(keID_map[match(KERdata$Upstream.event, keID_map$name),"KE_EID"]),
  to=as.character(keID_map[match(KERdata$Downstream.event, keID_map$name),"KE_EID"]),
  stringsAsFactors=FALSE
)

# add KER ID as edge attribute
E(AOPg, P=as.vector(t(ker_by_eid[,2:3])))$KER_ID<-ker_by_eid$ID  

# edge color
E(AOPg)$col<-"grey"

#add direct vs indirect as edge attribute (assume all unspeficied are "directly leads to")
E(AOPg)$DvsI<-"directly leads to"       #assign all as "directly leads to" to cover undefined edges
E(AOPg)$DvsI[E(AOPg)$KER_ID%in%KERWdata$KER_ID] <- KERWdata$DIRECT_INDIRECT[match(E(AOPg)$KER_ID[E(AOPg)$KER_ID%in%KERWdata$KER_ID],KERWdata$KER_ID)]

# of direct (adjacent) and indirect (non-adjacent) edges, as defined by users
sum(E(AOPg)$DvsI=="directly leads to")
# 893 direct (adjacent) KERs
sum(E(AOPg)$DvsI=="indirectly leads to")
# 165 indirect (non-adjacent) KERs



### PLOT AOPg

# generate plot coordinates
set.seed(1)
layout.AOP<-layout.fruchterman.reingold(AOPg, weight=rep(0.4,ecount(AOPg)))

# map plot coordinates to vertices (so that "subgraphs" can easily have the same plot coords)
V(AOPg)$plotX<-layout.AOP[,1]
V(AOPg)$plotY<-layout.AOP[,2]

# plot
plotLay<-cbind(V(AOPg)$plotX,V(AOPg)$plotY) #defining before plot to make it easier to copy and paste for other "subgraphs"
par(mar=c(0,0,0,0))
plot(AOPg, vertex.size=2.5, vertex.color=V(AOPg)$col, edge.width=1, edge.arrow.size=0.2, vertex.label=NA, layout=plotLay)
#with V names
#plot(AOPg, vertex.size=2.5, vertex.color=V(AOPg)$col, edge.width=1, edge.arrow.size=0.2, vertex.label=V(AOPg)$name,vertex.label.cex=0.2, layout=plotLay)



### identify linear AOPs

# Linear AOPs
laops_AOPg<-linear.AOPs(AOPg, use_KE_PD=FALSE)
sum(sapply(laops_AOPg,length))
# 43252 "linear AOPs" ( 45380 when using "origin" and "terminus" instead of MIE to AO (ie use_KE_PD=TRUE))

# MIE pair with largest number of simple paths
m<-max(sapply(laops_AOPg,length))
# 3426
names(laops_AOPg)[sapply(laops_AOPg,length)==m]
# [1] "424 563"   MIE=424 AO=563 has 3426 linear AOPs between them



### Component analysis

# weak components
AOPg_wComps<-components(AOPg, mode="weak")
# 34 weak components

# strong components
AOPg_sComps<-components(AOPg, mode="strong")
ntcomps<-which(AOPg_sComps$csize>1) # non-trivial ccs (i.e. with more than 1 node)
length(ntcomps)
#7 strong comps

# assign the attribute scc to nodes based on their membership 
V(AOPg)$scc<-AOPg_sComps$membership 

# color and width attributes of edges in scc
E(AOPg)$scc_col<-E(AOPg)$col 
E(AOPg)$scc_width<-1
sccPal<-brewer.pal(n=length(ntcomps), name="Dark2")
for(i in 1:length(ntcomps)){
  subG<-induced_subgraph(AOPg, V(AOPg)[V(AOPg)$scc==ntcomps[i]])
  E<-as.character(as.vector(t(as_edgelist(subG))))
  E(AOPg, P=E)$scc_col<-sccPal[i]
  E(AOPg, P=E)$scc_width<-2
}

#plot
plotLay<-cbind(V(AOPg)$plotX,V(AOPg)$plotY) 
par(mar=c(0,0,0,0))
plot(AOPg, vertex.size=2.5, vertex.color=V(AOPg)$col, vertex.label=NA, edge.width=E(AOPg)$scc_width, edge.arrow.size=0.2, edge.color=E(AOPg)$scc_col, layout=plotLay)



#contract strong components -AOPg
AOPg.con<-contract.scc(AOPg)
#plot
plotLay<-cbind(V(AOPg.con)$plotX,V(AOPg.con)$plotY) 
par(mar=c(0,0,0,0))
plot(AOPg.con, vertex.size=V(AOPg.con)$size, vertex.color=V(AOPg.con)$col, vertex.label=NA, edge.width=E(AOPg.con)$width, edge.arrow.size=0.2, edge.curved=0, edge.color=E(AOPg.con)$col, layout=plotLay)

# linear AOPs in contracted network (AOPg.con)
laops_AOPg.con<-linear.AOPs(AOPg.con)
sum(sapply(laops_AOPg.con,length))
# 20397 linear AOPs in contracted network (almost halved!)

# MIE pair with largest number of simple paths
m<-max(sapply(laops_AOPg.con,length))
m
# 1887
names(laops_AOPg.con)[sapply(laops_AOPg.con,length)==m]
# [1] "424 563"   MIE=424 AO=563 has 1887 linear AOPs between them when SCCs are contracted (almost halved!)

