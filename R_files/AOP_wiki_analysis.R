#### Visualize Adverse Outcome Pathway (AOP) WIKI Data ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##  Set working directory and import key event relationships
library(igraph)
workingDir <-"C://Users/NPollesc/Desktop/GitHub/AOPwiki/" #EPA Dir
# workingDir<- "C://Users/Nathan Pollesch/Documents/GitHub/AOPWiki/" #Personal Dir
setwd(workingDir)
KERimport <- "data/all-KERs.txt"
KEimport <- "data/all-KEs.txt"
KEplus <- "data/all_KEs_plus.txt"

source(paste(workingDir,"R_files/AOP_net_functions.R",sep="")) #imports custom functions
KERdata<-read.table(paste(workingDir, KERimport, sep=""), sep="\t", stringsAsFactors=FALSE, header=TRUE)
KEdata<-read.table(paste(workingDir, KEimport, sep=""), sep="\t", stringsAsFactors=FALSE, header=TRUE)
KEPdata<-read.table(paste(workingDir, KEplus, sep=""), sep="\t", stringsAsFactors=FALSE, header=TRUE)

#Format KEPdata for easier handling later
KEPdata[,1]<-as.numeric(substring(KEPdata[,1],5)) #strips the characters Aop: from AOPID column and turns the result numeric
KEPdata[,3]<-as.numeric(substring(KEPdata[,3],7)) #strips the characters Event: from the Event ID column and turns the result numeric

##  Identify all unique KEs
allKEs<-c(KERdata[,2],KERdata[,4])
uniqueKEs<-unique(allKEs)
keID<-data.frame(ID=1:length(uniqueKEs),KE=uniqueKEs)

KERs<-cbind(match(KERdata[,2], keID[,2]),match(KERdata[,4], keID[,2]))
AOPg<-graph_from_edgelist(KERs, directed=T)

## Add names for key event nodes
V(AOPg)$KE_name<-as.character(keID$KE)
V(AOPg)$name<-keID$KE
V(AOPg)$KE_EID<-KEdata[match(V(AOPg)$KE_name,KEdata[,2]),1] # adds event ID number
V(AOPg)$AOP_ID<-KEPdata[match(V(AOPg)$KE_EID,KEPdata[,3]),1] # finds AOPID to add to V(AOPg) data

# This identifies which KEs are included in KERs, but are not themselves included in the KE event listings.
# V(AOPg)$KE_name[which(is.na(V(AOPg)$KE_EID))]

# Plot the AOP wiki, colored by AOP
length(unique(V(AOPg)$AOP_ID))
length(E(AOPg))
length(V(AOPg))

acols=topo.colors(length(unique(V(AOPg)$AOP_ID)))
for(i in 1:length(unique(V(AOPg)$AOP_ID))){
  V(AOPg)[which(V(AOPg)$AOP_ID==unique(V(AOPg)$AOP_ID)[i])]$acol<-acols[i]
  }
par(bg="black",xpd=FALSE)
set.seed(1)
plot(AOPg,vertex.color=V(AOPg)$acol,vertex.label=NA, vertex.size=2, edge.arrow.size=.1)

AOP_freqs<-table(V(AOPg)$AOP_ID)
bp_wcc<-barplot(table(V(AOPg)$AOP_ID),col.axis="white", xlab="AOP ID",ylab="# Key Events",col.lab="white")
abline(h=mean(AOP_freqs),col="red")



#  TASK: WORK ON EDGE COLORING FOR AOP ID
#  edgecombcc<-expand.grid(V(gr)[which(comps$membership==ntcomps[i])],V(gr)[which(comps$membership==ntcomps[i])]) #creates a pairwise list of all nodes in the cc
#  edgecombflat<-as.vector(rbind(edgecombcc[[1]],edgecombcc[[2]])) #flattens the pairwise list to a vector where entries are read pairwise
#  edges.in.cc<-get.edge.ids(gr,edgecombflat,directed=TRUE)
#  E(gr)$color[edges.in.cc]<-cols[[i]]

  
#### CONNECTED COMPONENTS ANALYSIS ####

## Color vertices and edges by their weakly or strongly connected components.
V(AOPg)$cc_color<-unlist(color.comps(AOPg,"strong")$vcol)  #color.comps is a custom function stored in the AOP_net_functions.R file
E(AOPg)$cc_color<-unlist(color.comps(AOPg,"strong")$ecol)  #color.comps is a custom function stored in the AOP_net_functions.R file

# When the "strong" option is passed to color.comps, vsize and ewidth are calculated and can be used within plot
V(AOPg)$cc_size<-unlist(color.comps(AOPg,"strong")$vsize)
E(AOPg)$cc_width<-unlist(color.comps(AOPg,"strong")$ewidth)

# plot of connected components
par(bg="black")
set.seed(1)
plot(AOPg, vertex.size=V(AOPg)$cc_size, edge.width=E(AOPg)$cc_width, vertex.color=V(AOPg)$cc_color, edge.color=E(AOPg)$cc_color, edge.arrow.size=.1, vertex.label=NA)

## barplot for size of weakly connected components
wcomps<-components(AOPg,mode="weak")
wcc_freqs<-table(wcomps$csize)
bp_wcc<-barplot(table(wcomps$csize),col.axis="white", xlab="Component size",ylab="Frequency",col.lab="white")

# This points out how many of the feedback loops/cycles are contained within the same AOP and how many are a result of the network
scomps<-components(AOPg,mode="strong")
ntcomps<-which(scomps$csize>1) # non-trivial ccs (i.e. with more than 1 node)
V(AOPg)$scc<-scomps$membership # assign the attribute scc to nodes based on their membership 
for(i in 1:length(ntcomps)){
  print(V(AOPg)[which(V(AOPg)$scc==ntcomps[i])]$AOP_ID)
  }



#### Level of biological organization plot ####

## Add level of biological organization for key event nodes
V(AOPg)$lobo<-KEdata[[4]][match(V(AOPg)$KE_name,KEdata[[2]])]
V(AOPg)$lobo[which(is.na(V(AOPg)$lobo))]<-"" #assigns blank to NA data
tcols=rainbow(length(unique(V(AOPg)$lobo))) #creates a color scheme for visualization
lobo_order=c("Molecular","Cellular","Tissue","Organ","Individual","Population","") #creates an ordering of biological organization
V(AOPg)$lobo_o<-match(V(AOPg)$lobo,lobo_order) #assigns a value of biological organization instead of string.  1=molecular, 2=cellular, ...
lobo_freqs<-table(V(AOPg)$lobo_o)

## Barplot of lobo frequency
par(bg="black")
xx<- barplot(table(V(AOPg)$lobo_o), col=tcols, axes=F,names.arg=NA)
text(x=xx, y=10, label=lobo_freqs, cex=.75)
legend('topright',c("Molecular","Cellular","Tissue","Organ","Individual","Population","Not Specified"), pch=22,
       col="#777777", xjust=1,yjust=1, pt.bg=tcols, pt.cex=2, cex=.8, bty="n", ncol=1, y.intersp=.5, box.col="black", text.col="white")

#plot the AOP wiki by lobo info
V(AOPg)$lobo_col<-tcols[V(AOPg)$lobo_o]
set.seed(1)
plot(AOPg, vertex.size=2, edge.color="gray", edge.arrow.size=.1, vertex.label=NA, vertex.color=V(AOPg)$lobo_col)


# #### Graph condensation ####
# 
# condense.map(AOPg)
# condense.graph(AOPg,condense.map(AOPg))
# plot(AOPg.cond,vertex.size=1, edge.arrow.size=.1, vertex.label=NA)
# 
# is.dag(AOPg)
# is.dag(AOPg.cond)
# 
# #Attempt a topological sorting of the entire AOPwiki (Hint, it is unreadable :) )
# plot(AOPg.cond,vertex.size=2, edge.arrow.size=.1,layout=topo.layout(AOPg.cond))

#### CENTRALITY MEASURES FOR THE AOPWIKI ####

#### Degree centrality ####

V(AOPg)$cent_size<-1
V(AOPg)$cent_col<-"white"

# Which key event has the most incident nodes?
sort(degree(AOPg, mode="in"))
  V(AOPg)$KE_name[which(V(AOPg)$name==449)]
  V(AOPg)$cent_size[which(V(AOPg)$name==345)]<-3
  V(AOPg)$cent_col[which(V(AOPg)$name==345)]<-"blue"
  
#global degree coloring for network plot
V(AOPg)$deg_col<-deg.col(AOPg)
#colored by degree
plot(AOPg, vertex.size=2, vertex.color=V(AOPg)$deg_col, edge.arrow.size=.1, vertex.label=NA, edge.color="gray",edge.width=1)
#colored by degree and sized by degree
plot(AOPg, vertex.size=500*degree(AOPg,mode="all",normalized=TRUE), vertex.color=V(AOPg)$deg_col, edge.arrow.size=.1, vertex.label=NA, edge.color="gray",edge.width=1)
#barplot to show degree histogram
barplot(table(degree(AOPg)),col=rev(heat.colors(max(degree(AOPg,mode="all")))))

#### Betweenness Centrality ####

# Which key event is involved in the most shortest paths between other key events?
sort(betweenness(AOPg))
which(V(AOPg)$name==345)
V(AOPg)$cent_size[which(V(AOPg)$name==345)]<-3
V(AOPg)$cent_col[which(V(AOPg)$name==345)]<-"green"
V(AOPg)$deg_col<-deg.col(AOPg)
#colored by betweenness

wbpal=colorRampPalette(c("white","blue"))
V(AOPg)$bet_col<-wbpal(10)[as.numeric(cut(betweenness(AOPg),breaks = 10))]

par(bg="black")
plot(AOPg, vertex.size=2, vertex.color=V(AOPg)$bet_col, edge.arrow.size=.1, vertex.label=NA, edge.color="gray",edge.width=1)
#colored by degree and sized by degree
plot(AOPg, vertex.size=1000*betweenness(AOPg,normalized=TRUE), vertex.color=V(AOPg)$bet_col, edge.arrow.size=.1, vertex.label=NA, edge.color="gray",edge.width=1)
#barplot to show degree histogram
barplot(table(degree(AOPg)),col=rev(heat.colors(max(degree(AOPg,mode="all")))))




# Which key event is the closest to the rest?
sort(closeness(AOPg))
V(AOPg)$KE_name[which(V(AOPg)$name==711)]
V(AOPg)$cent_size[which(V(AOPg)$name==711)]<-3
V(AOPg)$cent_col[which(V(AOPg)$name==711)]<-"purple"

set.seed(1)
par(bg="black")
plot(AOPg, vertex.size=V(AOPg)$cent_size, vertex.color=V(AOPg)$cent_col, edge.arrow.size=.1, vertex.label=NA, edge.color="white")#, vertex.color="orange",edge.color="gray")

#### Thyroid AOP network ####
THYimport <- "data/thyroid_net_direct.txt"
THYdata<-read.table(paste(workingDir, THYimport, sep=""), sep="\t", stringsAsFactors=FALSE, header=TRUE)

THYKERdata<-match(THYdata[[2]],KERdata[[1]])
TKERdata<-KERdata[THYKERdata,]
allTKEs<-c(TKERdata[,2],TKERdata[,4])
uniqueTKEs<-unique(allTKEs)
tkeID<-data.frame(ID=1:length(uniqueTKEs),TKE=uniqueTKEs)
TKERs<-cbind(match(TKERdata[,2],tkeID[,2]),match(TKERdata[,4],tkeID[,2]))

TAOPg<-graph_from_edgelist(TKERs,directed=T)
dev.off()
plot(TAOPg, vertex.size=2,vertex.label=NA,edge.arrow.size=.5)

