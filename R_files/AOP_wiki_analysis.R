
#### LOAD PACKAGES AND IMPORT AOPWIKI DATA ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

## Needs igraph package
library(igraph)
#localLibDir<-"C:\\Program Files\\R\\R-3.2.0\\library\\"  #Jason's local library directory
#library(igraph, lib.loc = localLibDir)

##  Set working directory
workingDir<-"C://Users/NPollesc/Desktop/GitHub/AOPwiki/" ## Nate's EPA working directory
#workingDir<-"C://Users/Nathan Pollesch/Documents/GitHub/AOPWiki/" ## Nate's personal comp working directory
# workingDir<-"C:\\Users\\obrienja\\Documents\\GitHub\\AOPWiki\\" 
setwd(workingDir)

## Identifies location of data files
KERimport <- "data/all-KERs.txt"
KEimport <- "data/all-KEs.txt"
KEplus <- "data/all-KEs-plus.txt" # Additional ontology information file
KEked <- "data/all-KEs-KED.txt" # Additional ontology information file
KERwoe<- "data/all-ke-ker-woe.txt" # WOE for KER

## source() imports custom functions from associated file
source(paste(workingDir,"R_files/AOP_net_functions.R",sep="")) #imports custom functions
KERdata<-read.table(paste(workingDir, KERimport, sep=""), sep="\t", stringsAsFactors=FALSE, header=TRUE)
KEdata<-read.table(paste(workingDir, KEimport, sep=""), sep="\t", stringsAsFactors=FALSE, header=TRUE)
KEPdata<-read.table(paste(workingDir, KEplus, sep=""), sep="\t", stringsAsFactors=FALSE, header=TRUE)
KEKdata<-read.table(paste(workingDir, KEked, sep=""), sep="\t", stringsAsFactors=FALSE, header=F)

## Format KEPdata for easier handling later
KEPdata[,1]<-as.numeric(substring(KEPdata[,1],5)) #strips the characters Aop: from AOPID column and turns the result numeric
KEPdata[,3]<-as.numeric(substring(KEPdata[,3],7)) #strips the characters Event: from the Event ID column and turns the result numeric
KEKdata[,2]<-as.numeric(substring(KEKdata[,2],7)) #strips the characters Event: from the Event ID column and turns the result numeric

## Identify all unique KEs
allKEs<-c(KERdata[,2],KERdata[,4])
#fullKEs<-c(KERdata[,2],KERdata[,4],KERdata[,1])
#length(allKEs)
uniqueKEs<-unique(allKEs)
#uniqueFullKEs<-fullKEs[!duplicated(fullKEs, incomparables = fullKEs[,3])]
#length(uniqueKEs)
#length(uniqueFullKEs)
keID<-data.frame(ID=1:length(uniqueKEs),KE=uniqueKEs)
KERs<-cbind(match(KERdata[,2], keID[,2]),match(KERdata[,4], keID[,2]))
AOPg<-graph_from_edgelist(KERs, directed=T)

## Add names for key event nodes
V(AOPg)$KE_name<-as.character(keID$KE)
## Add event and AOP IDs
V(AOPg)$KE_EID<-KEdata[match(V(AOPg)$KE_name,KEdata[,2]),1] # adds event ID number
V(AOPg)$AOP_ID<-KEPdata[match(V(AOPg)$KE_EID,KEPdata[,3]),1] # finds AOPID to add to V(AOPg) data
V(AOPg)$name<-V(AOPg)$KE_EID

## Add key event designation data
V(AOPg)$KE_KED<-KEKdata[match(V(AOPg)$KE_EID,KEKdata[,2]),4] # finds KED (Key Event Designator) to add to V(AOPg) data
V(AOPg)$KE_KED[which(is.na(V(AOPg)$KE_KED))]<-"KE" # ALL KEs without KED (NA values from file) are assigned as generic KE

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##
## Builds AOPwiki network from Edgelist with weights ##
KERWdata<-read.table(paste(workingDir, KERwoe, sep=""), sep="\t", stringsAsFactors=FALSE, header=T)
KERWdata[,2]<-as.numeric(substring(KERWdata[,2],7)) #strips the characters Event: from the Event ID column and turns the result numeric
KERWdata[,3]<-as.numeric(substring(KERWdata[,3],7)) #strips the characters Event: from the Event ID column and turns the result numeric
KERWdata[,4]<-as.numeric(substring(KERWdata[,4],14))
##Create igraph object
wKERs<-cbind(as.character(KERWdata[,2]),as.character(KERWdata[,3]))
AOPw<-graph_from_edgelist(wKERs, directed=T)
E(AOPw)$KER_ID<-KERWdata[,7]
E(AOPw)$evidence<-KERWdata[,6]
E(AOPw)$evidence[which(is.na(E(AOPw)$evidence))]<-3
E(AOPw)$quant<-KERWdata[,7]
E(AOPw)$quant[which(is.na(E(AOPw)$quant))]<-3
## Remove multiple edges
AOPws<-simplify(AOPw,remove.multiple=T,edge.attr.comb="min")
plot(AOPws, vertex.label=NA,vertex.size=2,edge.arrow.size=.05, edge.width=3/E(AOPws)$evidence)
## Use event ID # as vertex name attribute
V(AOPws)$name<-as.integer(V(AOPws)$name)
## Add KED data for MIE to AO analysis using weighted edges
V(AOPws)$KE_KED<-KEKdata[match(V(AOPws)$name,KEKdata[,2]),4] # finds KED (Key Event Designator) to add to V(AOPg) data
V(AOPws)$KE_KED[which(is.na(V(AOPws)$KE_KED))]<-"KE" # ALL KEs without KED (NA values from file) are assigned as generic KE
## Adds KED coloring to AOPws
V(AOPws)$ked_color<-"Yellow"
V(AOPws)$ked_color[which(V(AOPws)$KE_KED=="MIE")]<-"Green"
V(AOPws)$ked_color[which(V(AOPws)$KE_KED=="AO")]<-"Red"
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

### Import edge weight values from AOPws to AOPg ###

## Builds list of edge names that is used to match edges between AOPws and AOPg
wsEnames<-cbind(names(head_of(AOPws,E(AOPws))),names(tail_of(AOPws,E(AOPws))))
AOPgEnames<-cbind(names(head_of(AOPg,E(AOPg))),names(tail_of(AOPg,E(AOPg))))

## Creates a string from both edges to use in 'match' 
wsPairs<-paste(wsEnames[,1],wsEnames[,2])
AOPgPairs<-paste(AOPgEnames[,1],AOPgEnames[,2])

## Assigns evidence values from AOPws to AOPg graph
E(AOPg)$evidence<-E(AOPws)$evidence[charmatch(AOPgPairs,wsPairs)]
# Assigns NAs a value of 3 for evidence
E(AOPg)$evidence[is.na(E(AOPg)$evidence)]<-3
## Assigns quant values from AOPws to AOPg graph
E(AOPg)$quant<-E(AOPws)$quant[charmatch(AOPgPairs,wsPairs)]
# Assigns NAs a value of 3 for evidence
E(AOPg)$quant[is.na(E(AOPg)$quant)]<-3


## Set default plotting background color to black 
##!! Evaluate as either T or F or plots will not display properly
set.bg.black(T)


## Identifies which KEs are included in KERs, but are not themselves included in the KE event listings.
# V(AOPg)$KE_name[which(is.na(V(AOPg)$KE_EID))]



####~ AOPWIKI BY AOP ID VISUALIZATION ####

## Define a color palette for AOPwiki by AOP ID (usually lots of colors needed)
acols=colorRampPalette(c("green","red","cyan","orange","magenta","yellow","blue"))

##  Assign colors
for(i in 1:length(unique(V(AOPg)$AOP_ID))){
  V(AOPg)[which(V(AOPg)$AOP_ID==unique(V(AOPg)$AOP_ID)[i])]$acol<-acols(length(unique(V(AOPg)$AOP_ID)))[i]
  }

sort(table(V(AOPg)$AOP_ID))

## Highlight a given AOP for identification by using size on the network plot
V(AOPg)$exsize<-2
V(AOPg)[which(V(AOPg)$AOP_ID==130)]$exsize<-3 #Note: must specify V(AOPg)$exsize as vertex.size in plot for this to work

## Plot
set.seed(1)
plot(AOPg,layout=layout.fruchterman.reingold(AOPg), vertex.color=V(AOPg)$acol,vertex.label=NA, vertex.size=V(AOPg)$exsize, edge.arrow.size=.08)

## Calculates number of KE per unique AOP ID
AOP_freqs<-table(V(AOPg)$AOP_ID)
sort(AOP_freqs)
## Histogram of number of KE per unique AOP ID
hist(AOP_freqs,xlab="# Key Events",ylab="Frequency",col.axis=plotlabcol,col.lab=plotlabcol,col="gray")

## Barplot of number of KE per unique AOP ID with red line to show mean
bp_wcc<-barplot(table(V(AOPg)$AOP_ID),xaxt='n',xlab="AOP",ylab="# Key Events",col.axis=plotlabcol,col.lab=plotlabcol)
abline(h=mean(AOP_freqs),col="red")

#  TASK: WORK ON EDGE COLORING FOR AOP ID
#  edgecombcc<-expand.grid(V(gr)[which(comps$membership==ntcomps[i])],V(gr)[which(comps$membership==ntcomps[i])]) #creates a pairwise list of all nodes in the cc
#  edgecombflat<-as.vector(rbind(edgecombcc[[1]],edgecombcc[[2]])) #flattens the pairwise list to a vector where entries are read pairwise
#  edges.in.cc<-get.edge.ids(gr,edgecombflat,directed=TRUE)
#  E(gr)$color[edges.in.cc]<-cols[[i]]


####~ LEVEL OF BIOLOGICAL ORGANIZATION VISUALIZATION ####

## Add level of biological organization for key event nodes
V(AOPg)$lobo<-KEdata[[4]][match(V(AOPg)$KE_name,KEdata[[2]])]
V(AOPg)$lobo[which(is.na(V(AOPg)$lobo))]<-"" #assigns blank to NA data
tcols=rainbow(length(unique(V(AOPg)$lobo))) #creates a color scheme for visualization
lobo_list=c("Molecular","Cellular","Tissue","Organ","Individual","Population","") #creates an ordering of biological organization
V(AOPg)$lobo_o<-match(V(AOPg)$lobo,lobo_list) #assigns a value of biological organization instead of string.  1=molecular, 2=cellular, ...
lobo_freqs<-table(V(AOPg)$lobo_o)

## Plot the AOP wiki by lobo info
V(AOPg)$lobo_col<-tcols[V(AOPg)$lobo_o]
set.seed(1)
plot(AOPg ,vertex.size=2, edge.color="gray", edge.arrow.size=.1, vertex.label=NA, vertex.color=V(AOPg)$lobo_col)


## Plot the AOP wiki using a standard left to right lobo layout.
plot(AOPg, layout=lobo.layout(AOPg),vertex.size=2,  edge.curved=.3, edge.color="gray", edge.arrow.size=.1, vertex.label=NA, vertex.color=V(AOPg)$lobo_col)
# legend('topright',c("Molecular","Cellular","Tissue","Organ","Individual","Population","Not Specified"), pch=22,
# col="#777777", xjust=1,yjust=1, pt.bg=tcols, pt.cex=2, cex=.8, bty="n", ncol=1, y.intersp=.5, box.col="black", text.col="white")

## Barplot of lobo frequency
par(bg="white")
xx<- barplot(table(V(AOPg)$lobo_o), ylab="# of Key Events", xlab="Level of biological organization",col.axis=plotlabcol, col.lab=plotlabcol, col=tcols, axes=F,names.arg=NA)
text(x=xx, y=10, label=lobo_freqs, cex=.75)
legend('topright',c("Molecular","Cellular","Tissue","Organ","Individual","Population","Not Specified"), pch=22,
       col="#777777", xjust=1,yjust=1, pt.bg=tcols, pt.cex=2, cex=.8, bty="n", ncol=1, y.intersp=.5, box.col="white", text.col="black")

#### ~ MIE TO AO ANALYSES ####


## MIE, KE, and AO coloring assignments

V(AOPg)$ked_color<-"white"
V(AOPg)$ked_color[which(V(AOPg)$KE_KED=="MIE")]<-"Green"
V(AOPg)$ked_color[which(V(AOPg)$KE_KED=="AO")]<-"Red"


#### ~~ AOP edge connectivity ####
# This extends the graph wide edge connectivity calculation to
# only consider paths between MIE and AO as designated by V(graph)$KE_KED
ecl<-aop.edge.connectivity(AOPg,kelist=V(AOPg)$KE_KED)

# How many simple paths are there between MIEs and AOs in the wiki?
length(ecl)

ecl_sort<-ecl[order(ecl[,3],decreasing=T),]
colnames(ecl_sort)<-c("MIE","AO","ECL")

ecl_short<-head(ecl_sort,10)

ecl_short_name<-cbind(V(AOPg)$KE_name[ecl_short[,1]],V(AOPg)$KE_name[ecl_short[,2]],ecl_short[,3])
colnames(ecl_short_name)<-c("MIE","AO","ECL")
#name file to output to and create R variable for it
output<-file("results/ecl_results.txt")
#tell what to write into "output" file variable
write.table(ecl_short_name, output,col.names=T, sep="\t", row.names=F)
ecl_short
ecl_short_name

#### ~~~ MIE to AO Paths ####

## This part of the code can be used to identify all simple paths between
# a specified node-node pair.  These node-node pairs are usually MIE to AO
# as found in the aop.edge.connectivity analysis 
## Assign color to simple paths
E(AOPg)$asp_clr<-"gray"
E(AOPg)$asp_clr[simple.path.coloring(AOPg,57,677)]<-"purple"
## Assign sizes to simple paths
E(AOPg)$asp_size<-1
E(AOPg)$asp_size[simple.path.size(AOPg,57,677)]<-2
## Plot full network with simple paths between MIE and AO highlighted
set.seed(1)
plot(AOPg ,vertex.size=2, edge.width=E(AOPg)$asp_size, edge.color=E(AOPg)$asp_clr, edge.arrow.size=.1, vertex.label=NA, vertex.color=V(AOPg)$ked_color)
## Export network plot
jpeg.netplot(plot(AOPg ,vertex.size=2, edge.width=E(AOPg)$asp_size, edge.color=E(AOPg)$asp_clr, edge.arrow.size=.3, vertex.label=NA, vertex.color=V(AOPg)$ked_color),"f57t677",seedval=1)

#### ~~~~ MIE to AO Simple Path Subgraph ####
f57t677<-induced_subgraph(AOPg,unique(unlist(all_simple_paths(AOPg, from=57, to=677, mode="out"))))

## Plot subgraph (using ASP colors and sizes)
set.seed(1)
plot(f57t677,vertex.size=10, edge.width=E(f57t677)$asp_size, edge.color=E(f57t677)$asp_clr, edge.arrow.size=.15,
     vertex.label=V(f57t677)$name, vertex.color=V(f57t677)$ked_color)

## Export plot for subgraph
jpeg.netplot(plot(f57t677, vertex.size=10, edge.width=E(f57t677)$asp_size, edge.color=E(f57t677)$asp_clr, edge.arrow.size=1, vertex.label=V(f57t677)$name, vertex.color=V(f57t677)$ked_color),
             "f57t677_ev_sp",seedval=1,maii=c(0,0,0,0))
## Export topo plot for subgraph
jpeg.netplot(plot(f57t677,layout=topo.layout(f57t677),
                  vertex.label.degree=0, vertex.label.dist=1.18, edge.curved=1, vertex.label.color="white",
                  vertex.size=7, edge.width=E(f57t677)$asp_size, edge.color=E(f57t677)$asp_clr, edge.arrow.size=1,
                  vertex.label=V(f57t677)$KE_name, vertex.color=V(f57t677)$ked_color)
             ,"f57t677_topo",seedval=1,maii=c(0,0,0,2.1))

### Path analysis for subgraph ###

## Need to find the vertex ID's of the nodes from AOPg within the
## subgraph created in order to conduct path analyses
V(AOPg)[57]$name
V(AOPg)[677]$name
which(V(f57t677)$KE_EID==V(AOPg)[57]$KE_EID)
which(V(f57t677)$KE_EID==V(AOPg)[677]$KE_EID)



E(f57t677)$sp_cols<-"gray"
E(f57t677)$sp_cols[shortest.path.coloring(f57t677,f=4,t=17,weight=E(f57t677)$evidence)]<-"green"
E(f57t677)$sp_size<-1
E(f57t677)$sp_size[shortest.path.coloring(f57t677,f=4,t=17,weight=E(f57t677)$evidence)]<-2
set.seed(1)
plot(f57t677, vertex.size=10, edge.color=E(f57t677)$sp_cols, edge.arrow.size=.15, vertex.label=V(f57t677)$name, vertex.color=V(f57t677)$ked_color,edge.width=E(f57t677)$sp_size)

## Export plot for evidence weighted paths
jpeg.netplot(plot(f57t677, vertex.size=10, edge.color=E(f57t677)$sp_cols, edge.arrow.size=1, vertex.label=V(f57t677)$name, vertex.color=V(f57t677)$ked_color,edge.width=E(f57t677)$sp_size),
             "f57t677_ev_sp",seedval=1,maii=c(0,0,0,0))

## Export topo.plot for evidence weighted paths 
jpeg.netplot(plot(f57t677,layout=topo.layout(f57t677),vertex.label.degree=0, vertex.label.cex=1,vertex.label.dist=1.2, edge.curved=1, vertex.label.color="white", vertex.size=7, edge.color=E(f57t677)$sp_cols, edge.arrow.size=1, vertex.label=V(f57t677)$KE_name, vertex.color=V(f57t677)$ked_color,edge.width=E(f57t677)$sp_size),
             "f57t677_ev_sp_topo",seedval=1,maii=c(0,0,0,2.1))

### Quantitative Understanding Short

E(f57t677)$sp_cols<-"gray"
E(f57t677)$sp_cols[shortest.path.coloring(f57t677,f=4,t=17,weight=E(f57t677)$quant)]<-"blue"
E(f57t677)$sp_size<-1
E(f57t677)$sp_size[shortest.path.coloring(f57t677,f=4,t=17,weight=E(f57t677)$quant)]<-2
set.seed(1)
plot(f57t677, vertex.size=10, edge.color=E(f57t677)$sp_cols, edge.arrow.size=.15, vertex.label=V(f57t677)$name, vertex.color=V(f57t677)$ked_color,edge.width=E(f57t677)$sp_size)

## Export plot for quantitative understanding weighted paths
jpeg.netplot(plot(f57t677, vertex.size=10, edge.color=E(f57t677)$sp_cols, edge.arrow.size=1, vertex.label=V(f57t677)$name, vertex.color=V(f57t677)$ked_color,edge.width=E(f57t677)$sp_size),
             "f57t677_qu_sp",seedval=1,maii=c(0,0,0,0))

## Export topo.plot for evidence weighted paths 
jpeg.netplot(plot(f57t677,layout=topo.layout(f57t677),vertex.label.degree=0, vertex.label.cex=1,vertex.label.dist=1.2, edge.curved=1, vertex.label.color="white", vertex.size=7, edge.color=E(f57t677)$sp_cols, edge.arrow.size=1, vertex.label=V(f57t677)$KE_name, vertex.color=V(f57t677)$ked_color,edge.width=E(f57t677)$sp_size),
             "f57t677_qu_sp_topo",seedval=1,maii=c(0,0,0,2.1))

## Non edge-weighted shortest path
E(f57t677)$sp_cols<-"gray"
E(f57t677)$sp_cols[shortest.path.coloring(f57t677,f=4,t=17)]<-"orange"
E(f57t677)$sp_size<-1
E(f57t677)$sp_size[shortest.path.coloring(f57t677,f=4,t=17)]<-2
set.seed(1)
plot(f57t677, vertex.size=10, edge.color=E(f57t677)$sp_cols, edge.arrow.size=.15, vertex.label=V(f57t677)$name, vertex.color=V(f57t677)$ked_color,edge.width=E(f57t677)$sp_size)

## Export plot for quantitative understanding weighted paths
jpeg.netplot(plot(f57t677, vertex.size=10, edge.color=E(f57t677)$sp_cols, edge.arrow.size=1, vertex.label=V(f57t677)$name, vertex.color=V(f57t677)$ked_color,edge.width=E(f57t677)$sp_size),
             "f57t677_nw_sp",seedval=1,maii=c(0,0,0,0))


## Export topo.plot for evidence weighted paths 
jpeg.netplot(plot(f57t677,layout=topo.layout(f57t677),vertex.label.degree=0, vertex.label.cex=1,vertex.label.dist=1.2, edge.curved=1, vertex.label.color="white", vertex.size=7, edge.color=E(f57t677)$sp_cols, edge.arrow.size=1, vertex.label=V(f57t677)$KE_name, vertex.color=V(f57t677)$ked_color,edge.width=E(f57t677)$sp_size),
             "f57t677_nw_sp_topo",seedval=1,maii=c(0,0,0,2.1))

#### ~~ AOP Path Occurrence Visualization NEEDS REFINEMENT ####
V(AOPg)$aopp<-aop.paths(AOPg,kelist=V(AOPg)$KE_KED)
heatgrad=rev(heat.colors(n=max(V(AOPg)$aopp)))
heatgrad<-append(heatgrad, "#FFFFFF", after=0) #in the instance that a KE is not included in paths between MIE and AO...

wbpal=colorRampPalette(c("white","blue"))(n=max(V(AOPg)$aopp)+1)

V(AOPg)$po_col<-wbpal[V(AOPg)$aopp+1]

set.seed(1)
plot(AOPg ,vertex.size=2, edge.color=E(AOPg)$asp_clr, edge.arrow.size=.05, vertex.label=NA, vertex.color=V(AOPg)$po_col)

####~ CONNECTED COMPONENTS ANALYSIS####

####~~ Weakly connected component analysis and plotting ####
## Color vertices and edges by their weakly connected components 
V(AOPg)$cc_color<-unlist(color.comps(AOPg,"weak")$vcol)  #color.comps is a custom function stored in the AOP_net_functions.R file
E(AOPg)$cc_color<-unlist(color.comps(AOPg,"weak")$ecol)  #color.comps is a custom function stored in the AOP_net_functions.R file

## Plot
set.seed(1)
plot(AOPg,vertex.size=2, edge.arrow.size=.1,vertex.color=V(AOPg)$cc_color, edge.color=E(AOPg)$cc_color,  vertex.size=2,vertex.label=NA)


## barplot for size of weakly connected components
wcomps<-components(AOPg,mode="weak")
wcc_freqs<-table(wcomps$csize)
bp_wcc<-barplot(table(wcomps$csize),col=plotlabcol, col.axis=plotlabcol, xlab="Component size",ylab="Frequency",col.lab=plotlabcol)


####~~ Strongly connected component analysis and plotting ####

## When the "strong" option is passed to color.comps, vsize and ewidth are calculated and can be used within plot
V(AOPg)$cc_color<-unlist(color.comps(AOPg,"strong")$vcol)  #color.comps is a custom function stored in the AOP_net_functions.R file
E(AOPg)$cc_color<-unlist(color.comps(AOPg,"strong")$ecol)  #color.comps is a custom function stored in the AOP_net_functions.R file

V(AOPg)$cc_size<-unlist(color.comps(AOPg,"strong")$vsize)
E(AOPg)$cc_width<-unlist(color.comps(AOPg,"strong")$ewidth)

## Plot of connected components with strong sizing option
par(bg="white")
set.seed(1)
plot(AOPg,vertex.size=V(AOPg)$cc_size, edge.width=E(AOPg)$cc_width, vertex.color=V(AOPg)$cc_color, edge.color=E(AOPg)$cc_color,  vertex.size=2, edge.arrow.size=.1, vertex.label=NA)

## This points out how many of the feedback loops/cycles are contained within the same AOP and how many are a result of the network
scomps<-components(AOPg,mode="strong")
ntcomps<-which(scomps$csize>1) # non-trivial ccs (i.e. with more than 1 node)
V(AOPg)$scc<-scomps$membership # assign the attribute scc to nodes based on their membership 
for(i in 1:length(ntcomps)){
  print(V(AOPg)[which(V(AOPg)$scc==ntcomps[i])]$AOP_ID)
}
for(i in 1:length(ntcomps)){
  print(V(AOPg)[which(V(AOPg)$scc==ntcomps[i])]$name)
}



####~ CENTRALITY MEASURES FOR THE AOPWIKI ####



####~~ Degree centrality ####




####~~~ Total Degree ####

## Which key events have the most incident nodes?
sort(degree(AOPg, mode="all"))

## Prints the top ten key events based on total degree 
rev(as.integer(tail(sort(degree(AOPg, mode="all")),10)))
as.integer(names(tail(sort(degree(AOPg, mode="all")),10)))
 rev(V(AOPg)$KE_name[match(as.integer(names(tail(sort(degree(AOPg, mode="all")),10))),V(AOPg)$name)])

## Assigns node coloring by degree
V(AOPg)$deg_col<-deg.col.grad(AOPg,dmode="all",totdegcol)

## Plot colored and nodes sized by degree
set.seed(1)
plot(AOPg, vertex.size=500*degree(AOPg,mode="all",normalized=TRUE), vertex.color=V(AOPg)$deg_col, edge.arrow.size=.1, vertex.label=NA, edge.color="gray",edge.width=1)

## Barplot to show degree histogram
b2gpal=colorRampPalette(totdegcol)
barplot(table(degree(AOPg,mode="all")), xlab="Total Degree", ylab="Frequency",col.axis=plotlabcol, col.lab=plotlabcol,col=b2gpal(max(degree(AOPg,mode="all"))))
legend('topright',legend=rev(seq(min(degree(AOPg,mode="all")),max(degree(AOPg,mode="all")),3)), pch=22,
       col="#777777", xjust=1,yjust=1, pt.bg=rev(b2gpal(length(rev(seq(min(degree(AOPg,mode="all")),max(degree(AOPg,mode="all")),3))))), pt.cex=2, cex=.8, bty="n", ncol=1, y.intersp=.5, box.col="white", text.col=plotlabcol)




####~~~ In-Degree ####

## Prints the top-ten key event names by in-degree value
sort(degree(AOPg,mode="in"))
rev(as.integer(tail(sort(degree(AOPg, mode="in")),10)))
rev(V(AOPg)$KE_name[match(as.integer(names(tail(sort(degree(AOPg, mode="in")),10))),V(AOPg)$name)])

## Assigns node coloring by degree
V(AOPg)$deg_col<-deg.col.grad(AOPg,dmode="in",indegcol)

## Plot colored and nodes sized by degree
set.seed(1)
plot(AOPg, vertex.size=500*degree(AOPg,mode="in",normalized=TRUE), vertex.color=V(AOPg)$deg_col, edge.arrow.size=.1, vertex.label=NA, edge.color="gray",edge.width=1)

## Barplot to show degree histogram
b2ppal=colorRampPalette(indegcol)
barplot(table(degree(AOPg,mode="in")), xlab="In-Degree", ylab="Frequency",col.axis=plotlabcol, col.lab=plotlabcol,col=b2ppal(max(degree(AOPg,mode="in"))))
legend('topright',legend=rev(seq(min(degree(AOPg,mode="in")),max(degree(AOPg,mode="in")),3)), pch=22,
       col="#777777", xjust=1,yjust=1, pt.bg=rev(b2ppal(length(rev(seq(min(degree(AOPg,mode="in")),max(degree(AOPg,mode="in")),3))))), pt.cex=2, cex=.8, bty="n", ncol=1, y.intersp=.5, box.col="white", text.col=plotlabcol)




####~~~ Out-Degree ####

## Prints the top-ten key event names by out-degree value
sort(degree(AOPg,mode="out"))
rev(as.integer(tail(sort(degree(AOPg, mode="out")),10)))
rev(V(AOPg)$KE_name[match(as.integer(names(tail(sort(degree(AOPg, mode="out")),10))),V(AOPg)$name)])

## Colors nodes by out-degree value
V(AOPg)$deg_col<-deg.col.grad(AOPg,dmode="out",outdegcol)

## Plot colored and nodes sized by degree
set.seed(1)
plot(AOPg, vertex.size=500*degree(AOPg,mode="out",normalized=TRUE), vertex.color=V(AOPg)$deg_col, edge.arrow.size=.1, vertex.label=NA, edge.color="gray",edge.width=1)

## Barplot to show degree histogram
b2opal=colorRampPalette(outdegcol)
barplot(table(degree(AOPg,mode="out")), xlab="Out-Degree", ylab="# of KEs with Out-Degree",col.axis=plotlabcol, col.lab=plotlabcol,col=b2opal(max(degree(AOPg,mode="out"))))
legend('topright',legend=rev(seq(min(degree(AOPg,mode="out")),max(degree(AOPg,mode="out")),3)), pch=22,
       col="#777777", xjust=1,yjust=1, pt.bg=rev(b2opal(length(rev(seq(min(degree(AOPg,mode="out")),max(degree(AOPg,mode="out")),3))))), pt.cex=2, cex=.8, bty="n", ncol=1, y.intersp=.5, box.col="white", text.col=plotlabcol)




####~~ Betweenness Centrality ####

## Prints the top-ten key event names by betweenness value
sort(betweenness(AOPg,directed=TRUE))
rev(as.integer(tail(sort(betweenness(AOPg,directed=TRUE)),10)))
rev(V(AOPg)$KE_name[match(as.integer(names(tail(sort(betweenness(AOPg,directed=TRUE)),10))),V(AOPg)$name)])

## Colors nodes by betweenness values
b2upal=colorRampPalette(betcol)
V(AOPg)$bet_col<-b2upal(20)[as.numeric(cut(betweenness(AOPg,directed=TRUE),breaks = 20))]

## Plot colored and nodes sized by degree
set.seed(1)
plot(AOPg, vertex.size=3000*betweenness(AOPg,normalized=TRUE,directed=T), vertex.color=V(AOPg)$bet_col, edge.arrow.size=.1, vertex.label=NA, edge.color="gray",edge.width=1)

## Betweenness histogram
hist(log10(betweenness(AOPg)),breaks=20,col=b2upal(20), col.axis=plotlabcol, col.lab=plotlabcol, xlab=expression(paste("Log"[10],"(Betweeness"["dir"],")","")), main="")
## Barplot of betweenness totals with 20 bins
barplot(table(as.numeric(cut(betweenness(AOPg),breaks = 20))),col=b2upal(20),col.axis=plotlabcol, col.lab=plotlabcol)

####~~~ Undirected Betweenness ####

## Prints the top-ten key event names by betweenness value
sort(betweenness(AOPg,directed=F))
rev(as.integer(tail(sort(betweenness(AOPg,directed=F)),10)))
rev(V(AOPg)$KE_name[match(as.integer(names(tail(sort(betweenness(AOPg,directed=F)),10))),V(AOPg)$name)])

## Colors nodes by betweenness values
b2upal=colorRampPalette(betcol)
V(AOPg)$bet_col<-b2upal(30)[as.numeric(cut(betweenness(AOPg,directed=F),breaks = 30))]

## Plot colored and nodes sized by degree
set.seed(1)
plot(AOPg, vertex.size=50*betweenness(AOPg,normalized=TRUE,directed=F), vertex.color=V(AOPg)$bet_col, edge.arrow.size=.1, vertex.label=NA, edge.color="gray",edge.width=1)

## Betweenness histogram
hist(log10(betweenness(AOPg, directed=F)),breaks=30,col=b2upal(30), col.axis=plotlabcol, col.lab=plotlabcol, xlab=expression(paste("Log"[10],"(Betweeness"["undir"],")","")), main="")

hist(betweenness(AOPg, directed=F),breaks=20,col=b2upal(20), col.axis=plotlabcol, col.lab=plotlabcol, xlab=expression(paste("Log"[10],"(Betweeness)","")), main="")
## Barplot of betweenness totals with 20 bins
barplot(table(as.numeric(cut(betweenness(AOPg,directed=F),breaks = 30))),col=b2upal(30),col.axis=plotlabcol, col.lab=plotlabcol)




####~~ Closeness Analysis ####
## Prints the top-ten key event names by closeness value
sort(closeness(AOPg,mode="all"))
mostclsallval<-round(rev((tail(sort(closeness(AOPg,mode="all")*10^6),5))),5)
mostclsall<-rev(V(AOPg)$KE_name[match(as.integer(names(tail(sort(closeness(AOPg,mode="all")),5))),V(AOPg)$name)])

clsall<-paste(mostclsall," (",mostclsallval,")",sep="")
print(clsall)
#name file to output to and create R variable for it
output<-file("results/cls_all_results.txt")
#tell what to write into R variable
writeLines(clsall, output)
#Stop close file being edited
close(output)

## Colors nodes based on cloesness values
b2mpal=colorRampPalette(clscol)
V(AOPg)$close_col<-b2mpal(1500)[as.numeric(cut(closeness(AOPg,mode="all",norm=T),breaks = 1500))]

## Plots color and nodes sized by closeness
set.seed(1)
plot(AOPg, vertex.size=1500*closeness(AOPg,normalized=T,mode="all"), vertex.color=V(AOPg)$close_col, edge.arrow.size=.1, vertex.label=NA)#, vertex.color="orange",edge.color="gray")

## Histogram of closeness (including a change of units)
hist(closeness(AOPg,mode="all")*10^6,breaks=20,col=b2mpal(20))
## Scatterplot of closeness values
plot(closeness(AOPg,mode="all"), xlab="Key Event", col.axis=plotlabcol, col.lab=plotlabcol, xaxt='n', ylab="Closeness Value",main="KE Total Closeness in AOPwiki",col=V(AOPg)$close_col)

####~~~In-path Closeness ####

clsmode="in"
## Prints the top-ten key event names by closeness value
sort(closeness(AOPg,mode=clsmode))
mostclsinval<-round(rev(tail(sort(closeness(AOPg,mode=clsmode)*10^6),5)),5)
mostclsin<-rev(V(AOPg)$KE_name[match(as.integer(names(tail(sort(closeness(AOPg,mode=clsmode)),5))),V(AOPg)$name)])

clsin<-paste(mostclsin," (",mostclsinval,")",sep="")
print(clsin)
#name file to output to and create R variable for it
output<-file("results/cls_in_results.txt")
#tell what to write into R variable
writeLines(clsin, output)
#Stop close file being edited
close(output)


## Colors nodes based on cloesness values
b2mpal=colorRampPalette(clscol)
V(AOPg)$close_col<-b2mpal(1500)[as.numeric(cut(closeness(AOPg,mode=clsmode,norm=TRUE),breaks = 1500))]

## Plots color and nodes sized by closeness
set.seed(1)
plot(AOPg, vertex.size=1500*closeness(AOPg,normalized=TRUE,mode=clsmode), vertex.color=V(AOPg)$close_col, edge.arrow.size=.1, vertex.label=NA)#, vertex.color="orange",edge.color="gray")

## Histogram of closeness (including a change of units)
hist(closeness(AOPg,mode=clsmode)*10^6,breaks=20,col=b2mpal(20),col.axis=plotlabcol, col.lab=plotlabcol)
## Scatterplot of closeness values
plot(closeness(AOPg,mode=clsmode), xlab="Key Event", col.axis=plotlabcol, col.lab=plotlabcol, xaxt='n', ylab="Closeness Value",main="KE In-Closeness in AOPwiki",col=V(AOPg)$close_col)

####~~~Out-path Closeness ####

clsmode="out"
## Prints the top-ten key event names by closeness value
sort(closeness(AOPg,mode=clsmode))
mostclsoutval<-round(rev(tail(sort(closeness(AOPg,mode=clsmode)*10^6),5)),5)
mostclsout<-rev(V(AOPg)$KE_name[match(as.integer(names(tail(sort(closeness(AOPg,mode=clsmode)),5))),V(AOPg)$name)])

clsout<-paste(mostclsout," (",mostclsoutval,")",sep="")
print(clsout)
#name file to output to and create R variable for it
output<-file("results/cls_out_results.txt")
#tell what to write into R variable
writeLines(clsout, output)
#Stop close file being edited
close(output)

## Colors nodes based on cloesness values
b2mpal=colorRampPalette(clscol)
V(AOPg)$close_col<-b2mpal(1500)[as.numeric(cut(closeness(AOPg,mode=clsmode,norm=T),breaks = 1500))]

## Plots color and nodes sized by closeness
set.seed(1)
plot(AOPg, vertex.size=1500*closeness(AOPg,normalized=T,mode=clsmode), vertex.color=V(AOPg)$close_col, edge.arrow.size=.1, vertex.label=NA)#, vertex.color="orange",edge.color="gray")

## Histogram of closeness (including a change of units)
hist(closeness(AOPg,mode=clsmode)*10^6,breaks=20,col=b2mpal(20),col.axis=plotlabcol, col.lab=plotlabcol)
## Scatterplot of closeness values
plot(closeness(AOPg,mode=clsmode), xlab="Key Event", col.axis=plotlabcol, col.lab=plotlabcol, xaxt='n', ylab="Closeness Value",main="KE Out-Closeness in AOPwiki",col=V(AOPg)$close_col)

## Summary File ##

cls<-paste(mostclsall," (",mostclsallval,") & ",mostclsin," (",mostclsinval,") & ",mostclsout," (",mostclsoutval,")",sep="")
print(cls)
#name file to output to and create R variable for it
output<-file("results/cls_results.txt")
#tell what to write into R variable
writeLines(cls, output)
#Stop close file being edited
close(output)



####~~ Eccentricity ####


####~~~ Total Eccentricity ####
eccmode="all"
## Assigns eccentricity values as a node attribute
V(AOPg)$ecc<-eccentricity(AOPg,mode =eccmode)

## Prints the top-ten key event names by eccentricity value
sort(eccentricity(AOPg,mode =eccmode))
mosteccallval<-rev(as.integer(tail(sort(eccentricity(AOPg,mode=eccmode)),10)))
mosteccall<-rev(V(AOPg)$KE_name[match(as.integer(names(tail(sort(eccentricity(AOPg,mode=eccmode)),10))),V(AOPg)$name)])
eccall<-paste(mosteccall," (",mosteccallval,")",sep="")
print(eccall)
#name file to output to and create R variable for it
output<-file("results/ecc_all_results.txt")
#tell what to write into R variable
writeLines(eccall, output)
#Stop close file being edited
close(output)


## Assigns color based on eccentricity value
b2cpal=colorRampPalette(ecccol)
V(AOPg)$ecc_col<-b2cpal(max(V(AOPg)$ecc)+1)[V(AOPg)$ecc+1]

## Plot colors and node sizes based on eccentricity value
set.seed(1)
plot(AOPg, vertex.size=4*(V(AOPg)$ecc/max(V(AOPg)$ecc)), vertex.color=V(AOPg)$ecc_col, edge.arrow.size=.1, vertex.label=NA, edge.color="gray",edge.width=1)

## Histogram of eccentricity values
barplot(table(eccentricity(AOPg,mode=eccmode)), xlab="Total Eccentricity", ylab="Frequency",col.axis=plotlabcol, col.lab=plotlabcol,col=b2cpal(max(eccentricity(AOPg,mode =eccmode))))
legend('topright',legend=rev(seq(min(eccentricity(AOPg,mode =eccmode)),max(eccentricity(AOPg,mode =eccmode)),3)), pch=22,
       col="#777777", xjust=1,yjust=1, pt.bg=rev(b2cpal(length(rev(seq(min(eccentricity(AOPg,mode =eccmode)),max(eccentricity(AOPg,mode =eccmode)),3))))), pt.cex=2, cex=.8, bty="n", ncol=1, y.intersp=.5, box.col="white", text.col=plotlabcol)

####~~~ In-Eccentricity ####

eccmode="in"
## Assigns eccentricity values as a node attribute
V(AOPg)$ecc<-eccentricity(AOPg,mode =eccmode)

## Prints the top-ten key event names by eccentricity value
sort(eccentricity(AOPg,mode =eccmode))
mosteccinval<-rev(as.integer(tail(sort(eccentricity(AOPg,mode=eccmode)),10)))
mosteccin<-rev(V(AOPg)$KE_name[match(as.integer(names(tail(sort(eccentricity(AOPg,mode=eccmode)),10))),V(AOPg)$name)])
eccin<-paste(mosteccin," (",mosteccinval,")",sep="")
print(eccin)

#name file to output to and create R variable for it
output<-file("results/ecc_in_results.txt")
#tell what to write into R variable
writeLines(eccin, output)
#Stop close file being edited
close(output)

## Assigns color based on eccentricity value
b2cpal=colorRampPalette(ecccol)
V(AOPg)$ecc_col<-b2cpal(max(V(AOPg)$ecc)+1)[V(AOPg)$ecc+1]

## Plot colors and node sizes based on eccentricity value
set.seed(1)
plot(AOPg, vertex.size=4*(V(AOPg)$ecc/max(V(AOPg)$ecc)), vertex.color=V(AOPg)$ecc_col, edge.arrow.size=.1, vertex.label=NA, edge.color="gray",edge.width=1)

## Histogram of eccentricity values
barplot(table(eccentricity(AOPg,mode=eccmode)), xlab="In-Eccentricity", ylab="Frequency",col.axis=plotlabcol, col.lab=plotlabcol,col=b2cpal(max(eccentricity(AOPg,mode =eccmode)+1)))
legend('topright',legend=rev(seq(min(eccentricity(AOPg,mode =eccmode)),max(eccentricity(AOPg,mode =eccmode))+1,3)), pch=22,
       col="#777777", xjust=1,yjust=1, pt.bg=rev(b2cpal(length(rev(seq(min(eccentricity(AOPg,mode =eccmode)),max(eccentricity(AOPg,mode =eccmode)+1),3))))), pt.cex=2, cex=.8, bty="n", ncol=1, y.intersp=.5, box.col="white", text.col=plotlabcol)

####~~~ Out-Eccentricity ####

eccmode="out"
## Assigns eccentricity values as a node attribute
V(AOPg)$ecc<-eccentricity(AOPg,mode =eccmode)

## Prints the top-ten key event names by eccentricity value
sort(eccentricity(AOPg,mode =eccmode))
mosteccoutval<-rev(as.integer(tail(sort(eccentricity(AOPg,mode=eccmode)),10)))
mosteccout<-rev(V(AOPg)$KE_name[match(as.integer(names(tail(sort(eccentricity(AOPg,mode=eccmode)),10))),V(AOPg)$name)])
eccout<-paste(mosteccout," (",mosteccoutval,")",sep="")
print(eccout)

#name file to output to and create R variable for it
output<-file("results/ecc_out_results.txt")
#tell what to write into R variable
writeLines(eccout, output)
#Stop close file being edited
close(output)
## Assigns color based on eccentricity value
b2cpal=colorRampPalette(ecccol)
V(AOPg)$ecc_col<-b2cpal(max(V(AOPg)$ecc)+1)[V(AOPg)$ecc+1]

## Plot colors and node sizes based on eccentricity value
set.seed(1)
plot(AOPg, vertex.size=4*(V(AOPg)$ecc/max(V(AOPg)$ecc)), vertex.color=V(AOPg)$ecc_col, edge.arrow.size=.1, vertex.label=NA, edge.color="gray",edge.width=1)

## Histogram of eccentricity values
barplot(table(eccentricity(AOPg,mode=eccmode)), xlab="Out-Eccentricity", ylab="Frequency",col.axis=plotlabcol, col.lab=plotlabcol,col=b2cpal(max(eccentricity(AOPg,mode =eccmode))+1))
legend('topright',legend=rev(seq(min(eccentricity(AOPg,mode =eccmode)),max(eccentricity(AOPg,mode =eccmode)),3)), pch=22,
       col="#777777", xjust=1,yjust=1, pt.bg=rev(b2cpal(length(rev(seq(min(eccentricity(AOPg,mode =eccmode)),max(eccentricity(AOPg,mode =eccmode)),3))))), pt.cex=2, cex=.8, bty="n", ncol=1, y.intersp=.5, box.col="white", text.col=plotlabcol)

## Summary File ##
ecc<-paste(mosteccall," (",mosteccallval,") & ",mosteccin," (",mosteccinval,") & ",mosteccout," (",mosteccoutval,")",sep="")
print(ecc)

#name file to output to and create R variable for it
output<-file("results/ecc_results.txt")
#tell what to write into R variable
writeLines(ecc, output)
#Stop close file being edited
close(output)


####~ ADDITIONAL ANALYSES ####

####~~ Reciprocity ####
# The measure of reciprocity defines the proportion 
#of mutual connections, in a directed graph. 
#It is most commonly defined as the probability that the 
#opposite counterpart of a directed edge is also included 
#in the graph. Or in adjacency matrix notation: 
#sum(i, j, (A.*A')ij) / sum(i, j, Aij), where A.*A' 
#is the element-wise product of matrix A and its transpose. 
#This measure is calculated if the mode argument is default.
#In the AOP framework this will identify the proportion of nodes with
#KE to KE pair feedbacks.

r<-reciprocity(AOPg)
## Gives number of 2-cycles in the graph
length(E(AOPg))*reciprocity(AOPg)/2
## This can be used to find which nodes are involved in the 2cycles.
AA<-as_adjacency_matrix(AOPg)
A<-as.matrix(AA)
At<-t(A)
which(At*A!=0, arr.ind=T)
# rho is a modification of reciprocity presented by 
# Garlaschelli and Loffredo, 2004 to improve on standard 
# reciprocity calculations, -1<rho<1 rho=1=> perfectly reciprocal
# rho=0=>areciprocal, rho=-1=> perfectly antireciprocal
L<-length(which(A!=0))
N<-length(V(AOPg))
ab<-L/(N*(N-1))
rho<-(r-ab)/(1-ab)
rho




####~~ Full Network Edge Connectivity ####
# # The edge connectivity of a pair of vertices (source and target) 
# # is the minimum number of edges needed to remove to eliminate 
# # all (directed) paths from source to target.
# 
# str(V(AOPg))
# # edge_connectivity doesn't work on a full network unless it is connected
# edge_connectivity(AOPg, source = NULL, target = NULL, checks = TRUE)
# 
# #In order to run this analysis, the graph can be decomposed
# dg<-decompose.graph(AOPg)
# edge_connectivity(dg[[3]],source=NULL, target=NULL, checks=TRUE)
# 
# wcc1<-induced_subgraph(AOPg,subcomponent(AOPg,1))
# edge_connectivity(wcc1,checks=TRUE)


####~~ Graph condensation ####
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
