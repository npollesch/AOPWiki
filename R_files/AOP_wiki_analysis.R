
#### LOAD PACKAGES AND IMPORT AOPWIKI DATA ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

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

## Identify all unique KEs by looking at all names involved in KERs
allKEs<-append(KERdata[,2],KERdata[,4])
uniqueKEs<-unique(allKEs)


keID<-data.frame(ID=1:length(uniqueKEs),KE=uniqueKEs)
KERs<-cbind(match(KERdata[,2], keID[,2]),match(KERdata[,4], keID[,2]))
AOPg<-graph_from_edgelist(KERs, directed=T) #This construction has unnessary steps

## Add names for key event nodes
V(AOPg)$KE_name<-as.character(keID$KE)
## Add event IDS
V(AOPg)$KE_EID<-KEdata[match(V(AOPg)$KE_name,KEdata[,3]),1] # adds event ID number
V(AOPg)$name<-V(AOPg)$KE_EID #changes default 'name' object to be the KE_EID

## Add AOP IDs. Note: Each KE may belong to more than one AOP, so the AOP_IDs object is a list not just a single AOP ID
for(i in 1:length(V(AOPg))){
  if(length(which(!is.na(match(KEPdata[,3],V(AOPg)$KE_EID[i]))))>0){
    V(AOPg)[i]$AOP_ID<-list(unique(KEPdata[which(!is.na(match(KEPdata[,3],V(AOPg)[i]$KE_EID))),1]))}
else{V(AOPg)[i]$AOP_ID<-NA}
}

## Add key event designation data
V(AOPg)$KE_KED<-KEKdata[match(V(AOPg)$KE_EID,KEKdata[,2]),4] # finds KED (Key Event Designator) to add to V(AOPg) data
length(V(AOPg)$KE_KED[which(is.na(V(AOPg)$KE_KED))]) #The number of KEs without KEDs
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
#E(AOPw)$evidence[which(is.na(E(AOPw)$evidence))]<-3
E(AOPw)$quant<-KERWdata[,7]
#E(AOPw)$quant[which(is.na(E(AOPw)$quant))]<-3

## Remove multiple edges
AOPws<-simplify(AOPw,remove.multiple=T,edge.attr.comb="min")
#plot(AOPws, vertex.label=NA,vertex.size=2,edge.arrow.size=.05, edge.width=3/E(AOPws)$evidence)

## Use event ID # as vertex name attribute
V(AOPws)$name<-as.integer(V(AOPws)$name)

## Add KED data for MIE to AO analysis using weighted edges
V(AOPws)$KE_KED<-KEKdata[match(V(AOPws)$name,KEKdata[,2]),4] # finds KED (Key Event Designator) to add to V(AOPg) data
length(V(AOPws)$KE_KED[which(is.na(V(AOPws)$KE_KED))]) #The number of KEs without KEDs
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
length(E(AOPg)$evidence[!is.na(E(AOPg)$evidence)])
mean(E(AOPg)$evidence[!is.na(E(AOPg)$evidence)])
# Assigns NAs a value of 3 for evidence
E(AOPg)$evidence[is.na(E(AOPg)$evidence)]<-3
## Assigns quant values from AOPws to AOPg graph
E(AOPg)$quant<-E(AOPws)$quant[charmatch(AOPgPairs,wsPairs)]
length(E(AOPg)$quant[!is.na(E(AOPg)$quant)])
mean(E(AOPg)$quant[!is.na(E(AOPg)$quant)])
# Assigns NAs a value of 3 for evidence
E(AOPg)$quant[is.na(E(AOPg)$quant)]<-3

### AOP ID from AOPws 
V(AOPws)$KE_KED<-KEKdata[match(V(AOPws)$name,KEKdata[,2]),4] # finds KED (Key Event Designator) to add to V(AOPg) data
length(V(AOPws)$KE_KED[which(is.na(V(AOPws)$KE_KED))]) # The number of KEs without KEDs
V(AOPws)$KE_KED[which(is.na(V(AOPws)$KE_KED))]<-"KE" # ALL KEs without KED (NA values from file) are assigned as generic KE

KEKdata[,1]<-as.numeric(substring(KEKdata[,1],5))

## MIE, KE, and AO coloring assignments
V(AOPg)$ked_color<-"white"
V(AOPg)$ked_color[which(V(AOPg)$KE_KED=="MIE")]<-"Green"
V(AOPg)$ked_color[which(V(AOPg)$KE_KED=="AO")]<-"Red"


## ASSIGN AOP IDS TO AOPWS VERTICES #NEEDS TO BE FIXED# ##
#V(AOPws)$AOP_ID<-KEKdata[match(V(AOPws)$name,KEKdata[,2]),1]
#V(AOPws)$AOP_ID



## Set default plotting background color to black 
##!! Evaluate as either T or F or plots will not display properly
set.bg.black(F)

#### NETWORK SUMMARY ####

## AOP ID Quant

dev.off()
jpeg(file = paste("images/","aop_id_occ",".jpeg",sep=""),
     width=800, height=800, bg =bg_col, quality=100)
plot(table(sort(unlist(V(AOPg)$AOP_ID))),xlab="AOP ID",ylab="Number of KEs")
abline(h=6.75, col="blue")
mtext("Mean \n 6.75",side=2,col="blue", at=6.75)
dev.off()
     
head(sort(table(unlist(V(AOPg)$AOP_ID))))
tail(sort(table(unlist(V(AOPg)$AOP_ID))))

## KED Quant
table(V(AOPg)$KE_KED)

V(AOPg)$ked_color<-"white"
V(AOPg)$ked_color[which(V(AOPg)$KE_KED=="MIE")]<-"Green"
V(AOPg)$ked_color[which(V(AOPg)$KE_KED=="AO")]<-"Red"

## Average number of AOP IDS per KE
length(unlist(V(AOPg)$AOP_ID))/750

##PLOT OF WIKI BASED ON KED
aplot(AOPg,vcol=V(AOPg)$ked_color,vsize=2)
legend("bottomright",c("MIE","KE","AO"),pch=c(16,1,16),col=c("Green","black","red"),y.intersp =1)

## Attempting to get at connectivity growth and scaling of the AOP wiki (similar to Strogatz and other nature/science paper (Albert barbasi??))
plot(table(degree(AOPg)))
comps<-components(AOPg)

plot(rev(log10(sort(comps$csize)/750)))

sum(match(comps$membership,1,nomatch=0))

for(i in 1:34){
  csize[i]*sum(match(comps$membership,i,nomatch=0))
}

table(degree(AOPg))/sum(degree(AOPg))
plot(log10(degree.distribution((AOPg))))
curve(x^(-3),from=0, to=22)


####~ Quality Control ####

## All the nodes without event ID Numbers
as.numeric(V(AOPg)[is.na(V(AOPg)$KE_EID)])

## KED is unique in that if an actual designation doesn't exist in the WIKI, we can
## simply apply the generic "KE" label, however it is useful to know how
## many KEs do not have KEDs from the wiki 
## NOTE: This code is used in buliding the network above
V(AOPg)$KE_KED<-KEKdata[match(V(AOPg)$KE_EID,KEKdata[,2]),4] # finds KED (Key Event Designator) to add to V(AOPg) data
length(V(AOPg)$KE_KED[which(is.na(V(AOPg)$KE_KED))]) #The number of KEs without KEDs
V(AOPg)$KE_KED[which(is.na(V(AOPg)$KE_KED))]<-"KE" # ALL KEs without KED (NA values from file) are assigned as generic KE

## Find the KEs without any AOP designation
## Event IDs for KEs with no AOP
V(AOPg)[which(is.na(V(AOPg)$AOP_ID))]
## Number of KEs without AOP IDs
length(V(AOPg)[which(is.na(V(AOPg)$AOP_ID))])

## Plot AOP WIKI emphasizing KEs without AOP_IDS
## These, upon a brief investigation, are those KEs that have been DELETED mostly, but remain in the KER table.
V(AOPg)$size<-1
V(AOPg)$size[which(is.na(V(AOPg)$AOP_ID))]<-3
jpeg.netplot(plot(AOPg,layout=layout.fruchterman.reingold(AOPg),vertex.size=1.5,vertex.label=NA, edge.arrow.size=.5),"AOPwiki")
set.seed(1)
plot(AOPg,layout=layout.fruchterman.reingold(AOPg),vertex.size=2,vertex.label=NA, edge.arrow.size=.08)
aplot(AOPg)
## Export KEs without AOPs Plot
# jpeg.netplot(plot(AOPg,layout=layout.fruchterman.reingold(AOPg),vertex.label=NA, edge.arrow.size=.08),"NO_AOP_IDs",seedval=1)

####~ Level of Biological Organization Viz ####

## Add level of biological organization for key event nodes
V(AOPg)$lobo<-KEdata[[4]][match(V(AOPg)$KE_name,KEdata[[2]])]
V(AOPg)$lobo[which(is.na(V(AOPg)$lobo))]<-"" #assigns blank to NA data
tcols=rainbow(length(unique(V(AOPg)$lobo))) #creates a color scheme for visualization
lobo_list=c("Molecular","Cellular","Tissue","Organ","Individual","Population","") #creates an ordering of biological organization
V(AOPg)$lobo_o<-match(V(AOPg)$lobo,lobo_list) #assigns a value of biological organization instead of string.  1=molecular, 2=cellular, ...
lobo_freqs<-table(V(AOPg)$lobo_o)

length(V(AOPg)$lobo[which(V(AOPg)$lobo=="Molecular")])

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
       col="#777777", xjust=1,yjust=1, pt.bg=tcols, pt.cex=2, cex=.8, bty="n", ncol=1, y.intersp=.75, box.col="white", text.col="black")

#### MIE TO AO ANALYSES ####
##~~~~~~~~~~~~~~~~~~~~~~~~##

####~~ Full Simple Path Analysis ####
## Combinations of MIE and AOs: If each MIE was connected to every possible AO with a single path
length(V(AOPg)[which(V(AOPg)$KE_KED=="MIE")])*length(V(AOPg)[which(V(AOPg)$KE_KED=="AO")])
#result: 8918

## All Linear AOPs in AOP-wiki
head(rev(sort(aop_paths)))
V(AOPg)$aop_size<-1
V(AOPg)[205]$aop_size<-4
aplot(AOPg,vsize=V(AOPg)$aop_size)

## Total Simple Paths between MIEs and AOs
aop.paths(AOPg,tot=T)
#result: 43252

## Total simple paths by MIE/AO pairs
mie2aopths<-aop.paths(AOPg,mie=T)
m2apths<-as.data.frame(mie2aopths)
m2apths[order(-m2apths$paths),]
head(sort(mie2aopths$paths,decr=T))

## calculate and plot normalized path counts
aop.paths(AOPg,normalized=T)
plot(sort(AOPg.aop.path.cts))
dev.off()
plot(m2apths$paths)

####~~~ Simple Path Subgraph ####

## Create a subgraph for the MIE AO pair with the most simple paths between them
sg<-induced_subgraph(AOPg,unique(unlist(all_simple_paths(AOPg, from=160, to=678, mode="out"))))

length(V(sg))
length(E(sg))

## Plot subgraph (using ASP colors and sizes)
jpeg("images/424_to_563.jpeg",hei=800,wid=800,qual=100)
set.seed(1)
plot(sg,vertex.size=8, edge.arrow.size=.5,
     vertex.label=V(sg)$name, vertex.color=V(sg)$ked_color,
     main="MIE/AO Subgraph for MIE EventID:424 to AO EventID:563") ##Notice, an error, since there is a cycle in the sg
dev.off()
#Verify strongly connected components (cycles)
components(sg,mode="strong") ## Two exist, both with 4 nodes in them

##Visualization for strongly connected components 

V(sg)$cc_color<-unlist(color.comps(sg,"strong")$vcol)  #color.comps is a custom function stored in the AOP_net_functions.R file
E(sg)$cc_color<-unlist(color.comps(sg,"strong")$ecol)  #color.comps is a custom function stored in the AOP_net_functions.R file
V(sg)$cc_color[which(V(sg)$KE_EID==563)]<-"red" ##colors AO
V(sg)$cc_color[which(V(sg)$KE_EID==424)]<-"green" #colors MIE 

V(sg)$cc_size<-unlist(color.comps(sg,"strong")$vsize)
E(sg)$cc_width<-unlist(color.comps(sg,"strong")$ewidth)

## Plot of connected components with strong sizing option
set.seed(1)
jpeg("images/424_to_563_SCC.jpeg",hei=800,wid=800,qual=100)
plot(sg,vertex.size=8, vertex.color=V(sg)$cc_color, edge.color=E(sg)$cc_color, 
     main="MIE/AO Subgraph for MIE EventID:424 to AO EventID:563 - Strongly Connected Components",
     edge.arrow.size=.75, edge.width=E(sg)$cc_width*.9,vertex.label=V(sg)$name)
dev.off()

f57t677<-induced_subgraph(AOPg,unique(unlist(all_simple_paths(AOPg, from=57, to=677, mode="out"))))


## Export plot for subgraph
jpeg.netplot(plot(f57t677, vertex.size=10, edge.width=E(f57t677)$asp_size, edge.color=E(f57t677)$asp_clr, edge.arrow.size=1, vertex.label=V(f57t677)$KE_EID, vertex.color=V(f57t677)$ked_color),
             "f57t677_sp",seedval=1,maii=c(0,0,0,0))

## Export topo plot for subgraph
jpeg.netplot(plot(f57t677,layout=topo.layout(f57t677),
                  vertex.label.degree=0, vertex.label.dist=1.18, edge.curved=1, vertex.label.color="Black",
                  vertex.size=7, edge.width=E(f57t677)$asp_size, edge.color=E(f57t677)$asp_clr, edge.arrow.size=1,
                  vertex.label=V(f57t677)$KE_EID, vertex.color=V(f57t677)$ked_color)
             ,"f57t677_topo_EID",seedval=1,maii=c(0,0,0,2.1))

### Shortest path analyses for subgraphs ###

## Need to find the vertex ID's of the nodes from AOPg within the
## subgraph created in order to conduct path analyses
V(AOPg)[57]$name
V(AOPg)[677]$name
which(V(f57t677)$KE_EID==V(AOPg)[57]$KE_EID)
which(V(f57t677)$KE_EID==V(AOPg)[677]$KE_EID)

## Provide coloring for unweighted shortest path
E(f57t677)$sp_cols<-"gray"
E(f57t677)$sp_cols[shortest.path.coloring(f57t677,f=4,t=17,weight=E(f57t677)$evidence)]<-"green"
E(f57t677)$sp_size<-1
E(f57t677)$sp_size[shortest.path.coloring(f57t677,f=4,t=17,weight=E(f57t677)$evidence)]<-2
set.seed(1)
plot(f57t677, vertex.size=10, edge.color=E(f57t677)$sp_cols, edge.arrow.size=.15, vertex.label=V(f57t677)$name, vertex.color=V(f57t677)$ked_color,edge.width=E(f57t677)$sp_size)

## Coloring for non-weighted shortest path
E(f57t677)$sp_cols<-"gray"
E(f57t677)$sp_cols[shortest.path.coloring(f57t677,f=4,t=17)]<-"orange"
E(f57t677)$sp_size<-1
E(f57t677)$sp_size[shortest.path.coloring(f57t677,f=4,t=17)]<-2

##Coloring for Quantitative Understanding shortest path
E(f57t677)$sp_cols<-"gray"
E(f57t677)$sp_cols[shortest.path.coloring(f57t677,f=4,t=17,weight=E(f57t677)$quant)]<-"blue"
E(f57t677)$sp_size<-1
E(f57t677)$sp_size[shortest.path.coloring(f57t677,f=4,t=17,weight=E(f57t677)$quant)]<-2

## Export topo.plot for shortest paths
## NOTE: Must evaluate the appropriate weight determination above
## to store as sp_size and sp_color attributes
jpeg.netplot(plot(f57t677,layout=topo.layout(f57t677),
vertex.label.degree=0, vertex.label.cex=1,vertex.label.dist=1.2, 
edge.curved=1, vertex.label.color="black", vertex.size=7, 
edge.color=E(f57t677)$sp_cols, edge.arrow.size=1, vertex.label=V(f57t677)$KE_name,
vertex.color=V(f57t677)$ked_color,edge.width=E(f57t677)$sp_size),
"f57t677_ev_topo",seedval=1,maii=c(0,0,0,2.1))

### LONGEST PATH ANALYSIS ###
##~~~~~~~~~~~~~~~~~~~~~~~~~##
## THere are some algorithms for finding the longest simple path, especially when the graph is DAG
## However in this case, it is being found by finding all simple paths and determining the max

asp<-all_simple_paths(f57t677,f=4,t=17,mode="out")

## Find the longest simple paths
aspl<-c()
for(i in 1:length(asp)){
  aspl[i]<-length(asp[[i]])}

alp<-asp[which(aspl==max(aspl))]

E(f57t677)$p_clrs<-"gray"
E(f57t677)$p_size<-1
#set color and size of edges in longest paths
for(i in 1:length(alp)){
  E(f57t677,path=alp[[i]],dir=T)$p_clrs<-"red"
  E(f57t677,path=alp[[i]],dir=T)$p_size<-2}


## Plot Simple Paths in Subgraph with topological layout
jpeg.netplot(plot(f57t677,layout=topo.layout(f57t677),
                  vertex.label.degree=0, vertex.label.cex=1,vertex.label.dist=1.2, 
                  edge.curved=1, vertex.label.color="black", vertex.size=7, 
                  edge.color=E(f57t677)$p_clrs, edge.arrow.size=1, vertex.label=V(f57t677)$KE_name,
                  vertex.color=V(f57t677)$ked_color,edge.width=E(f57t677)$p_size),
             "f57t677_avgqu_topo",seedval=1,maii=c(0,0,0,2.1))

asp<-all_simple_paths(f57t677,f=4,t=17,mode="out")

## KE NORMALIZED PATH DETECTION ##
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

## Find average evidence and average quantitative understanding in 
## all the simple paths
avgev<-c()
avgqu<-c()
for(i in 1:length(asp)){
  avgev[i]<-mean(E(f57t677,path=asp[[i]],dir=T)$evidence)
  avgqu[i]<-mean(E(f57t677,path=asp[[i]],dir=T)$quant)
}
min(avgev)
min(avgqu)
evopt<-asp[which(avgev==min(avgev))]
quopt<-asp[which(avgqu==min(avgqu))]
length(evopt)
length(quopt)


## Assign network colors and sizes based on optimal evidence paths
E(f57t677)$p_clrs<-"gray"
E(f57t677)$p_size<-1
for(i in 1:length(evopt)){
  E(f57t677,path=evopt[[i]],dir=T)$p_clrs<-"DarkGreen"
  E(f57t677,path=evopt[[i]],dir=T)$p_size<-"2"
}

## Assign network colors and sizes based on optimal quantitaitve understanding paths
E(f57t677)$p_clrs<-"gray"
E(f57t677)$p_size<-1
for(i in 1:length(quopt)){
  E(f57t677,path=quopt[[i]],dir=T)$p_clrs<-"DarkBlue"
  E(f57t677,path=quopt[[i]],dir=T)$p_size<-"2"
}


## Plot Simple Paths in Subgraph with topological layout
jpeg.netplot(plot(f57t677,layout=topo.layout(f57t677),
                  vertex.label.degree=0, vertex.label.cex=1,vertex.label.dist=1.2, 
                  edge.curved=1, vertex.label.color="black", vertex.size=7, 
                  edge.color=E(f57t677)$p_clrs, edge.arrow.size=1, vertex.label=V(f57t677)$KE_name,
                  vertex.color=V(f57t677)$ked_color,edge.width=E(f57t677)$p_size),
             "f57t677_avgqu_topo",seedval=1,maii=c(0,0,0,2.1))


####~~~ Simple Path Visualization####

## This part of the code can be used to visualize in the graph, all 
## simple paths between a specified node-node pair.  

## Specify from and to nodes by vertex number (not name/KE_ID)
from<-57
to<-677

## Assign color to simple paths
E(AOPg)$asp_clr<-"gray"# Set color for all edges
E(AOPg)$asp_clr[simple.path.coloring(AOPg,from,to)]<-"purple" #set color for edges in simple paths

## Assign sizes to simple paths
E(AOPg)$asp_size<-1 #set size for all edges
E(AOPg)$asp_size[simple.path.sizing(AOPg,from,to)]<-2 #set size for edges in simple paths

## Plot full network with simple paths between MIE and AO highlighted
set.seed(1)
plot(AOPg ,vertex.size=2, edge.width=E(AOPg)$asp_size, edge.color=E(AOPg)$asp_clr, edge.arrow.size=.1, vertex.label=NA, vertex.color=V(AOPg)$ked_color)
## Export network plot
jpeg.netplot(plot(AOPg ,vertex.size=2, edge.width=E(AOPg)$asp_size, edge.color=E(AOPg)$asp_clr, edge.arrow.size=.3, vertex.label=NA, vertex.color=V(AOPg)$ked_color),paste(from," to ",to,sep=""),seedval=1)


## This section was used to highlight the 5 MIE to AO pairs with the higest edge connectivity.
# E(AOPg)$asp_clr<-"gray"
# E(AOPg)$asp_clr[simple.path.coloring(AOPg,233,535)]<-"yellow"
# ## Assign sizes to simple paths
# E(AOPg)$asp_size<-1
# E(AOPg)$asp_size[simple.path.sizing(AOPg,233,535)]<-2
# ## Plot full network with simple paths between MIE and AO highlighted
# set.seed(1)
# plot(AOPg ,vertex.size=2, edge.width=E(AOPg)$asp_size, edge.color=E(AOPg)$asp_clr, edge.arrow.size=.1, vertex.label=NA, vertex.color=V(AOPg)$ked_color)
# ## Export network plot
# jpeg.netplot(plot(AOPg ,vertex.size=2, edge.width=E(AOPg)$asp_size, edge.color=E(AOPg)$asp_clr, edge.arrow.size=.3, vertex.label=NA, vertex.color=V(AOPg)$ked_color)
#              ,"f233t535",seedval=1)
# 
# E(AOPg)$asp_clr<-"gray"
# E(AOPg)$asp_clr[simple.path.coloring(AOPg,5,677)]<-"blue"
# ## Assign sizes to simple paths
# E(AOPg)$asp_size<-1
# E(AOPg)$asp_size[simple.path.sizing(AOPg,5,677)]<-2
# ## Plot full network with simple paths between MIE and AO highlighted
# set.seed(1)
# plot(AOPg ,vertex.size=2, edge.width=E(AOPg)$asp_size, edge.color=E(AOPg)$asp_clr, edge.arrow.size=.1, vertex.label=NA, vertex.color=V(AOPg)$ked_color)
# ## Export network plot
# 
# E(AOPg)$asp_clr[simple.path.coloring(AOPg,128,670)]<-"orange"
# E(AOPg)$asp_size[simple.path.sizing(AOPg,128,670)]<-2
# E(AOPg)$asp_clr[simple.path.coloring(AOPg,3,63)]<-"cyan"
# E(AOPg)$asp_size[simple.path.sizing(AOPg,3,63)]<-2
# 
# 
# jpeg.netplot(plot(AOPg ,vertex.size=2, edge.width=E(AOPg)$asp_size, edge.color=E(AOPg)$asp_clr, edge.arrow.size=.3, vertex.label=NA, vertex.color=V(AOPg)$ked_color)
#              ,"AOP_ecl_top5",seedval=1)

####~~ Within AOP Simple Path Analysis ####
### A look at the MIEs and AOs specified in the network and their relationships
### to AOP ids to determine within AOP simple path potentials
mies<-V(AOPg)[which(V(AOPg)$KE_KED=="MIE")]
length(mies)
aos<-V(AOPg)[which(V(AOPg)$KE_KED=="AO")]
length(aos)

## This identifies the set of AOP_IDs that have both an MIE and AO in them 
isecsall<-c()
for(i in 1:length(mies)){
  for(j in 1:length(aos)){
    isecs<-intersect(unlist(mies[i]$AOP_ID),unlist(aos[j]$AOP_ID))
    isecsall<-append(isecsall,isecs)}}

## View Result
sort(unique(isecsall))

## The number of AOPs have the potential for an MIE to AO path, however
## this does not guarantee that there are KEs connecting them.  
length(unique(isecsall))

## All simple paths within a specified AOP (by AOP ID)
## example is AOP_ID=1
in.aop.paths(AOPg,1)

## All simple paths within the entire AOP wiki.  This function loops
## over all unique AOP IDs within the AOP wiki and sums the number of paths
in.aop.paths.all(AOPg)
#result: 151

## This option of the in.aop.paths.all() function provides the paths instead
## of the total number of paths. However this paths are in 
#pathlist<-list(in.aop.paths.all(AOPg,paths=T))

####~~ Simple Path N/A AOP_ID REMOVED ####
## This analyis is being conducted in order to estimate the influence of KEs
## in the wiki that do not have AOP IDs.  These types of KEs can arise in multiple ways
## including test KERs or AOPs that have since been deleted or simply KEs that don't have
## a proper AOPID encoded.

## Create new graph object for analyses
AOPn<-induced.subgraph(AOPg,V(AOPg)[which(!is.na(V(AOPg)$AOP_ID))])

V(AOPg)$KE_KED[which(is.na(V(AOPg)$AOP_ID))]
## Plot to visualize
jpeg.netplot(aplot(AOPn),"AOPWiki_no_aop_id_nas")

## Combinations of MIE and AOs: If each MIE was connected to every possible AO with a single path
length(V(AOPn)[which(V(AOPn)$KE_KED=="MIE")])*length(V(AOPn)[which(V(AOPn)$KE_KED=="AO")])
#result: 8178

in.aop.paths.all(AOPn)
#result: 151

aop.paths(AOPn,tot=T)
#result: 32562
####~~ AOPwiki Growth Calculations ####

## Construct sequential AOPWiki by AOP ID and find MIE to AO paths in it.
aopids<-sort(unique(unlist(V()$AOP_ID)))
AOPg
## Create and run aop.paths analysis on sequential AOPwiki subgraphs
# vlist<-list()
# aopps<-list()
# isgs<-list()
# for(i in 1:length(unique(unlist(V(AOPg)$AOP_ID)))){
#   vlist[[i+1]]<-append(vlist[i],as.numeric(V(AOPg)[which(unlist(lapply(V(AOPg)$AOP_ID, function(x) is.element(aopids[i],x))))]))
#   kes<-unique(unlist(vlist))
#   isgs[[i+1]]<-induced.subgraph(AOPg,kes)
#   aopps[[i]]<-aop.paths(isg,tot=T)
# }

## Import aop.paths data from previous run 
aoppths<-read.table(paste(workingDir,"data/mie2aopaths.csv", sep=""), sep=",", stringsAsFactors=FALSE, header=F)
AOPno<-as.vector(unlist(aoppths))

## Create subgraph objects without aop.paths analysis
vlists<-list() #sequential list of vertices to create subgraphs
isgs<-list() #subgraph outputs
for(i in 1:length(unique(unlist(V(AOPg)$AOP_ID)))){
  if(i==1){vlists[[i]]<-as.numeric(V(AOPg)[which(unlist(lapply(V(AOPg)$AOP_ID, function(x) is.element(aopids[i],x))))])
  kes<-unique(unlist(vlists))
  isgs[[i]]<-induced.subgraph(AOPg,kes)}
  else{
  vlists[[i]]<-append(vlists[i-1],as.numeric(V(AOPg)[which(unlist(lapply(V(AOPg)$AOP_ID, function(x) is.element(aopids[i],x))))]))
  kes<-unique(unlist(vlists))
  isgs[[i]]<-induced.subgraph(AOPg,kes)
  }}


##Create and analyze AOPwiki subgraphs KE inclusion
isgs<-list()
for(i in 1:length(V(AOPg))){
  isgs[[i]]<-induced.subgraph(AOPg,V(AOPg)[1:i])
}

aops<-c()
st=748
for(i in st:length(isgs)){
  aops[i]<-aop.paths(isgs[[i]],tot=T)
  write.csv(aops,file=paste("results/aoppthcts",st,".csv",sep=""))
}

## Components Analysis for each subgraph
cmpnost<-c()
cmpnowk<-c()
cmpszwk<-c()
for(i in 1:length(isgs)){
  cmpnost[i]<-length(which(components(isgs[[i]],mode="strong")$csize>1))
  cmpnowk[i]<-components(isgs[[i]])$no
  cmpszwk[i]<-max(components(isgs[[i]])$csize)
}

## Number of KES and KERs in each subgraph
KEno<-c()
KERno<-c()
MIEno<-c()
AOno<-c()
MIEinMCno<-c()
AOinMCno<-c()
for(i in 1:length(isgs)){
  KEno[i]<-length(V(isgs[[i]]))
  KERno[i]<-length(E(isgs[[i]]))
  MIEno[i]<-length(which(V(isgs[[i]])$KE_KED=="MIE"))
  AOno[i]<-length(which(V(isgs[[i]])$KE_KED=="AO"))
  MIEinMCno[i]<-length(which(V(isgs[[i]])[which(which(components(isgs[[i]])$csize==max(components(isgs[[i]])$csize))==components(isgs[[i]])$membership)]$KE_KED=="MIE"))
  AOinMCno[i]<-length(which(V(isgs[[i]])[which(which(components(isgs[[i]])$csize==max(components(isgs[[i]])$csize))==components(isgs[[i]])$membership)]$KE_KED=="AO"))
  }
aoppths<-read.table(paste(workingDir,"data/aoppthcts.csv", sep=""), sep=",", stringsAsFactors=FALSE, header=F)
AOPno<-as.vector(unlist(aoppths))


##summary data frame
sgsum<-as.data.frame(cbind(no=1:length(isgs),keno=KEno,kerno=KERno,mieno=MIEno,aono=AOno,cns=cmpnost,cnw=cmpnowk,mcsw=cmpszwk))

####~~~ Visualization of AOPwiki Growth by Attribute####

#KE vs NO AOP LM
kelm<-lm(keno ~ no, data = sgsum)
summary(kelm)
#KEno
jpeg(file = paste("images/aopkeno.jpeg",sep=""),
     width=800, height=800, bg =bg_col, quality=100)
par(mar = c(5,5,2,5))
plot(AOPno,pch=2,col="red",ylab="Number of MIE to AO Paths",xlab="Number of User-Specified AOPs in AOPwiki Network")
par(new=T)
plot(KEno,pch=1,col="blue",axes=F,xlab=NA,ylab=NA,cex=1.2)
axis(side=4)
mtext(side=4,line=3,"Number of KEs")
#abline(36.11,4.06,col="blue")
legend(x=c(0),y=c(max(KEno)),legend=c("# AOPs","# KEs"),bty="n",pch=c(2,1),col=c("red","blue"))
dev.off()
#KER vs NO AOP LM
kerlm<-lm(kerno ~ no, data = sgsum)
summary(kerlm)
#KERno
jpeg(file = paste("images/aopkerno.jpeg",sep=""),
     width=800, height=800, bg =bg_col, quality=100)
par(mar = c(5,5,2,5))
plot(AOPno,pch=2,col="red",ylab="Number of MIE to AO Paths",xlab="Number of User-Specified AOPs in AOPwiki Network")
par(new=T)
plot(KERno,pch=1,col="purple",axes=F,xlab=NA,ylab=NA,cex=1.2)
axis(side=4)
mtext(side=4,line=3,"Number of KERs")
#abline(44.52,5.65,col="blue")
legend(x=c(0),y=c(max(KERno)),legend=c("# AOPs","# KERs"),bty="n",pch=c(2,1),col=c("red","purple"))
dev.off()

#KER/KE ratio
jpeg(file = paste("images/aopkekerrat.jpeg",sep=""),
     width=800, height=800, bg =bg_col, quality=100)
par(mar = c(5,5,2,5))
plot(AOPno,pch=2,col="red",ylab="Number of MIE to AO Paths",xlab="Number of User-Specified AOPs in AOPwiki Network")
par(new=T)
plot(KERno/KEno,pch=1,col="LightBlue",axes=F,xlab=NA,ylab=NA,cex=1.2)
axis(side=4)
mtext(side=4,line=3,"Ratio of KERs/KEs")
#abline(44.52,5.65,col="blue")
legend(x=c(0),y=c(max(KERno/KEno)),legend=c("# AOPs","# KERs/KEs"),bty="n",pch=c(2,1),col=c("red","lightblue"))
dev.off()

#Diff KER/KE ratio
jpeg(file = paste("images/aopkekerrat_diff.jpeg",sep=""),
     width=800, height=800, bg =bg_col, quality=100)
par(mar = c(5,5,2,5))
plot(diff(AOPno),pch=2,col="red",type="h",ylab="Difference in Number of MIE to AO Paths",xlab="Number of User-Specified AOPs in AOPwiki Network")
par(new=T)
plot(abs(diff(KERno/KEno)),pch=1,col="LightGreen",axes=F,xlab=NA,ylab=NA,cex=1.2,type="h")
axis(side=4)
mtext(side=4,line=3,"Magnitude Difference in Ratio of KERs/KEs")
#abline(44.52,5.65,col="blue")
legend(x=c(0),y=c(max(KERno/KEno)),legend=c("Diff AOPs","|delta(KERs/KEs)|"),bty="n",pch=c(2,1),col=c("red","lightgreen"))
dev.off()

#max weakly connected components size
jpeg(file = paste("images/aopmaxccsz.jpeg",sep=""),
     width=800, height=800, bg =bg_col, quality=100)
par(mar = c(5,5,2,5))
plot(AOPno,pch=2,col="red",ylab="Number of MIE to AO Paths",xlab="Number of User-Specified AOPs in AOPwiki Network")
par(new=T)
plot(cmpszwk,pch=9,col="green",axes=F,xlab=NA,ylab=NA,cex=1.2)
axis(side=4)
mtext(side=4,line=3,"Size of Maximum Weakly Connected Component")
#abline(44.52,5.65,col="blue")
legend(x=c(0),y=c(max(cmpszwk)),legend=c("# AOPs","Max Size, Weak CC"),bty="n",pch=c(2,9),col=c("red","green"))
dev.off()

#Number of weakly connected components
jpeg(file = paste("images/aopnoccwk.jpeg",sep=""),
     width=800, height=800, bg =bg_col, quality=100)
par(mar = c(5,5,2,5))
plot(AOPno,pch=2,col="red",ylab="Number of MIE to AO Paths",xlab="Number of User-Specified AOPs in AOPwiki Network")
par(new=T)
plot(cmpnowk,pch=9,col="orange",axes=F,xlab=NA,ylab=NA,cex=1.2)
axis(side=4)
mtext(side=4,line=3,"Number of Weakly Connected Component")
#abline(44.52,5.65,col="blue")
legend(x=c(0),y=c(max(cmpnowk)),legend=c("# AOPs","# Weak CC"),bty="n",pch=c(2,9),col=c("red","orange"))
dev.off()

#Number of strongly connected components
jpeg(file = paste("images/aopnoccst.jpeg",sep=""),
     width=800, height=800, bg =bg_col, quality=100)
par(mar = c(5,5,2,5))
plot(AOPno,pch=2,col="red",ylab="Number of MIE to AO Paths",xlab="Number of User-Specified AOPs in AOPwiki Network")
par(new=T)
plot(cmpnost,pch=10,col="LightGreen",axes=F,xlab=NA,ylab=NA,cex=1.2)
axis(side=4)
mtext(side=4,line=3,"Number of Strongly Connected Component")
#abline(44.52,5.65,col="blue")
legend(x=c(0),y=c(max(cmpnowk)),legend=c("# AOPs","# Strong CC"),bty="n",pch=c(2,10),col=c("red","LightGreen"))
dev.off()

#Number of strongly connected components
jpeg(file = paste("images/aopnoccst_diff.jpeg",sep=""),
     width=800, height=800, bg =bg_col, quality=100)
par(mar = c(5,5,2,5))
plot(AOPno,pch=2,col="red",ylab="Number of MIE to AO Paths",xlab="Number of User-Specified AOPs in AOPwiki Network")
par(new=T)
plot(cmpnost,pch=10,col="LightGreen",axes=F,xlab=NA,ylab=NA,cex=1.2)
axis(side=4)
mtext(side=4,line=3,"Number of Strongly Connected Component")
par(new=T)
plot(c(0,diff(cmpnost))*cmpnost,pch=10,col="LightGreen",axes=F,xlab=NA,ylab=NA,cex=1.2,type="h")
legend(x=c(0),y=c(max(cmpnost)),legend=c("# AOPs","# Strong CC"),bty="n",pch=c(2,10),col=c("red","LightGreen"))
dev.off()



#KEno & KERno
jpeg(file = paste("images/kenokerno.jpeg",sep=""),
     width=800, height=800, bg =bg_col, quality=100)
par(mar = c(5,5,2,5))
plot(KEno,pch=1,col="blue",ylab="Number of KEs",xlab="Number of User-Specified AOPs in AOPwiki Network")
par(new=T)
plot(KERno,pch=16,col="purple",axes=F,xlab=NA,ylab=NA,cex=1.2)
axis(side=4)
mtext(side=4,line=3,"Number of KERs")
#abline(44.52,5.65,col="blue")
legend(x=c(0),y=c(max(KERno)),legend=c("# KEs","# KERs"),bty="n",pch=c(1,16),col=c("blue","purple"))
dev.off()

####~~~ Visualization of AOPwiki KE Growth by Attribute####


#KE vs NO AOP LM
kelm<-lm(keno ~ no, data = sgsum)
summary(kelm)
#KEno
jpeg(file = paste("images/ke_aopkeno.jpeg",sep=""),
     width=800, height=800, bg =bg_col, quality=100)
par(mar = c(5,5,2,5))
plot(AOPno,pch=2,col="red",ylab="Number of MIE to AO Paths",xlab="Number of KEs in AOPwiki Network")
par(new=T)
plot(KEno,pch=1,col="blue",axes=F,xlab=NA,ylab=NA,cex=1.2)
axis(side=4)
mtext(side=4,line=3,"Number of KEs")
#abline(36.11,4.06,col="blue")
legend(x=c(0),y=c(max(KEno)),legend=c("# AOPs","# KEs"),bty="n",pch=c(2,1),col=c("red","blue"))
dev.off()
#KER vs NO AOP LM
kerlm<-lm(kerno ~ no, data = sgsum)
summary(kerlm)
#KERno
jpeg(file = paste("images/ke_aopkerno.jpeg",sep=""),
     width=800, height=800, bg =bg_col, quality=100)
par(mar = c(5,5,2,5))
plot(AOPno,pch=2,col="red",ylab="Number of MIE to AO Paths",xlab="Number of KEs in AOPwiki Network")
par(new=T)
plot(KERno,pch=1,col="purple",axes=F,xlab=NA,ylab=NA,cex=1.2)
axis(side=4)
mtext(side=4,line=3,"Number of KERs")
#abline(44.52,5.65,col="blue")
legend(x=c(0),y=c(max(KERno)),legend=c("# AOPs","# KERs"),bty="n",pch=c(2,1),col=c("red","purple"))
dev.off()

#KER/KE ratio
jpeg(file = paste("images/aopkekerrat.jpeg",sep=""),
     width=800, height=800, bg =bg_col, quality=100)
par(mar = c(5,5,2,5))
plot(AOPno,pch=2,col="red",ylab="Number of MIE to AO Paths",xlab="Number of KEs in AOPwiki Network")
par(new=T)
plot(KERno/KEno,pch=1,col="LightBlue",axes=F,xlab=NA,ylab=NA,cex=1.2)
axis(side=4)
mtext(side=4,line=3,"Ratio of KERs/KEs")
#abline(44.52,5.65,col="blue")
legend(x=c(0),y=c(max(KERno/KEno)),legend=c("# AOPs","# KERs/KEs"),bty="n",pch=c(2,1),col=c("red","lightblue"))
dev.off()

#Diff KER/KE ratio
jpeg(file = paste("images/aopkekerrat_diff.jpeg",sep=""),
     width=800, height=800, bg =bg_col, quality=100)
par(mar = c(5,5,2,5))
plot(diff(AOPno),pch=2,col="red",type="h",ylab="Difference in Number of MIE to AO Paths",xlab="Number of User-Specified AOPs in AOPwiki Network")
par(new=T)
plot(abs(diff(KERno/KEno)),pch=1,col="LightGreen",axes=F,xlab=NA,ylab=NA,cex=1.2,type="h")
axis(side=4)
mtext(side=4,line=3,"Magnitude Difference in Ratio of KERs/KEs")
#abline(44.52,5.65,col="blue")
legend(x=c(0),y=c(max(KERno/KEno)),legend=c("Diff AOPs","|delta(KERs/KEs)|"),bty="n",pch=c(2,1),col=c("red","lightgreen"))
dev.off()

#max weakly connected components size
jpeg(file = paste("images/ke_aopmaxccsz.jpeg",sep=""),
     width=800, height=800, bg =bg_col, quality=100)
par(mar = c(5,5,2,5))
plot(AOPno,pch=2,col="red",ylab="Number of MIE to AO Paths",xlab="Number of KEs in AOPwiki Network")
par(new=T)
plot(cmpszwk,pch=1,col="green",axes=F,xlab=NA,ylab=NA,cex=1.2)
axis(side=4)
mtext(side=4,line=3,"Size of Maximum Weakly Connected Component")
#abline(44.52,5.65,col="blue")
legend(x=c(0),y=c(max(cmpszwk)),legend=c("# AOPs","Max Size, Weak CC"),bty="n",pch=c(2,1),col=c("red","green"))
dev.off()

#Number of weakly connected components
jpeg(file = paste("images/ke_aopnoccwk.jpeg",sep=""),
     width=800, height=800, bg =bg_col, quality=100)
par(mar = c(5,5,2,5))
plot(AOPno,pch=2,col="red",ylab="Number of MIE to AO Paths",xlab="Number of KEs in AOPwiki Network")
par(new=T)
plot(cmpnowk,pch=1,col="orange",axes=F,xlab=NA,ylab=NA,cex=1.2)
axis(side=4)
mtext(side=4,line=3,"Number of Weakly Connected Component")
#abline(44.52,5.65,col="blue")
legend(x=c(0),y=c(max(cmpnowk)),legend=c("# AOPs","# Weak CC"),bty="n",pch=c(2,1),col=c("red","orange"))
dev.off()

#Number of strongly connected components
jpeg(file = paste("images/ke_aopnoccst.jpeg",sep=""),
     width=800, height=800, bg =bg_col, quality=100)
par(mar = c(5,5,2,5))
plot(AOPno,pch=2,col="red",ylab="Number of MIE to AO Paths",xlab="Number of KEs in AOPwiki Network")
par(new=T)
plot(cmpnost,pch=10,col="LightGreen",axes=F,xlab=NA,ylab=NA,cex=1.2)
axis(side=4)
mtext(side=4,line=3,"Number of Strongly Connected Component")
#abline(44.52,5.65,col="blue")
legend(x=c(0),y=c(max(cmpnowk)),legend=c("# AOPs","# Strong CC"),bty="n",pch=c(2,10),col=c("red","LightGreen"))
dev.off()

#Number of strongly connected components
jpeg(file = paste("images/ke_aopnoccst_diff.jpeg",sep=""),
     width=800, height=800, bg =bg_col, quality=100)
par(mar = c(5,5,2,5))
plot(AOPno,pch=2,col="red",ylab="Number of MIE to AO Paths",xlab="Number of KEs in AOPwiki Network")
par(new=T)
plot(cmpnost,pch=10,col="LightGreen",axes=F,xlab=NA,ylab=NA,cex=1.2)
axis(side=4)
mtext(side=4,line=3,"Number of Strongly Connected Component")
par(new=T)
plot(c(0,diff(cmpnost))*cmpnost,pch=10,col="LightGreen",axes=F,xlab=NA,ylab=NA,cex=1.2,type="h")
legend(x=c(0),y=c(max(cmpnost)),legend=c("# AOPs","# Strong CC"),bty="n",pch=c(2,10),col=c("red","LightGreen"))
dev.off()

#KEno & KERno
jpeg(file = paste("images/ke_kenokerno.jpeg",sep=""),
     width=800, height=800, bg =bg_col, quality=100)
par(mar = c(5,5,2,5))
plot(KEno,pch=1,col="blue",ylab="Number of KEs",xlab="Number of KEs in AOPwiki Network")
par(new=T)
plot(KERno,pch=16,col="purple",axes=F,xlab=NA,ylab=NA,cex=1.2)
axis(side=4)
mtext(side=4,line=3,"Number of KERs")
#abline(44.52,5.65,col="blue")
legend(x=c(0),y=c(max(KERno)),legend=c("# KEs","# KERs"),bty="n",pch=c(1,16),col=c("blue","purple"))
dev.off()

####~~~ Animation of AOPwiki Growth ####
library(magick)
library(dplyr)
library(purrr)
##Formatting outputs
library(stringr) #supplies a function to add leading zeros to strings (used in the plotting commands)
library(latex2exp) #Allows for TeX commands in plots 
##Library for 3d plotting
library(plot3D)

##Constructive plot of AOPwiki network 
#define layout for network
set.seed(1)
fulllo<-layout.fruchterman.reingold(AOPg)
## Loop for creating successive plots
for(i in 1:153){
  V(isgs[[153]])$color<-NA
  V(isgs[[153]])[!is.na(match(V(isgs[[153]])$KE_EID,V(isgs[[i]])$KE_EID))]$color<-'blue'
  E(isgs[[153]])$color<-NA
  E(isgs[[153]])[match(as_ids(E(isgs[[i]])),as_ids(E(isgs[[153]])))]$color<-'white'
  png(paste("images/networks/net_",str_pad(i, width=3, side="left", pad="0"),".png",sep=""),wid=800,hei=800)
  par(col.main="white",bg="black")
  plot(isgs[[153]],margin=c(-.05,-.05,-.05,-.05),layout=fulllo,vertex.label=NA, vertex.frame.color=NA, edge.arrow.size=.5,vertex.size=2,main="AOPwiki Network Evolution")
text(0,y=-1.1,paste("-- Network Summary -- \n Key Events:",KEno[i]," Key Event Relationships:",KERno[i]," Linear AOPs:",AOPno[i]),font=2,col="white")
dev.off()
}


list.files(path=paste(getwd(),"/images/networks",sep=""),pattern="*.png",full.names=T) %>%
  map(image_read) %>%
  image_join() %>%
  image_animate(fps=4,loop=1) %>%
  image_write("images/net_evolution.gif")

## Plot of all AOPwiki sub-networks, does not maintain layout...previous plots do maintain layout
for(i in 1:153){
  png(paste("images/networks/net_",str_pad(i, width=3, side="left", pad="0"),".png",sep=""))
  aplot(isgs[[i+1]],vsize=2)
  dev.off()
}
list.files(path=paste(getwd(),"/images/networks",sep=""),pattern="*.png",full.names=T) %>%
  map(image_read) %>%
  image_join() %>%
  image_animate(fps=4) %>%
  image_write("images/nets.gif")

####~~ Edge Connectivity ####
# This extends the graph wide edge connectivity calculation to
# only consider paths between MIE and AO as designated by V(graph)$KE_KED

edge.connectivity(AOPg,source=V(AOPg)[which(V(AOPg)$name==408)],target=V(AOPg)[which(V(AOPg)$name==406)])

ecl<-aop.edge.connectivity(AOPg,kelist=V(AOPg)$KE_KED)

ecl_sort<-ecl[order(ecl[,3],decreasing=T),]
colnames(ecl_sort)<-c("MIE","AO","ECL")

ecl_short<-head(ecl_sort,10)

ecl_short_name<-cbind(ecl_short[,1],V(AOPg)$KE_name[ecl_short[,1]],ecl_short[,2],V(AOPg)$KE_name[ecl_short[,2]],ecl_short[,3])
colnames(ecl_short_name)<-c("MIE (Event ID)","MIE (Short Name)","AO (Event ID)","AO (Short Name)","ECL Value")
#name file to output to and create R variable for it
output<-file("results/ecl_results.txt")
#tell what to write into "output" file variable
write.table(ecl_short_name, output,col.names=T, sep="\t", row.names=F)
ecl_short
ecl_short_name

####~~ AOP Occurrence NEEDS REFINEMENT ####
V(AOPg)$aopp<-aop.paths(AOPg,kelist=V(AOPg)$KE_KED)
heatgrad=rev(heat.colors(n=max(V(AOPg)$aopp)))
heatgrad<-append(heatgrad, "#FFFFFF", after=0) #in the instance that a KE is not included in paths between MIE and AO...

head(sort(V(AOPg)$aopp,decr=T))
which(V(AOPg)$aopp==max(V(AOPg)$aopp))
V(AOPg)[205]$KE_name
mean(V(AOPg)$aopp)
plot(table(V(AOPg)$aopp))
median(V(AOPg)$aopp)
wbpal=colorRampPalette(c("white","blue"))(n=max(V(AOPg)$aopp)+1)

V(AOPg)$po_col<-wbpal[V(AOPg)$aopp+1]
V(AOPg)$po_size<-(V(AOPg)$aopp/max(V(AOPg)$aopp))*7
set.seed(1)
jpeg("images/AOPwiki_po.jpeg",height=800,wid=800,qual=100)
set.seed(1)
plot(AOPg ,vertex.size=V(AOPg)$po_size, edge.color="gray", edge.arrow.size=.05, vertex.label=NA, vertex.color=V(AOPg)$po_col)
dev.off()
aop.paths(AOPg,kelist=V(AOPg)$KE_KED)

#### CONNECTED COMPONENTS ANALYSIS####

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

wcomps
####~~ Strongly connected component analysis and plotting ####

comps<-components(AOPg, mode="strong")
## count number of KEs involved in SCCs
cyc<-which(comps$csize>1)
length(na.exclude(match(comps$membership,cyc)))

## When the "strong" option is passed to color.comps, vsize and ewidth are calculated and can be used within plot
V(AOPg)$cc_color<-unlist(color.comps(AOPg,"strong")$vcol)  #color.comps is a custom function stored in the AOP_net_functions.R file
E(AOPg)$cc_color<-unlist(color.comps(AOPg,"strong")$ecol)  #color.comps is a custom function stored in the AOP_net_functions.R file

V(AOPg)$cc_size<-unlist(color.comps(AOPg,"strong")$vsize)
E(AOPg)$cc_width<-unlist(color.comps(AOPg,"strong")$ewidth)

## Plot of connected components with strong sizing option
par(bg="white")
set.seed(1)
plot(AOPg,vertex.size=V(AOPg)$cc_size, edge.width=E(AOPg)$cc_width, vertex.color=V(AOPg)$cc_color, edge.color=E(AOPg)$cc_color,  vertex.size=2, edge.arrow.size=.35, vertex.label=NA)

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

head(rev(sort(unlist(V(AOPg)$AOP_ID))))

#### CENTRALITY MEASURES ####
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~##

####~~ Degree centrality ####


####~~~ Total Degree ####

## Which key events have the most incident nodes?
sort(degree(AOPg, mode="all"),decr=T)

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
head(sort(degree(AOPg,mode="in"),decr=T))
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

####~~~ In-path Closeness ####

clsmode="in"
## Prints the top-ten key event names by closeness value
cls<-closeness(AOPg,mode=clsmode)
sort(cls)[1]
which(cls==min(cls))
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

####~~~ Out-path Closeness ####

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
sort(eccentricity(AOPg,mode =eccmode),decr=T)
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


#### ADDITIONAL ANALYSES ####
##~~~~~~~~~~~~~~~~~~~~~~~~~##
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
