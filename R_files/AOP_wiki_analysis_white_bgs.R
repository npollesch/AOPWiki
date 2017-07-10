#### Visualize Adverse Outcome Pathway (AOP) WIKI Data ####
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
##  Set working directory and import key event relationships
library(igraph)
par(bg="white")
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

length(E(AOPg))

# This identifies which KEs are included in KERs, but are not themselves included in the KE event listings.
# V(AOPg)$KE_name[which(is.na(V(AOPg)$KE_EID))]

# Plot the AOP wiki, colored by AOP
#acols=topo.colors(length(unique(V(AOPg)$AOP_ID)))
acols=colorRampPalette(c("green","red","cyan","orange","magenta","yellow","blue"))

for(i in 1:length(unique(V(AOPg)$AOP_ID))){
  V(AOPg)[which(V(AOPg)$AOP_ID==unique(V(AOPg)$AOP_ID)[i])]$acol<-acols(length(unique(V(AOPg)$AOP_ID)))[i]
  }
par(bg="white",xpd=FALSE)
set.seed(1)
plot(AOPg,layout=layout.fruchterman.reingold(AOPg),  vertex.color=V(AOPg)$acol,vertex.label=NA, vertex.size=2, edge.arrow.size=.1)
#Calculates number of KE per unique AOP ID
AOP_freqs<-table(V(AOPg)$AOP_ID)
#Histogram of number of KE per unique AOP ID
hist(AOP_freqs,col.axis="white",xlab="# Key Events",ylab="Frequency",col.lab="white",col="white")
hist(AOP_freqs,col.axis="black",xlab="Key Events per AOP",ylab="Frequency",col.lab="black",col="gray")
#Barplot of number of KE per unique AOP ID with red line to show mean
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
# vertex.size=V(AOPg)$cc_size, edge.width=E(AOPg)$cc_width,

# plot of connected components
par(bg="white")
set.seed(1)
plot(AOPg,vertex.size=V(AOPg)$cc_size,vertex.color=V(AOPg)$cc_color, edge.color=E(AOPg)$cc_color,  vertex.size=2, edge.arrow.size=.1, vertex.label=NA)

## barplot for size of weakly connected components
wcomps<-components(AOPg,mode="weak")
wcc_freqs<-table(wcomps$csize)
bp_wcc<-barplot(table(wcomps$csize),col.axis="black", xlab="Component size",ylab="Frequency",col.lab="black")

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
lobo_list=c("Molecular","Cellular","Tissue","Organ","Individual","Population","") #creates an ordering of biological organization
V(AOPg)$lobo_o<-match(V(AOPg)$lobo,lobo_list) #assigns a value of biological organization instead of string.  1=molecular, 2=cellular, ...
lobo_freqs<-table(V(AOPg)$lobo_o)

# a plot the AOP wiki using a standard left to right lobo layout.
plot(AOPg, layout=lobo.layout(AOPg),vertex.size=2,  edge.curved=.3, edge.color="gray", edge.arrow.size=.1, vertex.label=NA, vertex.color=V(AOPg)$lobo_col)
#legend('topright',c("Molecular","Cellular","Tissue","Organ","Individual","Population","Not Specified"), pch=22,
#       col="#777777", xjust=1,yjust=1, pt.bg=tcols, pt.cex=2, cex=.8, bty="n", ncol=1, y.intersp=.5, box.col="black", text.col="white")


## Barplot of lobo frequency
par(bg="white")
xx<- barplot(table(V(AOPg)$lobo_o), ylab="# of Key Events", col=tcols, axes=F,names.arg=NA)
text(x=xx, y=10, label=lobo_freqs, cex=.75)
legend('topright',c("Molecular","Cellular","Tissue","Organ","Individual","Population","Not Specified"), pch=22,
       col="#777777", xjust=1,yjust=1, pt.bg=tcols, pt.cex=2, cex=.8, bty="n", ncol=1, y.intersp=.5, box.col="white", text.col="black")

#plot the AOP wiki by lobo info
V(AOPg)$lobo_col<-tcols[V(AOPg)$lobo_o]
set.seed(1)
plot(AOPg ,vertex.size=2, edge.color="gray", edge.arrow.size=.1, vertex.label=NA, vertex.color=V(AOPg)$lobo_col)

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
sort(degree(AOPg, mode="all"))

rev(as.integer(tail(sort(degree(AOPg, mode="all")),10)))
as.integer(names(tail(sort(degree(AOPg, mode="all")),10)))
 rev(V(AOPg)$KE_name[match(as.integer(names(tail(sort(degree(AOPg, mode="all")),10))),V(AOPg)$name)])
 
V(AOPg)$KE_name[which(V(AOPg)$name==449)]
  V(AOPg)$cent_size[which(V(AOPg)$name==345)]<-3
  V(AOPg)$cent_col[which(V(AOPg)$name==345)]<-"blue"

## DEGREE ALL VISUALIZATIONS ##

#global degree coloring for network plot
V(AOPg)$deg_col<-deg.col.grad(AOPg,dmode="all",c("black","green"))
#colored by degree
plot(AOPg, vertex.size=2, vertex.color=V(AOPg)$deg_col, edge.arrow.size=.1, vertex.label=NA, edge.color="gray",edge.width=1)
#colored by degree and sized by degree
set.seed(1)
plot(AOPg, vertex.size=500*degree(AOPg,mode="all",normalized=TRUE), vertex.color=V(AOPg)$deg_col, edge.arrow.size=.1, vertex.label=NA, edge.color="gray",edge.width=1)
#barplot to show degree histogram
bgpal=colorRampPalette(c("black","green"))
barplot(table(degree(AOPg,mode="all")), xlab="Degree (all)", ylab="# of KEs with Degree (all)",col.axis="white", col.lab="white",col=rev(heat.colors(max(degree(AOPg,mode="all"))+1)))
barplot(table(degree(AOPg,mode="all")), xlab="Total Degree", ylab="# of KEs with Total Degree",col.axis="black", col.lab="black",col=bgpal(max(degree(AOPg,mode="all"))))
legend('topright',title="Total Degree",title.adj=c(0.5,-1),legend=rev(seq(min(degree(AOPg,mode="all")),max(degree(AOPg,mode="all")),3)), pch=22,
       col="#777777", xjust=1,yjust=1, pt.bg=rev(bgpal(length(rev(seq(min(degree(AOPg,mode="all")),max(degree(AOPg,mode="all")),3))))), pt.cex=2, cex=.8, bty="n", ncol=1, y.intersp=.5, box.col="white", text.col="black")

max(degree(AOPg,mode="all"))
## DEGREE IN VISUALIZATIONS ##

#global degree coloring for network plot
sort(degree(AOPg,mode="in"))
rev(as.integer(tail(sort(degree(AOPg, mode="in")),10)))

rev(V(AOPg)$KE_name[match(as.integer(names(tail(sort(degree(AOPg, mode="in")),10))),V(AOPg)$name)])



V(AOPg)$KE_name[which(V(AOPg)$name==449)]


V(AOPg)$deg_col<-deg.col.grad(AOPg,dmode="in",c("black","purple"))
#colored by degree
plot(AOPg, vertex.size=2, vertex.color=V(AOPg)$deg_col, edge.arrow.size=.1, vertex.label=NA, edge.color="gray",edge.width=1)
#colored by degree and sized by degree
set.seed(1)
plot(AOPg, vertex.size=500*degree(AOPg,mode="in",normalized=TRUE), vertex.color=V(AOPg)$deg_col, edge.arrow.size=.1, vertex.label=NA, edge.color="gray",edge.width=1)
#barplot to show degree histogram
barplot(table(degree(AOPg,mode="in")), xlab="Degree (in)", ylab="# of KEs with Degree (in)",col.axis="white", col.lab="white",col=rev(heat.colors(max(degree(AOPg,mode="in"))+1)))
bgpal=colorRampPalette(c("black","purple"))
barplot(table(degree(AOPg,mode="in")), xlab="In-Degree", ylab="# of KEs with In-Degree",col.axis="black", col.lab="black",col=bgpal(max(degree(AOPg,mode="in"))))

legend('topright',title="In-Degree",title.adj=c(0.5,-1),legend=rev(seq(min(degree(AOPg,mode="in")),max(degree(AOPg,mode="in")),3)), pch=22,
       col="#777777", xjust=1,yjust=1, pt.bg=rev(bgpal(length(rev(seq(min(degree(AOPg,mode="in")),max(degree(AOPg,mode="in")),3))))), pt.cex=2, cex=.8, bty="n", ncol=1, y.intersp=.5, box.col="white", text.col="black")

## DEGREE OUT VISUALIZATIONS ##

#global degree coloring for network plot
sort(degree(AOPg,mode="out"))
rev(as.integer(tail(sort(degree(AOPg, mode="out")),10)))
rev(V(AOPg)$KE_name[match(as.integer(names(tail(sort(degree(AOPg, mode="out")),10))),V(AOPg)$name)])




V(AOPg)$KE_name[which(V(AOPg)$name==345)]

V(AOPg)$deg_col<-deg.col.grad(AOPg,dmode="out",c("black","Orange"))
#colored by degree
plot(AOPg, vertex.size=2, vertex.color=V(AOPg)$deg_col, edge.arrow.size=.1, vertex.label=NA, edge.color="gray",edge.width=1)
#colored by degree and sized by degree
set.seed(1)
plot(AOPg, vertex.size=500*degree(AOPg,mode="out",normalized=TRUE), vertex.color=V(AOPg)$deg_col, edge.arrow.size=.1, vertex.label=NA, edge.color="gray",edge.width=1)
#barplot to show degree histogram
barplot(table(degree(AOPg,mode="in")), xlab="Degree (out)", ylab="# of KEs with Degree (out)",col.axis="white", col.lab="white",col=rev(heat.colors(max(degree(AOPg,mode="out"))+1)))
bgpal=colorRampPalette(c("black","orange"))
barplot(table(degree(AOPg,mode="out")), xlab="Out-Degree", ylab="# of KEs with Out-Degree",col.axis="black", col.lab="black",col=bgpal(max(degree(AOPg,mode="out"))))

legend('topright',title="Out-Degree",title.adj=c(0.5,-1),legend=rev(seq(min(degree(AOPg,mode="out")),max(degree(AOPg,mode="out")),3)), pch=22,
       col="#777777", xjust=1,yjust=1, pt.bg=rev(bgpal(length(rev(seq(min(degree(AOPg,mode="out")),max(degree(AOPg,mode="out")),3))))), pt.cex=2, cex=.8, bty="n", ncol=1, y.intersp=.5, box.col="white", text.col="black")

#plot(degree.distribution(AOPg), col.axis="white", col.lab="white", pch=20,col="white", xlab="Degree", ylab="Percentage (%)")
#### Betweenness Centrality ####

# Which key event is involved in the most shortest paths between other key events?
sort(betweenness(AOPg,directed=TRUE))
table(betweenness(AOPg,directed=TRUE))

V(AOPg)$KE_name[which(V(AOPg)$name==345)]
V(AOPg)$cent_size[which(V(AOPg)$name==345)]<-3
V(AOPg)$cent_col[which(V(AOPg)$name==345)]<-"green"
V(AOPg)$deg_col<-deg.col(AOPg)
#colored by betweenness
primary.colors()

wbpal=colorRampPalette(c("black","blue"))
V(AOPg)$bet_col<-wbpal(20)[as.numeric(cut(betweenness(AOPg,directed=TRUE),breaks = 20))]
V(AOPg)$bet_col<-wbpal(20)[as.numeric(cut(betweenness(AOPg,directed=TRUE),breaks = 20))]
par(bg="white")
plot(AOPg, vertex.size=2, vertex.color=V(AOPg)$bet_col, edge.arrow.size=.1, vertex.label=NA, edge.color="gray",edge.width=1)
#colored by degree and sized by degree
set.seed(1)
plot(AOPg, vertex.size=3000*betweenness(AOPg,normalized=TRUE,directed=T), vertex.color=V(AOPg)$bet_col, edge.arrow.size=.1, vertex.label=NA, edge.color="gray",edge.width=1)
#barplot to show degree histogram
hist(log10(betweenness(AOPg)),breaks=20,col=wbpal(20), main=expression(paste("Histogram of Log"[10],"-Transformed Betweenness","")), xlab=expression(paste("Log"[10],"(Betweeness)","")))

barplot(table(as.numeric(cut(betweenness(AOPg),breaks = 20))),col=wbpal(20))

plot(sort(betweenness(AOPg)),col="blue", col.axis="white",col.lab="white", pch=10)

legend('topright',title="Eccentricity",title.adj=c(0.5,-1),legend=rev(seq(min(eccentricity(AOPg,mode ="all")),max(eccentricity(AOPg,mode ="all")),3)), pch=22,
       col="#777777", xjust=1,yjust=1, pt.bg=rev(bgpal(length(rev(seq(min(eccentricity(AOPg,mode ="all")),max(eccentricity(AOPg,mode ="all")),3))))), pt.cex=2, cex=.8, bty="n", ncol=1, y.intersp=.5, box.col="white", text.col="black")





#### Closeness Analysis ####
# Which key event is the closest to the rest?
sort(closeness(AOPg,mode="all"))
sort(closeness(AOPg,mode="all")*10^6)
V(AOPg)$KE_name[which(V(AOPg)$name==345)]
V(AOPg)$cent_size[which(V(AOPg)$name==711)]<-3
V(AOPg)$close_col[which(V(AOPg)$name==345)]<-"purple"

wbpal=colorRampPalette(c("black","magenta"))
V(AOPg)$close_col<-wbpal(1000)[as.numeric(cut(closeness(AOPg,mode="all",norm=TRUE),breaks = 1000))]

plot(closeness(AOPg,mode="all"), xlab="", xaxt='n', ylab="",main="KE Closeness in AOPwiki",col=V(AOPg)$close_col)

set.seed(1)
plot(AOPg, vertex.size=1000*closeness(AOPg,normalized=TRUE,mode="all"), vertex.color=V(AOPg)$close_col, edge.arrow.size=.1, vertex.label=NA)#, vertex.color="orange",edge.color="gray")

hist(closeness(AOPg,mode="all")*10^6,breaks=20,col=wbpal(20))



#### Eccentricity ####
# The eccentricity of a vertex is its shortest path distance from the farthest other node in the graph.
V(AOPg)$ecc<-eccentricity(AOPg,mode ="all")
sort(eccentricity(AOPg,mode ="all"))

#Finds the KE names for the three KEs with the largest eccentricity values
V(AOPg)$KE_name[which(V(AOPg)$name==35)]
V(AOPg)$KE_name[which(V(AOPg)$name==513)]
V(AOPg)$KE_name[which(V(AOPg)$name==590)]


wbpal=colorRampPalette(c("black","cyan"))
V(AOPg)$ecc_col<-wbpal(max(V(AOPg)$ecc)+1)[V(AOPg)$ecc+1]

set.seed(1)
plot(AOPg, vertex.size=4*(V(AOPg)$ecc/max(V(AOPg)$ecc)), vertex.color=V(AOPg)$ecc_col, edge.arrow.size=.1, vertex.label=NA, edge.color="gray",edge.width=1)

bgpal=colorRampPalette(c("black","cyan"))
barplot(table(eccentricity(AOPg,mode ="all")), xlab="Eccentricity", ylab="# of KEs with Eccentricity",col.axis="black", col.lab="black",col=bgpal(max(eccentricity(AOPg,mode ="all"))))

legend('topright',title="Eccentricity",title.adj=c(0.5,-1),legend=rev(seq(min(eccentricity(AOPg,mode ="all")),max(eccentricity(AOPg,mode ="all")),3)), pch=22,
       col="#777777", xjust=1,yjust=1, pt.bg=rev(bgpal(length(rev(seq(min(eccentricity(AOPg,mode ="all")),max(eccentricity(AOPg,mode ="all")),3))))), pt.cex=2, cex=.8, bty="n", ncol=1, y.intersp=.5, box.col="white", text.col="black")








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

