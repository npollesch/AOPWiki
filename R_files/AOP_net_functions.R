
# color.comps() is a function to color all the strongly connected components in a graph.  In theory it should do the same thing as cycle.color below
# however, I could not get cycle.color to work properly on the AOPwiki graph, so I made this one.  -Nate 5/11/17
color.comps<-function(gr, ccmode="strong"){ #function to color all non-trivial strongly connected components in a graph
  V(gr)$color<-"gray"
  E(gr)$color<-"gray"
  comps<-components(gr, mode=ccmode)
  ntcomps<-which(comps$csize>1) #non-trivial ccs (i.e. with more than 1 node)
  cols=rainbow(length(ntcomps))
  V(gr)$cc<-comps$membership
  for(i in 1:length(ntcomps)){
    V(gr)[which(V(gr)$cc==ntcomps[i])]$color<-cols[i]
    edgecombcc<-expand.grid(V(gr)[which(comps$membership==ntcomps[i])],V(gr)[which(comps$membership==ntcomps[i])]) #creates a pairwise list of all nodes in the cc
    edgecombflat<-as.vector(rbind(edgecombcc[[1]],edgecombcc[[2]])) #flattens the pairwise list to a vector where entries are read pairwise
    edges.in.cc<-get.edge.ids(gr,edgecombflat,directed=TRUE)
    E(gr)$color[edges.in.cc]<-cols[[i]]
  }
  if(ccmode=="strong"){
    V(gr)$size<-1 
    E(gr)$width<-1
    V(gr)$size[which(V(gr)$color!="gray")]<-3 # resizes all the highlighted nodes
    E(gr)$width[which(E(gr)$color!="gray")]<-3
    return(list(vcol=V(gr)$color,ecol=E(gr)$color,vsize=V(gr)$size,ewidth=E(gr)$width))}
  else {return(list(vcol=V(gr)$color,ecol=E(gr)$color))}
}


# cycle.color() function can be used to color SCC edges and vertices for visualization
# The output is a list of vertex colors[[1]] and node colors[[2]]

cycle.color = function(gr){
  gsccs<-components(gr, mode="strong")
  if(gsccs$no<length(V(gr))){ # is graph already acyclic?
    noccs<-which(gsccs$csize>1) # finds number of non-trivial sccs
    if(is.null(V(gr)$color)){V(gr)$color<-"Gray"}
    if(is.null(E(gr)$color)){E(gr)$color<-"Gray"}
    colors=list() #initializes colors list
    raincols=rainbow(length(noccs)+1) # 
    for(i in 1:length(noccs)){ # this loop goes through all sccs and assigns them the same node position within the vector map
      cc<-which(gsccs$csize>1)[[i]]
      V(gr)$color[which(gsccs$membership==cc)]<-raincols[[i+1]] #colors all nodes in cc with distinct rainbow color
      edgecombcc<-expand.grid(V(gr)[which(gsccs$membership==cc)],V(gr)[which(gsccs$membership==cc)]) #creates a pairwise list of all nodes in the cc
      edgecombflat<-as.vector(rbind(edgecombcc[[1]],edgecombcc[[2]])) #flattens the pairwise list to a vector where entries are read pairwise
      edges.in.cc<-get.edge.ids(gr,edgecombflat,directed=TRUE) #finds all edges between nodes included in pairwise vector
      E(gr)$color[edges.in.cc]<-raincols[[i+1]] #colors all edges in cc with distinct rainbow color corresponding to nodes
    }
    
    colors=list(vert.col=V(gr)$color,edge.col=E(gr)$color)
    return(colors)
  }
  
  else print("Graph is acyclic")
  
}


# condense.map() function determines condensation mapping needed to contract vertices in igraph
#o utput is a vector map to be used in 'condense.graph' function
condense.map = function(gr,suppress=FALSE){
  if(is.dag(gr)){print("Graph is acyclic")}
  else{gsccs<-components(gr, mode="strong")
  return(gsccs$membership)
  }
}

# condense.graph() function takes the condensation vector map and original graph and creates the condensed graph
# output is a new igraph object that has strongly connected components groups into supernodes 
condense.graph = function(gr,map){
  areColors <- function(x) { #function found at http://stackoverflow.com/questions/13289009/check-if-character-string-is-a-valid-color-representation by Josh O'Brien 
    sapply(x, function(X) {
      tryCatch(is.matrix(col2rgb(X)), 
               error = function(e) FALSE)
    })
  }
  
  cond.pre<-contract.vertices(gr,map,vertex.attr.comb=list(toString, label.cex="min")) 
  cond<-simplify(cond.pre, remove.multiple = TRUE, remove.loops = TRUE)
  V(cond)$color[which(areColors(V(cond)$color)==FALSE)]<-"Purple" # assigns condensed nodes the color purple
  a<-paste(deparse(substitute(gr)),"cond", sep = ".") # creates name for output from function input
  assign(a, cond, envir=globalenv()) # creates condensation map output
  print(paste("Condensed graph stored as:", a))
  
}

# topo.layout() function computes a topological ordering layout for visualizaiton of an igraph object
topo.layout = function(gr){
  tsort<-topo_sort(gr, mode = c("out"))
  tpos<-match(V(gr)$name,tsort$name) # finds position of node names in topoloigcal ordering
  tlay<-cbind(tpos,tpos) #define a plot layout that orders accoring to topological order
  a<-paste(deparse(substitute(gr)),"topo","lay", sep = ".") # creates name for output from function input
  x<-assign(a, tlay, envir=globalenv()) # creates layout global variable
  return(x) #makes the layout object the output
}

# lobo.layout() creates a layout in increasing order of level of biological organization
lobo.layout<-function(gr,lobos=lobo_list){ #assumes lobo data is stored as V(gr)$lobo 
  if(exists("lobos")) lobos else lobos=c("Molecular","Cellular","Tissue","Organ","Individual","Population","")
  
  for(i in 1:length(lobos)){ 
    lobo_len<-length(which(V(gr)$lobo==lobos[i]))
    V(gr)[which(V(gr)$lobo==lobos[i])]$lobo_locx<-rep(i,lobo_len)
    V(gr)[which(V(gr)$lobo==lobos[i])]$lobo_locy<-seq(1,lobo_len)
  }
  lobo_locs<-cbind(V(gr)$lobo_locx,V(gr)$lobo_locy)
  return(lobo_locs)
}

# aop.paths() function computes path counts for nodes involved in simple paths between MIEs and AOs.  
aop.paths= function(gr,normalized=FALSE,kelist = V(gr)$ked){ #kelist is a list of key event designation characters (MIE,KE, or AO) corresponding to nodes in the graph gr
  if(is.null(kelist)){print("Error: No key event designation list supplied")} 
  #is.character(V(sg.cond)$ked)
  else{
    path.counts<-data.frame(n=1:length(V(gr)),count=0) #creates a dataframe to store path counts
    # Creates node number list of MIEs to create a 'from list' for the paths
    path.counts$name<-V(gr)$name
    mie.list<-which(kelist=="MIE")
    ao.list<-which(kelist=="AO")
    for(fromnode in mie.list){
      for(tonode in ao.list){
        x<- all_simple_paths(gr,fromnode, to = tonode, mode="out")
        if(length(x)>0){
        #assign(paste("from",fromnode,"to",tonode,sep="."), x)  #Can uncomment to create lists of paths
        for(i in 1:length(x)){
          path.counts$count[as.integer(x[[i]])]<-path.counts$count[as.integer(x[[i]])]+1}}
      } 
    }
    if(normalized==FALSE){
      a<-paste(deparse(substitute(gr)),"aop","path","cts", sep = ".") # creates name for output from function input
      x<-assign(a, path.counts$count, envir=globalenv()) # creates layout global variable
      return(x) #makes the layout object the output
    }
    else{
      a<-paste(deparse(substitute(gr)),"aop","path","cts", sep = ".") # creates name for output from function input
      x<-assign(a, path.counts$count/max(path.counts$count), envir=globalenv()) # creates layout global variable
      return(x)
    }
    #set.vertex.attribute(gr,"aopp",path.counts$count)
    #set.vertex.attribute(gr,"aoppn",path.counts$count/max(path.counts$count))
  }
}

aop.edge.connectivity= function(gr,kelist = V(gr)$ked,names=F){ #kelist is a list of key event designation characters (MIE,KE, or AO) corresponding to nodes in the graph gr
  if(is.null(kelist)){print("Error: No key event designation list supplied")} 
  #is.character(V(sg.cond)$ked)
  else{
    mie.list<-which(kelist=="MIE")
    ao.list<-which(kelist=="AO")
     ec.list<-list()
     ec.listn<-list()
    i=0;
    for(fromnode in mie.list){
      for(tonode in ao.list){
        i=i+1
        x<- edge.connectivity(gr,source=fromnode,target=tonode)
        ec.list[[i]]<-c(fromnode,tonode,x)
        ec.listn[[i]]<-c(V(gr)$name[fromnode],V(gr)$name[tonode],x)
      }
    }
    if(names==F){return(matrix(unlist(ec.list),ncol=3,byrow=T))}
    else{return(matrix(unlist(ec.listn),ncol=3,byrow=T))
    }
}
}
# This function colors all simple paths between the 'fromnode' node to the 
# 'tonode' node using the color specificed by 'clr' stored as the asp_color 
# attribute of the edges
simple.path.coloring<-function(gr,fromnode,tonode,loc=T,clr="purple"){
  paths<-all_simple_paths(gr,fromnode,to=tonode)
  if(length(paths)==0){return("No simple paths between nodes")}
  else{
  if(loc){
  for(i in 1:length(paths)){
    E(gr,path=paths[[i]],dir=T)$clrs<-clr}
return(which(!is.na(E(gr)$clrs)))
  }
else{
  E(gr)$clrs<-"gray"
  for(i in 1:length(paths)){
    E(gr,path=paths[[i]],dir=T)$clrs<-clr}
  return(E(gr)$clrs)
}
  }
}
# SIMPLE PATH (NODE) REDUNDANCY
# Calculates the difference between the longest and shortest simple path
# between nodes.  This inidcates the number of nodes that can be removed
# the set of nodes included in paths between nodes while still maintaining
# a path between nodes.  
simple.path.redundancy<-function(gr,fromnode,tonode){
  paths<-all_simple_paths(gr,fromnode,to=tonode)
  spaths<-shortest_paths(gr,fromnode,to=tonode)
  allpathnodes<-unique(as.vector(unlist(paths)))
  allspathnodes<-as.vector(unlist(spaths))
  return(length(allpathnodes)-length(allspathnodes))
  }
# deg.col.grad() provides an option for coloring nodes based on degree and specified gradient colors
deg.col.grad<-function(gr,dmode="all",gradcols){
  gd<-degree(gr,mode=dmode)
  #heatcol=rev(heat.colors(max(gd)+1))
  #gdcols<-heatcol[gd+1]
  wbpal=colorRampPalette(gradcols)
  gdcols<-wbpal(max(gd)+1)[gd+1]
  
  return(gdcols)
}
# deg.col.heat() provides an option for coloring nodes based on degree and heat.colors()
deg.col.heat<-function(gr,dmode="all"){
  gd<-degree(gr,mode=dmode)
  heatcol=rev(heat.colors(max(gd)+1))
  gdcols<-heatcol[gd+1]
  return(gdcols)
}

# set.bg.black takes an arguement of true or false and will set the plotting background black and adjust colors for barplot visualizations.
set.bg.black<-function(x){
  if(x==T){
    plotlabcol<<-"white"
    totdegcol<<-c("gray30","green")
    indegcol<<-c("gray30","purple")
    outdegcol<<-c("gray30","orange")
    betcol<<-c("gray30","blue")
    clscol<<-c("gray30","magenta")
    ecccol<<-c("gray30","cyan")
    return(par(bg="black"))}
  else{
    plotlabcol<<-"black"
  totdegcol<<-c("black","green")
  indegcol<<-c("black","purple")
  outdegcol<<-c("black","orange")
  betcol<<-c("black","blue")
  clscol<<-c("black","magenta")
  ecccol<<-c("black","cyan")
  return(par(bg="white"))}
}