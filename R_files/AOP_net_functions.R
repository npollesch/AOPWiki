
# color.comps is a function to color all the strongly connected components in a graph.  In theory it should do the same thing as cycle.color below
# however, I could not get cycle.color to work properly on the AOPwiki graph, so I made this one.  -Nate 5/11/17
color.comps<-function(gr, ccmode="Strong"){ #function to color all non-trivial strongly connected components in a graph
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
  return(list(vcol=V(gr)$color,ecol=E(gr)$color))
}

# cycle.color function can be used to color SCC edges and vertices for visualization
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


#condense.map function determines condensation mapping needed to contract vertices in igraph
#output is a vector map to be used in 'condense.graph' function
condense.map = function(gr,suppress=FALSE){
  if(is.dag(gr)){print("Graph is acyclic")}
  else{gsccs<-components(gr, mode="strong")
  return(gsccs$membership)
  }
}

#condense.graph function takes the condensation vector map and original graph and creates the condensed graph
#output is a new igraph object that has strongly connected components groups into supernodes 
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

# this function computes a topological ordering layout for visualizaiton of an igraph object
topo.layout = function(gr){
  tsort<-topo_sort(gr, mode = c("out"))
  tpos<-match(V(gr)$name,tsort$name) # finds position of node names in topoloigcal ordering
  tlay<-cbind(tpos,tpos) #define a plot layout that orders accoring to topological order
  a<-paste(deparse(substitute(gr)),"topo","lay", sep = ".") # creates name for output from function input
  x<-assign(a, tlay, envir=globalenv()) # creates layout global variable
  return(x) #makes the layout object the output
}

# function computes path counts for nodes involved in simple paths between MIEs and AOs.  
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
        x<- all_simple_paths(gr,fromnode, to = tonode)
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

##NOTE: SET VERTEX ATTR NOT WORKING IN FUNCTION YET - PERHAPS JUST CREATE 2D OUTPUT LIKE COLORING FUNCTION ####

#~TASK: NEED TO FIX THIS FUNCTION - PERHAPS SOME PATHS DON"T EXIST BETWEEN MIE AND ALL AOS? ####
