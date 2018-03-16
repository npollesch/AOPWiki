#####################################################
## FUNCTION: find "origin" and "terminus" vertices
#####################################################

# Adds $KE_PD (KE Path Designator) vertex attribute: 
#    "origin"    if formally described as MIE, OR degree in = 0
#    "terminus"  if formally described as AO, or degree out = 0
#    ""          if not origin or terminus

### INPUT: Must be an igraph object
#   all vertices must have an attributes called "KE_KED" (Key Event Designator), with values of "MIE", "KE", or "AO" 

### OUTPUT: New igraph object that is identical to input object
#   but now how $KE_PD vertex attribute

add_KE_PD<-function(g){
  V(g)$KE_PD<-""
  V(g)$KE_PD[V(g)$KE_KED=="MIE"]<-"origin"
  V(g)$KE_PD[degree(g, mode="in")==0]<-"origin"
  V(g)$KE_PD[V(g)$KE_KED=="AO"]<-"terminus"
  V(g)$KE_PD[degree(g, mode="out")==0]<-"terminus"
  return(g)
}


#####################################################
## FUNCTION: DETERMINE LINEAR AOPS (ie simple paths) between all possible MIE and AO pairs (or "ORIGIN" and "TERMINUS" if desired)
#####################################################

### INPUT must be an igraph object
#   all vertices must have an attributes called "KE_KED" (Key Event Designator), with values of "MIE", "KE", or "AO"
#   all vertices must have a "names" attribute (this is used to identify vertices)
#   optionally, if use_KE_PD=TRUE, then all verticesmust have an attributes called "KE_PD" (Path Designator), with values of "origin", "", or "terminus" 

### OUTPUT is a "List of a List" of all paths (by vertex) between all MIE/AO pairs (or origin/terminus pairs)
#   List objects are named by the MIE/AO pair involved
#   List objects are Lists of all paths (by vertex) between the MIE and AO for which the List is named
#   if "remove.zero.paths=TRUE", all origin/terminus pairs with ZERO linear AOPs are removed from the summary, if FALSE they are left in

# takes about 3 mins to run on my computer for the full wiki

linear.AOPs<- function(g, use_KE_PD=FALSE, remove.zero.paths=TRUE){
  #libraries
  require(igraph)
  
  if(use_KE_PD==TRUE){
    #identify all "origins" and "termini"
    pStart<-V(g)[V(g)$KE_PD=="origin"]
    pEnd<-V(g)[V(g)$KE_PD=="terminus"]
  }else{  
    #identify all "MIEs" and "AOs"
    pStart<-V(g)[V(g)$KE_KED=="MIE"]
    pEnd<-V(g)[V(g)$KE_KED=="AO"]
  }
      
  #list and table to store results
  tempList<-list()
  
  #identify all simple paths between all possible MIE and AO pairs
  for(i in 1:length(pStart)){
    for(j in 1:length(pEnd)){
      tempList[[paste(names(pStart)[i],names(pEnd)[j])]]<-all_simple_paths(g, from=pStart[i], to=pEnd[j], mode="out")
    }
  }
  
  #if remove.zero.paths=TRUE, remove instances where #paths=0 from pathList and pathSummary
  if(remove.zero.paths==TRUE){
    pathList<-list()
    for(i in names(tempList)){
        if(length(tempList[[i]])>0){
          pathList[[i]]<-tempList[[i]]
        }
    }
  }else{
    pathList<-tempList
  }
  
  #result
  return(pathList)
}



####################################
## FUNCTION: Identify paths that have non-adjacent edges (KERs) based on MIE to AO paths
####################################

# algorithm description:
# A simple path "A", contains "non-adjacent" edges, if another longer simple,"B", contains all the same vertices as "A"

# INPUT  must be "List of a List" of Paths (as generated from "linear.AOPs" function)
# OUTPUT a "list of lists" of all simple paths between all "origin" and "terminus" pairs that ONLY CONTAIN adjacent edges

# return.nonAdj=TRUE will return the "list of lists" of simple paths containing NON ADJACENT edges
# return.summary=TRUE will return a data frame summarizing the number adj and non-adj paths between each "origin" to "terminus" pair

remove.nonAdj<-function(pathList, return.nonAdj=FALSE, return.summary=FALSE){
  
  #libraries
  require(igraph)
  
  # divides "pathList" into two lists:
  #   a list containg paths with only "adjacent" edges, and
  #   a list containg paths that have "non-adjacent" edges (though it does not identify WHICH edges are non-adjacent yet)
  
  allAdjList<-list()
  nonAdjList<-list()
  
  for(i in names(pathList)){
    allAdjList[[i]]<-list()
    nonAdjList[[i]]<-list()
    
    # sort all paths for MIE/AO pair i, in order from longest to shortest
    byL<-order(sapply(pathList[[i]], FUN=length), decreasing=TRUE) 
    testSet<-pathList[[i]][byL]   
    
    # moves through the "testSet" paths until they are all sorted into the "adjList" or "nonAdjlist"
    while(length(testSet)>0){
      
      # if only one path left, it is moved to adjList
      if(length(testSet)==1){
        allAdjList[[i]][[length(allAdjList[[i]])+1]]<-testSet[[1]]
        testSet<-testSet[-1]
      
      }else{
        hasNon<-vector()
        for(j in 2:length(testSet)){
          
          # if path "j" is shorter AND all if its vertices are contained within path "1", it contains non-adjacnet edges (though which edges specifically are not identified)
          if(length(testSet[[j]]) < length(testSet[[1]]) & all(testSet[[j]]%in%testSet[[1]])){  
            nonAdjList[[i]][[length(nonAdjList[[i]])+1]]<-testSet[[j]]
            hasNon<-c(hasNon,j)
          }
        }
        
        # move top path in testList to adjList and remove all nonAdj paths identifed from testList. Repeat until testList is empty
        allAdjList[[i]][[length(allAdjList[[i]])+1]]<-testSet[[1]]
        testSet<-testSet[-c(1,hasNon)]
      }
    }
  }
  
  # Create sumamry table
  adjSummary<-data.frame(
    Start=sub(" .*", "", names(allAdjList)),
    End=sub(".* ", "", names(allAdjList)),
    adjPaths=sapply(allAdjList,length),
    nonAdjPaths=sapply(nonAdjList,length),
    row.names=NULL)
  
  # result
  if(return.summary==TRUE){
    return(adjSummary)
  }else{
    if(return.nonAdj==TRUE){
      return(nonAdjList)      
    }else{
      return(allAdjList)
    }
  }
}
 


##############################################
## FUNCTION: CREATE EDGE LIST FROM PATHS
##############################################

## INPUT  is a vector that lists the sequence of vertices in a path (from start to end if a directed path)
## OUTPUT is a 2-column matrix, each row representing the start and end vetices of an edge

edge_from_path<-function(x, by.vertex.name=TRUE){
  require(igraph)
  
  edgeList<-vector()
  for(i in 1:(length(x)-1)){
    edgeList<-c(edgeList, x[i], x[i+1])
  }
  if(by.vertex.name==TRUE){
    edgeList<-matrix(names(edgeList), ncol=2, byrow=TRUE)
  }else{
    edgeList<-matrix(edgeList, ncol=2, byrow=TRUE)
  }
  return(edgeList)
}



##############################################################################
###  FUNCTION: Create igraph object from a "list of lists" of MIE to AO paths
###  (as generated from linear.AOPs() or remove.nonAdj() )
##############################################################################

graph_from_pathList<-function(pathList, by.vertex.name=TRUE){
  require(igraph)
  
  #use edge_from_path function to create edgelists
  pathEdges<-list()
  for(i in 1:length(pathList)){
    pathEdges[[i]]<-lapply(pathList[[i]], edge_from_path, by.vertex.name=by.vertex.name)  
  }
  names(pathEdges)<-names(pathList)
  
  #combine all edges into one edgelist, and remove redundant edges
  edgeList<-lapply(pathEdges, function(x) do.call(rbind,x))
  edgeList<-do.call(rbind,edgeList)
  edgeList<-unique(edgeList)
  
  #create graph from edgeList
  g<-graph_from_edgelist(edgeList)
  return(g)
}



###################################################################
## FUNCTION: if g2 is a subgraph of g1, this function identifies 
##           which edges have been removed from  g2 compared to g1
###################################################################

#INPUT: g2 must be a subgraph of g1
#OUTPUT: returns an edgelist of the edges in g1 that are not in g2

edge_difference<-function(g2,g1){
  require(igraph)
  require(prodlim)
  
  g2E<-as.data.frame(as_edgelist(g2))
  g1E<-as.data.frame(as_edgelist(g1))
  edgeDiff<-g1E[-row.match(g2E,g1E),]
  return (edgeDiff)
}



###################################################################
## FUNCTION: Condense all strongly connected components in a graph
###################################################################

contract.scc<-function(g, V.resize=TRUE){
  sccMap<-components(g, mode="strong")$membership
  ntcomps<-which(components(g, mode="strong")$csize>1)
  
  # identify if the strong components have MIE/AO or origin/terminus attributes
  has_M<-vector()
  has_A<-vector()
  no_MA<-vector()
  has_O<-vector()
  has_T<-vector()
  no_OT<-vector()
  for(i in ntcomps){
    if(any(V(g)$KE_KED[sccMap==i]=="MIE")){
      has_M<-c(has_M,i)
    }else{
      if(any(V(g)$KE_KED[sccMap==i]=="AO")){
        has_A<-c(has_A,i)
      }else{
        no_MA<-c(no_MA,i)
      }
    }
    if(any(V(g)$KE_PD[sccMap==i]=="origin")){
      has_O<-c(has_O,i)
    }else{
      if(any(V(g)$KE_PD[sccMap==i]=="terminus")){
        has_T<-c(has_T,i)
      }else{
        no_OT<-c(no_OT,i)
      }
    }
  }
  
  # tells contract() function that if merged vertices have plotting coordinates, take the mean of them
  vAttComb<-list(
    plotX = "mean",
    plotY = "mean",
    toString
  )
  
  #contract the vertices in strongly connected components
  g2<-contract(g,mapping=sccMap, vertex.attr.comb = vAttComb)
  # remove redundant edges and resulting loops
  newG<-simplify(g2, remove.multiple=TRUE, remove.loops=TRUE)
  
  # assign attributes to contracted vertices
  V(newG)$col[ntcomps]<-"purple"    
  V(newG)$size<-2.5
  V(newG)$size[ntcomps]<-3.5
  V(newG)$KE_KED[has_M]<-"MIE"
  V(newG)$KE_KED[has_A]<-"AO"
  V(newG)$KE_KED[no_MA]<-"KE"
  V(newG)$KE_PD[has_O]<-"origin"
  V(newG)$KE_PD[has_T]<-"terminus"
  V(newG)$KE_PD[no_OT]<-""
  return(newG)
}


