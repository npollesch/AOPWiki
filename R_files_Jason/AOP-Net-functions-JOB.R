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
#   but now with $KE_PD vertex attribute

add_KE_PD<-function(g){
  V(g)$KE_PD<-""
  V(g)$KE_PD[V(g)$KE_KED=="MIE"]<-"origin"
  V(g)$KE_PD[degree(g, mode="in")==0]<-"origin"
  V(g)$KE_PD[V(g)$KE_KED=="AO"]<-"terminus"
  V(g)$KE_PD[degree(g, mode="out")==0]<-"terminus"
  return(g)
}


####################################################################
## FUNCTION: DETERMINE LINEAR AOPS (ie simple paths) between all
## possible MIE and AO pairs (or "ORIGIN" and "TERMINUS" if desired)
####################################################################

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



##########################################################
### FUNCTION: Identify non-adjacent KERs in an AOP network
##########################################################

### Algorithm Description
#     STEP 1: Potenial non-adjacent KERs identified when distance between KEup and KEdown can be greater than two KEs
#     STEP 2: Longest unique paths between "origin" to "terminus" pairs identified
#     STEP 3: A KER is non-adjacent if its KEs have a possible distance >2 (Step 1) AND DOES NOT occur in any longest unique path (Step 2)
#
### Input: an igraph object with:
#     Vertex attribute called KE_KED with values "MIE", "KE", or "AO"
#
### Output: identical igraph object as input, with new edge attribute called "adjacency"


add_KER_adjacency<-function(g){
  
  # libraries
  require(igraph)

  
  ### STEP 1: check all edges for when possible simple path length >2
  
  allE<-as_edgelist(g, names=FALSE)
  maybeNon<-vector()
  for(i in 1:nrow(allE)){
    sPaths<-all_simple_paths(g, from=allE[i,1], to=allE[i,2])
    if(all(sapply(sPaths,length)<3)){
      maybeNon<-c(maybeNon,0)
    }else{
      maybeNon<-c(maybeNon,1)
    }
  }
  maybeNon_E<-allE[maybeNon==1,] 
  maybeNon_E<- matrix(V(g)[maybeNon_E]$name, ncol=2) # Edge list of "potentially" non-adjacent KERs


  ### Step 2: identify longest unique paths between "origin" to "terminus" pairs
  #     Step A: identify all simple paths between all origin and terminus pairs using linear.AOPs() function
  #     Step B: A path is a "longest unique path" if there is NO OTHER LONGER path that contains all the same KEs

  
  #~~ Step 2a ~~~~
  
  # add KE_PD vertex attributes (identifies "origin" and "terminus" KEs)
  gOT<-add_KE_PD(g)
  
  # identify all simple paths between Origin/Terminus pairs
  p<-linear.AOPs(gOT, use_KE_PD=TRUE)
  
  
  #~~ Step 2b ~~~
  
  # identif all longest unique paths:
  longPaths<-list()
  for(i in names(p)){
    longPaths[[i]]<-list()
    
    # sort all paths for O/T pair i, in order from longest to shortest
    byL<-order(sapply(p[[i]], FUN=length), decreasing=TRUE) 
    testSet<-p[[i]][byL]
    
    # moves through the "testSet" paths until they are all checked against shorter paths
    while(length(testSet)>0){
      
      # if only one path left, it is moved to longPaths list
      if(length(testSet)==1){
        longPaths[[i]][[length(longPaths[[i]])+1]]<-testSet[[1]]
        testSet<-testSet[-1]
      
      # else check if all KEs on shorter paths are contained within path "1"
      # if so, it is NOT a "longest unique path"
      }else{
        hasShort<-vector()
        for(j in 2:length(testSet)){
          if(length(testSet[[j]]) < length(testSet[[1]]) & all(testSet[[j]]%in%testSet[[1]])){  
            hasShort<-c(hasShort,j)
          }
        }
        
        # move top path in testList to longPAths and remove all short paths identifed from testList. Repeat until testList is empty
        longPaths[[i]][[length(longPaths[[i]])+1]]<-testSet[[1]]
        testSet<-testSet[-c(1,hasShort)]
      }
    }
  }
  
  #Convert longPaths into an edge list using edge_from_path() function
  lp_E_temp<-list()
  for(i in 1:length(longPaths)){
    lp_E_temp[[i]]<-lapply(longPaths[[i]], edge_from_path, by.vertex.name=TRUE)  
  }
  lp_E<-lapply(lp_E_temp, function(x) do.call(rbind,x))
  lp_E<-do.call(rbind,lp_E)
  lp_E<-unique(lp_E)  # Edge list of KERs occur longest unique paths
  
  
  ### Step 3
  
  # compare "maybeNon_E" to "lp_E" edgeList
  #   if edge from "maybeNon_E"  DOES NOT also occur in lp_E, then it is NON-ADJACENT
  lp_E<-as.data.frame(lp_E, stringsAsFactors=FALSE)
  maybeNon_E<-as.data.frame(maybeNon_E, stringsAsFactors=FALSE)
  nonAdjE<-maybeNon_E[is.na(row.match(maybeNon_E, lp_E)),]
  nonAdjE<-as.vector(as.character(t(nonAdjE))) #format so that list can be used to call edges using E(g)
  
  # FINAL: add edge attribute $adjacency, with value "adjacent" or "non-adjacent"
  newG<-g
  E(newG)$adjacency<-"adjacent"
  E(newG, P=nonAdjE)$adjacency<-"non-adjacent"
  return(newG)

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


