
#' FixedPop
#'
#' Designates a FixedPop, a structure whose size remains constant throughout the course of a simulation
#'
#' @param tree DiffTree object to add FixedPop to
#' @param name population name
#' @param size initial population size
#' @param label initial probabilty that a cell receives a unique barcode
#'
#' @export
#' @examples
#' \dontrun{
#' FixedPop(myTree, "Population A", 1000, 0.5)
#' }
FixedPop = function(tree, name, size, label = 0){
  # # If a tree doesn't exist, create one
  # if(!exists("myTree")){
  #   myTree <<- graph.empty(n = 0, directed = T)
  # }

  # If the tree already contains the population, alert the user, else: add it
  if(name %in% V(tree)$popName){
    stop(paste("Population", name, "already exists.  Please use unique population names."))
    #print(paste("Population", name, "already exists.  Please use unique population names."))
  } else {
    temp_tree <- tree +  vertex(name, type = "FixedPopCell", popName = name, root = FALSE, popsize = size, lbl = label)
  }
  eval.parent(substitute(tree<-temp_tree))
}

#' GrowingPop
#'
#' Designates a GrowingPop, a structure used to simulate differentiation via a branching process
#'
#' @param tree DiffTree object to add GrowingPop to
#' @param name population name
#' @param size initial population size
#' @param label initial probabilty that a cell receives a unique barcode
#'
#' @export
#' @examples
#' \dontrun{
#' GrowingPop(myTree, "Population A", 1000, 0.5)
#' }
GrowingPop = function(tree, name, size, label = 0.0){
  # # If a tree doesn't exist, create one
  # if(!exists("myTree")){
  #   myTree <<- graph.empty(n = 0, directed = T)
  # }

  # If the tree already contains the population, alert the user, else: add it
  if(name %in% V(tree)$popName){
    stop(paste("Population", name, "already exists.  Please use unique population names."))
  } else {
    temp_tree <- tree +  vertex(name, type = "CellPopulation", popName = name, root = FALSE, popsize = size, lbl = label)
  }
  eval.parent(substitute(tree<-temp_tree))
}

#' DiffTriangle
#'
#' Designates a DiffTriangle, a fixed-size structure used to model a population undergoing different levels of maturation
#'
#' @param tree DiffTree object to add DiffTriangle to
#' @param name population name
#' @param height number of levels in structure, corresponding to number of levels until full maturation
#' @param first_level number of triangles to be stacked side-by-side or the number of cells at the first level of maturation
#' @param label initial probabilty that a cell receives a unique barcode
#'
#' @export
#' @examples
#' \dontrun{
#' DiffTriangle(myTree, "Differentiated Population A", 12, 1000, 0.5)
#' }
DiffTriangle = function(tree, name, height, first_level = 1, label = 0.0){
  # # If a tree doesn't exist, create one
  # if(!exists("myTree")){
  #   myTree <<- graph.empty(n = 0, directed = T)
  # }

  # If the tree already contains the population, alert the user, else: add it
  if(name %in% V(tree)$popName){
    stop(paste("Population", name, "already exists.  Please use unique population names."))
  } else {
    temp_tree <- tree +  vertex(name, type = "DiffTriangle", height = height, firstlevel = first_level, popName = name, root = FALSE, lbl = label)
  }
  eval.parent(substitute(tree<-temp_tree))
}


#' DiffTree
#'
#' Creates an empty DiffTree object from which to build a hierarchy
#'
#' @export
#' @examples
#' \dontrun{
#' myTree = DiffTree()
#' }
DiffTree = function(){
  # # If a tree doesn't exist, create one
  # if(!exists("myTree")){
  #   myTree <<- graph.empty(n = 0, directed = T)
  # }

  # # If the tree already contains the population, alert the user, else: add it
  # if(name %in% V(myTree)$popName){
  #   stop(paste("Population", name, "already exists.  Please use unique population names."))
  # } else {
  #   myTree <<- myTree +  vertex(name, type = "DiffTree", popName = name, root = FALSE)
  # }
  return(graph.empty(n = 0, directed = T))
}

#' setRoot
#'
#' Used to designate furthest upstream population in the hierarchy
#'
#' @param tree DiffTree object to modify
#' @param popName name of furthest upstream population
#'
#' @export
#' @examples
#' \dontrun{
#' setRoot(myTree, "Population A")
#' }
setRoot = function(tree, popName){
  myTree = tree
  # If the population exists in the tree, set its Root attribute to True, else: alert user
  if(popName %in% V(myTree)$popName){
    for(v in V(myTree)){
      if (V(myTree)[v]$popName == popName){
        myTree <- set.vertex.attribute(graph = myTree, name = "root", index = v,value = TRUE)
        myTree <- set.graph.attribute(graph = myTree, name = "rootInd", value = v)
      } else {
        myTree <- set.vertex.attribute(graph = myTree, name = "root", index = v,value = FALSE)
      }
    }
  } else {
    stop(paste("Population", popName, "does not exist.  Please create it before setting root."))
  }
  eval.parent(substitute(tree<-myTree))
}

#' setFitnessDistribution
#'
#' Used to designate distribution from which changes in fitness for a new clone are drawn
#'
#' @param tree DiffTree object for which to set fitness distribution
#' @param distribution random distribution from which to draw ("normal", "doubleexp", "uniform")
#' @param is_random boolean variable set to TRUE if random fitness distribution is specified
#' @param alpha_fitness alpha parameter for fitness distribution
#' @param beta_fitness beta parameter for fitness distribution
#' @param pass_prob probability that new clone developed a passenger mutation that does not affect its fitness
#' @param upper_fitness upper limit on clone fitness
#' @param lower_fitness lower limit on clone fitness
#'
#' @export
#' @examples
#' \dontrun{
#' setFitnessDistribution(tree = myTree, distribution = "uniform", is_random = T, alpha_fitness = 0, beta_fitness = 1, pass_prob = 0, upper_fitness = 5, lower_fitness = 0)
#' }
setFitnessDistribution = function(tree, distribution = "normal", is_random = T, alpha_fitness = 0, beta_fitness = 1, pass_prob = 1, upper_fitness = NA, lower_fitness = 0){
  # If a tree doesn't exist, create one
  myTree = tree
  myTree <- set.graph.attribute(graph = myTree, name = "distribution", value = distribution)
  myTree <- set.graph.attribute(graph = myTree, name = "is_random", value = is_random)
  myTree <- set.graph.attribute(graph = myTree, name = "alpha_fitness", value = alpha_fitness)
  myTree <- set.graph.attribute(graph = myTree, name = "beta_fitness", value = beta_fitness)
  myTree <- set.graph.attribute(graph = myTree, name = "pass_prob", value = pass_prob)
  myTree <- set.graph.attribute(graph = myTree, name = "upper_fitness", value = upper_fitness)
  myTree <- set.graph.attribute(graph = myTree, name = "lower_fitness", value = lower_fitness)
  eval.parent(substitute(tree<-myTree))
}

#' addEdge
#'
#' Used to designate transition between populations
#'
#' @param tree DiffTree to add transition to
#' @param parent parent population name ("from" population)
#' @param child child population name ("to" population)
#' @param type transition or event type ("alpha", "beta", "gamma1", "gamma2", "delta", "zeta", "mutation")
#' @param rate expected number of events per cell per time unit
#'
#' @export
#' @examples
#' \dontrun{
#' addEdge(myTree, "Population A", "Population A", "alpha", 0.01)
#' addEdge(myTree, "Population A", "Population B", "gamma1", 0.05)
#' addEdge(myTree, "Population A", "Population A", "mutation", 1e-7)
#' }
addEdge = function(tree, parent, child, type, rate){
  myTree = tree
  c = checkEdge(myTree, parent, child, type) > 0
  if(rate < 0)
    stop("Rates must be non-negative")
  if(c == 0){
    myTree <- myTree + edge(parent, child, weight = rate, type = type)
  } else if (c == 1){
    stop("Edge already exists. Please use 'editEdit' to change its rate")
  }
  eval.parent(substitute(tree<-myTree))
}

# Returns -1 if populations don't exist
# returns 0 if edge does not exist
# returns 1 if edge exists
checkEdge = function(tree, parent, child, type){
  if(checkPopulation(tree, parent) & checkPopulation(tree, child)){
    el = getEdges(tree)
    edges = el[el$parent == parent & el$child == child & el$type == type,]
    #print(edges)
    return(nrow(edges) == 1)
  }
  return(-1)
}

#' @export
editEdge = function(tree, parent, child, type, rate){
  myTree = tree
  c = checkEdge(myTree, parent, child, type)
  if (c == 0){
    stop("Edge does not already exist.  Add it using 'addEdge'.")
  } else if(c == 1){
    for(e in E(myTree)){
      if(ends(myTree, e)[1] == parent & ends(myTree, e)[2] == child & E(myTree)[e]$type == type){
        myTree <- set_edge_attr(myTree, "weight", e, rate)
      }
    }
  }
  eval.parent(substitute(tree<-myTree))
}

getEdges = function(tree){
  myTree = tree
  g = myTree
  edge_list = data.frame()
  for (e in E(myTree)){
    edge_list = rbind(edge_list, cbind(ends(g,e), E(g)[e]$type, E(g)[e]$weight))
  }
  if(nrow(edge_list) > 0)
    colnames(edge_list) = c("parent", "child", "type", "rate")
  return(edge_list)

}

checkTree = function(){
  if(!exists("myTree")){
    stop("Tree does not exist.  Please create at least one cell population.")
    return(FALSE)
  }
  return(TRUE)
}

checkPopulation = function(tree, popName){
  myTree = tree
  if(popName %in% V(myTree)$popName){
    return(TRUE)
  }

  stop(paste("The population", popName, "does not exist.  Please create it first."))
  return(FALSE)
}

checkRoot = function(tree){
  myTree = tree
  if(!is.null(get.graph.attribute(graph = myTree, "rootInd"))){
    return(TRUE)
  }

  return(FALSE)
}

#' plotTree
#'
#' Used to plot hierarchy structure using igraph
#'
#' @param tree DiffTree object to plot
#' @param allEdge boolean variable set to TRUE to plot all transitions, set to FALSE to only show gamma1, gamma2, and beta edges
#'
#' @export
#' @examples
#' \dontrun{
#' plotTree(myTree, allEdge = TRUE)
#' }
plotTree = function(tree, allEdge = F){
  myTree = tree
  if(!allEdge){
    g.plot = delete.edges(myTree, subset(E(myTree), !E(myTree)$type %in% c("gamma1", "gamma2", "beta")))
    #g.plot = delete.edges(myTree, !which(E(myTree)$type %in% c("gamma1", "gamma2", "beta")))
  } else {
    g.plot = myTree
  }
  if(checkRoot(myTree)){
    rootInd = get.graph.attribute(graph = myTree, "rootInd")
    plot(g.plot, layout=layout_as_tree(g.plot, root = c(rootInd)))
  } else {
    plot(g.plot)
  }
}

#' writeTree
#'
#' Used to create tree files needed by C++ simulation
#'
#' @param tree DiffTree object to write
#' @param outdir output directory for files
#'
#' @export
#' @examples
#' \dontrun{
#' writeTree(myTree, outdir = "C:/diffpop/example1/")
#' }
writeTree = function(tree, outdir = "tree_files/"){
  options(scipen = 999)
  myTree = tree
  g = myTree
  cdir = getwd()

  directory = R.utils::getAbsolutePath(outdir)

  #outdir = paste(getwd(), outdir, sep = "")
  dir.create(directory, showWarnings = F, recursive = T)
  setwd(directory)
  #write.graph(g, "g.txt", format = "graphml")
  #write.table(cbind(V(g)$name, V(g)$type, V(g)$root, V(g)$popsize, V(g)$lbl, V(g)$height, V(g)$firstlevel), "nodes.txt", quote = F, row.names = F, col.names = F)
  write.table(cbind(V(g)$name, V(g)$type, V(g)$root, V(g)$popsize, V(g)$lbl, V(g)$height, V(g)$firstlevel), "nodes.csv", quote = F, row.names = F, col.names = F, sep = ",")
  # write.table(cbind(V(g)$name), "names.txt", quote = F, row.names = F, col.names = F)
  # write.table(cbind(V(g)$type), "types.txt", quote = F, row.names = F, col.names = F)
  # write.table(cbind(V(g)$root), "roots.txt", quote = F, row.names = F, col.names = F)
  # write.table(cbind(V(g)$popsize), "sizes.txt", quote = F, row.names = F, col.names = F)
  # write.table(cbind(V(g)$lbl), "labels.txt", quote = F, row.names = F, col.names = F)
  # write.table(cbind(V(g)$height), "heights.txt", quote = F, row.names = F, col.names = F)
  # write.table(cbind(V(g)$firstlevel), "firstlevels.txt", quote = F, row.names = F, col.names = F)
  write.table(row.names(cbind(bfs(myTree, which(V(myTree)$root == T))$order)), file = "bfs.txt", row.names = F, col.names = F, quote = F)

  edge_list = data.frame()
  for (e in E(g)){
    #print(typeof(e))
    #print(ends(g, e))
    #print(E(g)[e]$type)
    edge_list = rbind(edge_list, cbind(ends(g,e), E(g)[e]$type, E(g)[e]$weight))
  }

  #write.table(edge_list, "edges.txt", row.names = F, col.names = F, quote = F)
  write.table(edge_list, "edges.csv", row.names = F, col.names = F, quote = F, sep  = ",")
#
#   write.table(edge_list[,1], "parent.txt", row.names = F, col.names = F, quote = F)
#   write.table(edge_list[,2], "child.txt", row.names = F, col.names = F, quote = F)
#   write.table(edge_list[,3], "diff.txt", row.names = F, col.names = F, quote = F)
#   write.table(edge_list[,4], "rates.txt", row.names = F, col.names = F, quote = F)

  write.table(rbind(g$distribution, g$is_random, g$alpha_fitness, g$beta_fitness, g$pass_prob, g$upper_fitness, g$lower_fitness),
              "fitnessdist.txt", row.names = F, col.names = F, quote = F)
  setwd(cdir)
}

#' @export
generateExample = function(n_LT = 100){
  temp = DiffTree()
  FixedPop(temp, "LT", as.integer(n_LT), .11)
  FixedPop(temp, "ST", as.integer(2.9*n_LT), .65)
  FixedPop(temp, "MPP", as.integer(9*n_LT), .88)
  FixedPop(temp, "CLP", as.integer(13*n_LT), .88)
  FixedPop(temp, "CMP", as.integer(39*n_LT), 0.88)

  addEdge(temp, "LT", "ST", "gamma1", 0.009)
  addEdge(temp, "ST", "MPP", "gamma1", 0.045)
  addEdge(temp, "MPP", "CLP", "gamma1", 0.022)
  addEdge(temp, "MPP", "CMP", "gamma1", 3.992)
  addEdge(temp, "LT", "LT", "alpha", 0.009)
  addEdge(temp, "ST", "ST", "alpha", 0.042)
  addEdge(temp, "MPP", "MPP", "alpha", 4.000)
  addEdge(temp, "CLP", "CLP", "alpha", 3.980)

  setRoot(temp, "LT")
  return(temp)
  #plotTree()
  #writeTree()
}

removeTree = function(){
  #print(exists("myTree"))
  if(exists("myTree"))
  {
    objs <- ls(pos = ".GlobalEnv")
    rm(list = objs[grep("myTree", objs)], pos = ".GlobalEnv")
  }
}

#' simulateTree
#'
#' Used to plot hierarchy structure using igraph
#'
#' @param tree DiffTree object to simulate
#' @param fixed boolean.  TRUE if simulating a fixed population hierarchy using FixedPops and DiffTriangles, FALSE if simulating a dynamic-sized population hierarchy using GrowingPops
#' @param time number of time units to simulate
#' @param census how often to output entire system population to filesystem
#' @param indir input directory for tree files
#' @param outdir output directory
#' @param seed numeric seed for random number generator
#'
#' @export
#' @examples
#' \dontrun{
#' simulateTree(myTree, fixed = TRUE, time = 500, census = 1, indir = "input_dir", outdir = "output_dir")
#' }
simulateTree = function(tree, fixed = FALSE, time = 100, census = -1, indir = ".", outdir = ".", seed = NULL){
  # Write tree files
  writeTree(tree, indir)

  # Create directory for output files
  directory = R.utils::getAbsolutePath(outdir)
  dir.create(directory, showWarnings = F, recursive = T)

  # Call appropriate
  if(fixed){
    diffpop:::simulateFixedTreeCodeNew(nObs = time,
                                    traverseFrequency = census,
                                    indir = paste(R.utils::getAbsolutePath(indir), "/", sep = ""),
                                    outdir = paste(R.utils::getAbsolutePath(outdir), "/", sep = ""),
                                    seed = seed)
  } else {
    diffpop:::simulateTreeCodeNew(nObs = time,
                               traverseFrequency = census,
                               indir = paste(R.utils::getAbsolutePath(indir), "/", sep = ""),
                               outdir = paste(R.utils::getAbsolutePath(outdir), "/", sep = ""),
                               seed = seed)
  }
}



