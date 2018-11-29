# load libraries for parallelization
library(foreach)
library(doParallel)

cores = detectCores()
cl = makeCluster(cores[1]-1) 
registerDoParallel(cl)

ntrials = 1000

# Iterate over number of trials
foreach(i_=1:ntrials) %dopar%{
  print(i_)
  library(diffpop)
  
  # Simulation size and label parameter
  nLT = 5000
  LT_lbl = 0.01
  
  # Blank DiffTree object
  tree = DiffTree()
  
  # Add all populations to tree using population sizes estimated from Busch et al. 
  FixedPop(tree, "LT", nLT, LT_lbl)
  FixedPop(tree, "ST", as.integer(2.9*nLT), 0.0)
  FixedPop(tree, "MPP", as.integer(9*nLT), 0.0)
  FixedPop(tree, "CMP", as.integer(39*nLT), 0.0)
  FixedPop(tree, "CLP", as.integer(13*nLT), 0.0)
  FixedPop(tree, "GMP", as.integer(0.24*39*nLT), 0.0)
  FixedPop(tree, "MEP", as.integer(0.39*39*nLT), 0.0)
  FixedPop(tree, "proB", as.integer(108*13*nLT), 0.0)
  
  # Add self-renewal/mitosis events
  addEdge(tree, "LT", "LT", "alpha", 0.009)
  addEdge(tree, "ST", "ST", "alpha", 0.042)
  addEdge(tree, "MPP", "MPP", "alpha", 4)
  addEdge(tree, "CLP", "CLP", "alpha", 3.00)
  addEdge(tree, "CMP", "CMP", "alpha", 4)
  
  # Add differentiation events
  # Note: Busch et al. assume mitosis-independent differentiation
  addEdge(tree, "LT", "ST", "gamma1", 0.009)
  addEdge(tree, "ST", "MPP", "gamma1", 0.045)
  addEdge(tree, "MPP", "CLP", "gamma1", 0.022)
  addEdge(tree, "MPP", "CMP", "gamma1", 3.992)
  addEdge(tree, "CLP", "proB", "gamma1", 2.000)
  addEdge(tree, "CMP", "GMP", "gamma1", 2)
  addEdge(tree, "CMP", "MEP", "gamma1", 3)
  
  # Set LT population as root of tree
  setRoot(tree, "LT")
  
  # Simulate tree for 800 time units (days)
  # Note: we use FixedPops here because parameters were 
  #       estimated for steady state hematopoiesis
  simulateTree(tree = tree,
               fixed = TRUE,
               time = 800,
               indir = paste("app_input/", i_, sep = ""),
               outdir = "app_output/",
               census = -1)
  
}

stopCluster(cl)
