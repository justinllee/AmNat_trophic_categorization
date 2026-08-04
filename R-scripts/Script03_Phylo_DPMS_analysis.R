### Script objectives:
### Estimate dietary niche states of Elapids using phylogenetic implementation
### of Dirichlet Process Multinomial State (DPMS) model

## Prepare and install packages
library(macroevolution)
library(ape)
library(RColorBrewer)
library(phytools)

## Prepare directories
path <- "/PATH/"
path <- setwd(path)

## Import phylogeny
phy <- macroevolution::read.newick("/PATH/Datasets/All_Elapids_tree_v2.tre")
tree = read.tree("/PATH/Datasets/All_Elapids_tree_v2.tre")

## Prepare category schemes and associated datasets
schemes <- c("prey_traditional", "prey_phylo1", "prey_phylo2", 
             "prey_functional1", "prey_functional2", "prey_type")
datafiles <- sprintf("bhc_dataset_%s.csv", schemes)
datafiles <- datafiles[1:6]
setwd(path)

##-----Part 1: Begin model ------
for(CA in 1:length(datafiles)){
  setwd(path) #resets working directory after each categorization scheme is finished
  snakediet <- read.csv(datafiles[ CA ])
  ncat <- length(unique(snakediet$food))
  x = xtabs(snakediet$count ~ snakediet$species + snakediet$food)
  
  folder <- sprintf("BHC_%s", schemes[ CA ])
  dir.create(folder) #creates a new folder for each categorization scheme
  setwd(folder)
  
  #Perform method -------------------------------------------------------
  NSTATE = 1000
  BETA=1
  niter=10000 
  thin=10
  
  newdir <- sprintf("results_%s", schemes[CA]) #output file for BHC results
  dir.create(newdir)
  
  mcmc = macroevolution::make.rcm.dmm(phy, snakediet, NSTATE, beta.init=BETA)
  output.file = sprintf("results_%s/mcmc_%s_%s.out", schemes[CA], NSTATE, BETA)
  mcmc(
    niter=niter,
    thin=thin,
    output.file=output.file) 
  
  Exfile = sprintf("results_%s/Enodestate_%s_%s.rds", schemes[CA], NSTATE, BETA)
  resultfile = output.file
  
  out = macroevolution::read.rcm.dmm(resultfile)
  map = which.max(rowSums(out$pars[, 1:2]))
  obj = make.rcm.dmm.from.sample(map, resultfile)

  diet_states <- obj$tip.state
  diet_props <- obj$dens.map

  obj = macroevolution::make.rcm.dmm.from.sample(500, resultfile)
  
  # Extract multinomial state parameters -------------------------------------------------------
  ntip = macroevolution::Ntip(phy)
  nnode = macroevolution::Nnode(phy)
  
  X  = matrix(0, nnode, ncol(obj$dens.mean), dimnames=list(
    NULL, colnames(obj$dens.mean)))
  
  begin = nrow(out$pars)-2500
  end = nrow(out$pars)
  
  for (i in begin:end) {
    obj = macroevolution::make.rcm.dmm.from.sample(i, resultfile) 
    X[1:ntip, ] = X[1:ntip, ] + obj$dens.map[obj$tip.state, ] / end
    X[macroevolution::root(phy):nnode, ] =  
      X[macroevolution::root(phy):nnode, ] +
      t(apply(obj$asr, 2, function(mp) {
        colSums(sweep(obj$dens.map, 1, mp, "*"))})) / 2501
  }
  
  saveRDS(X, file=Exfile)  
  
  z = readRDS(sprintf("results_%s/Enodestate_%s_%s.rds", schemes[CA], NSTATE, BETA))
  nrow(z);
  ncol(z)  
  ncat <- ncol(z) #in case the categories change during the analysis
  
  # Get diet states for each tree tip -------------------------------------------------------
  out = read.rcm.dmm(resultfile)  
  phy = out$phy
  
  map = which.max(rowSums(out$pars[, 1:2]))
  
  obj = make.rcm.dmm.from.sample(map, resultfile)

  node.state = c(obj$tip.state, apply(obj$asr, 2, which.max))
  
  node.state = obj$smap(1)[,1]
  y = obj$dens.map[node.state, ]
  
  df <- cbind(macroevolution::tiplabels(phy), node.state[1:490])
  write.csv(df, sprintf("%s_spp_dietstates_490.csv", schemes[CA]), row.names = F)
  
  df1 <- cbind(macroevolution::tiplabels(phy), node.state[1:245])
  write.csv(df1, sprintf("%s_spp_dietstates_245_max_prob.csv", schemes[CA]), row.names = F)
  
  df2 <- cbind(macroevolution::tiplabels(phy), node.state[245:490])
  write.csv(df2, sprintf("%s_spp_dietstates_245_asr_prob.csv", schemes[CA]), row.names = F)
  
  # Gets diet niche states at tree tips
  states <- read.csv(sprintf("%s_spp_dietstates_245_max_prob.csv", schemes[CA]))
  
  # Gets ancestral diet niche states
  states2 <- read.csv(sprintf("%s_spp_dietstates_245_asr_prob.csv", schemes[CA]))
  
  # Gets proportions of each prey item for each diet niche state
  state_props <- as.data.frame(diet_props)
  write.csv(state_props, sprintf("%s_inferred_states.csv", schemes[CA]), row.names = F)
  
  # Get diet proportions for each species -------------------------------------------------------
  for (NSTATE in c(20L, 50L, 100L, 1000L)){
    mcmc = macroevolution::make.rcm.dmm(phy, snakediet, NSTATE)
    mcmc(niter=2^15, thin=2^4, output.file=sprintf("%s_mcmc-%s.out", schemes[CA], NSTATE))
  }
  
  X = vector("list", 4)
  w = numeric(4)
  names(w) = names(X) = c(20,50,100,1000)
  
  for (j in 1:4){
    mcmc = read.rcm.dmm(sprintf("%s_mcmc-%s.out", schemes[CA], names(X)[j]))
    
    D = xtabs(mcmc$dataset[,3] ~ mcmc$dataset[,1] + mcmc$dataset[,2])
    rownames(D) = macroevolution::tiplabels(mcmc$phy)[as.integer(rownames(D))+1L]
    colnames(D) = names(mcmc$dirichlet.prior)[as.integer(colnames(D))+1L]
    
    end = nrow(mcmc$pars)
    start = floor(end / 1.95)
    n = end-start+1
    
    w[j] = mean(rowSums(mcmc$pars[start:end,1:2]))
    
    post.states = mcmc$stateid[start:end,]
    
    X[[j]] = matrix(0, macroevolution::Ntip(mcmc$phy), ncat, #number of categories in the dataset
                    dimnames=list(macroevolution::tiplabels(mcmc$phy), names(mcmc$dirichlet.prior)))
    
    invisible(apply(post.states, 1, function(p) {
      up = unique(p)
      op = order(up)
      p = setNames(match(p, up[op]), names(p))
      x = t(sapply(lapply(split(D, p), matrix, ncol=ncat), colSums))+1
      x = sweep(x, 1, rowSums(x), "/")
      X[[j]] <<- X[[j]] + x[p,] / n
    }))
  }
  
  w = exp(w - max(w)) / sum(exp(w - max(w)))
  Y = matrix(0, macroevolution::Ntip(mcmc$phy), ncat, 
             dimnames=list(macroevolution::tiplabels(mcmc$phy), names(mcmc$dirichlet.prior)))
  
  for (k in 1:4){
    Y = Y + w[k] * X[[k]]
    Y = Y[,sort(colnames(Y))]
  }
  write.csv(Y, sprintf("diet-proportions-%s.csv", schemes[CA]))

}

## End script
