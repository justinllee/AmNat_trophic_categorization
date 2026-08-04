### Script objectives:
### Estimate diet tip rates of Elapids based on DPMS model analyses
### and based on diet breadth scores (plot for diet breadth is coded separately).
### Plot tip rates across phylogeny of Elapidae.
### Comparisons between diet breadth/proportions and their rates are also calculated.

## Prepare and install packages
library(remotes)
library(phylo) # remotes::install_github('blueraleigh/phylo')
library(bm) # remotes::install_github('blueraleigh/bm')
library(classInt)
library(ape)
library(scales)

## Prepare directories
wd <- "/PATH/"
setwd(wd)

## Import phylogeny
treefile <- '/PATH/Datasets/All_Elapids_tree_v2.tre' 
tree <- phylo::read.newick(treefile) #using phylo package
tree.a <- ape::read.tree(treefile) #using ape package

## Prepare diet category schemes and associated datasets
schemes <- c("prey_traditional", "prey_phylo1", "prey_phylo2", 
             "prey_functional1", "prey_functional2", "prey_type")
catnames <- c("traditional", "phylo1", "phylo2", 
              "functional1", "functional2", "type")
prop <- sprintf("/PATH/BHC_%s/diet-proportions-%s.csv", 
                schemes, schemes)
breadth <- sprintf("/PATH/BHC_%s/diet-breadth-%s.csv", 
                   schemes, schemes)

## Set up empty matrices to place tip rates
pr_ratemat <- matrix(0, Ntip(tree), length(prop), dimnames = list(tiplabels(tree), catnames))
br_ratemat <- matrix(0, Ntip(tree), length(prop), dimnames = list(tiplabels(tree), catnames))

##-----Part 1: Calculate tip rates ------

## Begin calculations (this portion of script adapted from Title et al. 2024)
for(TR in 1:length(prop)){
  
  # -------------------------
  ## Tip rates based on diet proportions from Markov Process model
  f10 = prop[TR]
  
  diet = data.matrix(read.csv(f10, row.names=1))
  subtree <- keep.tip(tree, rownames(diet))
  dietTR = matrix(0, Ntip(subtree), 4, dimnames=list(
    tiplabels(subtree), c("identity", "log", "alr", "clr")))
  alr = log(sweep(diet[,-ncol(diet)], 1, diet[,ncol(diet)], "/"))
  clr = sweep(log(diet), 1, rowMeans(log(diet)))
  dietTR[, 1] = sapply(bm.mvtiprate(diet, subtree), function(m) sum(diag(m)))
  dietTR[, 2] = sapply(bm.mvtiprate(log(diet), subtree), function(m) sum(diag(m)))
  dietTR[, 3] = sapply(bm.mvtiprate(alr, subtree), function(m) sum(diag(m)))
  dietTR[, 4] = sapply(bm.mvtiprate(clr, subtree), function(m) sum(diag(m)))
  
  # -----------------------------
  ## Tip rates based on diet breadth scores
  f11 <- breadth[TR]
  dietbreadth <- read.csv(f11, row.names = 1)
  dietbreadth <- setNames(dietbreadth[,3], dietbreadth[,1])
  subtree <- keep.tip(tree, names(dietbreadth))
  
  dietbreadthTR = matrix(0, Ntip(subtree), 2, dimnames=list(
    tiplabels(subtree), c("identity", "log")))
  
  dietbreadthTR[, 1] = bm.tiprate(dietbreadth, subtree)
  dietbreadthTR[, 2] = bm.tiprate(log(dietbreadth), subtree)
  
  # -----------------------------
  ## Write results to file
  taxa <- as.vector(rownames(diet))
  pr_ratemat[, TR] <- dietTR[, 4]
  br_ratemat[, TR] <- dietbreadthTR[, 2]
  rm(dietTR); rm(dietbreadthTR)
}

## Write results to .csv files
pr_ratemat <- as.data.frame(pr_ratemat)
outfile1 <- "Elapid_dietprop_tiprates.csv"
write.csv(pr_ratemat, file = outfile1, row.names = FALSE)

br_ratemat <- as.data.frame(br_ratemat)
outfile2 <- "Elapid_dietbreadth_tiprates.csv"
write.csv(br_ratemat, file = outfile2, row.names = FALSE)

##-----Part 2: Plot tip rates onto phylogeny ------

## Tip rates versus diet breadth for each category
breadths <- matrix(0, length(tree.a$tip.label), length(breadth), 
                   dimnames = list(tree.a$tip.label, catnames))

for(TR in 1:length(breadth)){
  f11 <- breadth[TR]
  db <- read.csv(f11, row.names = 1)
  breadths[,TR] <- db[,3]
}
breadths <- as.data.frame(breadths)

## log-transform rates prior to plotting
pr_ratemat.L <- log(pr_ratemat)
br_ratemat.L <- log(br_ratemat)

## First plot the diet proportion tip rates
quartz(width = 10, height = 8)
plot.phylo(tree.a, cex = 0.25, edge.width = 1.0)
pdf(file = "Fig4_Diet_Tip_Rates.pdf", width = 10, height = 8)
par(mar=c(0,1,0,0))
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv) 
#stores all the information of the tree that you just plotted

yy <- obj$yy[1:obj$Ntip] #subsets the number of nodes (yy) to the number of tips in tree
xx <- obj$xx[1:obj$Ntip]
#points(x = xx, y = yy, col = "red")
range(obj$xx)
plot.phylo(tree.a, x.lim = c(0, max(obj$xx) * 5), cex = 0.2, edge.width = 0.5) #changes the x-axis limits of plot space
#Begin loop
categories <- c(names(pr_ratemat.L))
colpalDiv <- colorRampPalette(rev(c('#d7191c', '#fdae61','#ffffbf', '#abd9e9', 'darkblue')), alpha = TRUE) 
origin <- 25 #hard codes the x-position so that every iteration adds a break
stretch <- 12.5 #hard codes x-position adjustment for subsequnt categories 
text(x = 10, y = 253.5, 'Diet tip rates (proportions)', pos = 1, cex = 0.75)
for(i in 1:length(categories)){
  traitdat <- pr_ratemat.L[,i]
  names(traitdat) <- tree.a$tip.label
  medianTraitDat <- median(traitdat, na.rm = TRUE)
  traitDiffs <- traitdat - medianTraitDat
  
  #tip rates for diet
  xvec <- c(0, max(obj$xx) * 5) #vector of the expanded x-axis
  yvec <- c(0, max(obj$yy)) #vector of the y-axis, should be 0 to 245 (the same no. as tips)
  xpos <- origin #position of tip rates on x-axis
  datRange <- range(traitDiffs)
  
  segments(x0 = xpos, x1 = xpos, y0 = 0, y1 = 245, xpd = NA, lend = 1, lwd = 1, lty = 1, col = "gray")
  text(x = xpos, y = 0, 'median', pos = 1, cex = 0.5)
  text(x = xpos, y = 252.5, categories[i], pos = 1, cex = 0.5)
  
  brks <- classIntervals(traitdat, n = 245, style = 'equal')
  cols <- findColours(brks, pal = colpalDiv(length(brks$brks) - 1))
  
  ref <- rescale(traitDiffs, from = datRange, to = datRange + xpos) #scales rates to median x-position
  vec <- c(1:length(tree.a$tip.label)) #vector so you can iterate through all the rates
  
  for(k in 1:length(tree.a$tip.label)){
    segments(x0 = xpos, x1 = ref[k], y0 = vec[k], y1 = vec[k], lwd = 0.5, xpd = NA, lend = 1, col = "gray")
    points(x = ref[k], y = vec[k], xpd = NA, pch = 21, cex = 0.4, bg = cols[k], col = "gray", lwd = 0.05)
  }
  origin <- c(origin) + stretch
}
dev.off()


##-----Part 3: Compare tip rates ------

## Plot diet breadth versus diet breadth rates first
quartz(width = 9.75, height = 6.5)
plot.new()
par(mfrow = c(2,3), 
    oma = c(2,2,2,2), 
    mar = c(3,3,3,2),
    cex.axis = 0.75,
    cex.lab = 0.75,
    mgp=c(1.5,0.5,0),
    xpd=TRUE)

plot(breadths$traditional, br_ratemat.L$traditional, 
     pch = 20, xlab = 'log(diet breadth)', 
     ylab = 'diet breadth tip rates', las = 1)
title(main = 'traditional', line = 0.75)

plot(breadths$phylo1, br_ratemat.L$phylo1, 
     pch = 20, xlab = 'log(diet breadth)', 
     ylab = 'diet breadth tip rates', las = 1)
title(main = 'phylo1', line = 0.75)

plot(breadths$phylo2, br_ratemat.L$phylo2, 
     pch = 20, xlab = 'log(diet breadth)', 
     ylab = 'diet breadth tip rates', las = 1)
title(main = 'phylo2', line = 0.75)

plot(breadths$functional1, br_ratemat.L$functional1, 
     pch = 20, xlab = 'log(diet breadth)', 
     ylab = 'diet breadth tip rates', las = 1)
title(main = 'functional', line = 0.75)

plot(breadths$functional2, br_ratemat.L$functional2, 
     pch = 20, xlab = 'log(diet breadth)', 
     ylab = 'diet breadth tip rates', las = 1)
title(main = 'habitat', line = 0.75)

plot(breadths$type, br_ratemat.L$type, 
     pch = 20, xlab = 'log(diet breadth)', 
     ylab = 'diet breadth tip rates', las = 1)
title(main = 'prey type', line = 0.75)

## Do the same between both rates
quartz(width = 9.75, height = 6.5)
plot.new()
par(mfrow = c(2,3), 
    oma = c(2,2,2,2), 
    mar = c(3,3,3,2),
    cex.axis = 0.75,
    cex.lab = 0.75,
    mgp=c(1.5,0.5,0),
    xpd=TRUE)

plot(br_ratemat.L$traditional, pr_ratemat.L$traditional,
     pch = 20, xlab = 'diet breadth tip rates', 
     ylab = 'diet prop. tip rates', las = 1)
title(main = 'taxonomy', line = 0.75)

plot(br_ratemat.L$phylo1, pr_ratemat.L$phylo1,
     pch = 20, xlab = 'log(diet breadth)', 
     ylab = 'diet breadth tip rates', las = 1)
title(main = 'phylo1', line = 0.75)

plot(br_ratemat.L$phylo2, pr_ratemat.L$phylo2,
     pch = 20, xlab = 'log(diet breadth)', 
     ylab = 'diet breadth tip rates', las = 1)
title(main = 'phylo2', line = 0.75)

plot(br_ratemat.L$functional1, pr_ratemat.L$functional1,
     pch = 20, xlab = 'diet breadth tip rates', 
     ylab = 'diet prop. tip rates', las = 1)
title(main = 'functional', line = 0.75)

plot(br_ratemat.L$functional2, pr_ratemat.L$functional2,
     pch = 20, xlab = 'diet breadth tip rates', 
     ylab = 'diet prop. tip rates', las = 1)
title(main = 'habitat', line = 0.75)

plot(br_ratemat.L$type, pr_ratemat.L$type,
     pch = 20, xlab = 'diet breadth tip rates', 
     ylab = 'diet prop. tip rates', las = 1)
title(main = 'prey type', line = 0.75)

## End script
