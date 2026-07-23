### Script objectives:
### Plot tip rates for diet breadth across phylogeny of Elapidae.

## Prepare and install packages
library(remotes)
library(phylo) # remotes::install_github('blueraleigh/phylo')
library(bm) # remotes::install_github('blueraleigh/bm')
library(classInt)
library(ape)
library(scales)

## Prepare directories (use same PATH from previous scripts)
wd <- "/PATH/Phylo BHC"
setwd(wd)

## Import phylogeny
treefile <- '/PATH/Datasets/All_Elapids_tree_v2.tre' 
tree <- phylo::read.newick(treefile) #using phylo package
tree.a <- ape::read.tree(treefile) #using ape package

## Import data
br_ratemat <- read.csv("Elapid_dietbreadth_tiprates.csv")

##----- Plot tip rates onto phylogeny ------

## log-transform rates prior to plotting
br_ratemat.L <- log(br_ratemat)

## Plot
plot.phylo(tree.a, cex = 0.25, edge.width = 0.5)
pdf(file = "FigS6_Diet_Tip_Rates_Breadth", width = 10, height = 8)
par(mar = c(0,1,0,0))
obj<-get("last_plot.phylo",envir=.PlotPhyloEnv) 
yy <- obj$xx[1:obj$Ntip]
xx <- obj$yy[1:obj$Ntip]

range(obj$xx)
plot.phylo(tree.a, x.lim = c(0, max(obj$xx) * 6), cex = 0.2, edge.width = 0.5)

categories <- c(names(br_ratemat.L))
colpalDiv <- colorRampPalette(rev(c('#d7191c', '#fdae61','#ffffbf', '#abd9e9', 'darkblue')), alpha = TRUE)
origin <- 30
stretch <- 15
text(x = 10, y = 252.5, 'Diet tip rates (breadth)', pos = 1, cex = 0.75)
for(j in 1:length(categories)){
  breadth_traitdat <- br_ratemat.L[,j]
  names(breadth_traitdat) <- tree.a$tip.label
  median_brTraitDat <- median(breadth_traitdat, na.rm = TRUE)
  br_traitDiffs <- breadth_traitdat - median_brTraitDat
  
  br_xvec <- c(0, max(obj$xx) * 3) #vector of the expanded x-axis
  br_yvec <- c(0, max(obj$yy)) #vector of the y-axis, should be 0 to 245 (the same no. as tips)
  br_xpos <- origin #position of tip rates on x-axis
  br_datRange <- range(br_traitDiffs)
  
  segments(x0 = br_xpos, x1 = br_xpos, y0 = 0, y1 = 245, xpd = NA, lend = 1, lwd = 1, lty = 1, col = "gray")
  text(x = br_xpos, y = 0, 'median', pos = 1, cex = 0.5)
  text(x = br_xpos, y = 252.5, categories[j], pos = 1, cex = 0.5)
  
  br_brks <- classIntervals(breadth_traitdat, n = 245, style = 'equal')
  br_cols <- findColours(br_brks, pal = colpalDiv(length(br_brks$brks) - 1))
  
  br_ref <- rescale(br_traitDiffs, from = br_datRange, to = br_datRange + br_xpos) #scales rates to median x-position
  br_vec <- c(1:length(tree.a$tip.label)) #vector so you can iterate through all the rates
  
  for(l in 1:length(tree.a$tip.label)){
    segments(x0 = br_xpos, x1 = br_ref[l], y0 = br_vec[l], y1 = br_vec[l], lwd = 0.5, xpd = NA, lend = 1, col = "gray")
    points(x = br_ref[l], y = br_vec[l], xpd = NA, pch = 21, cex = 0.4, bg = br_cols[l], col = "gray", lwd = 0.05)
  }
  origin <- c(origin) + stretch
}
dev.off()

## End script
