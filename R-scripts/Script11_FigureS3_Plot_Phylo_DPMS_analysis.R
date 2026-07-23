### Script objectives:
### Plot results from phylogenetic implementation of
### Dirichlet Process Multinomial State model of Elapid diet data

## Prepare and install packages
library(ape)
library(RColorBrewer)
library(phytools)

## Prepare directories, different method used here so that loop works
## (use the same PATH folder from previous scripts)
path <- "/PATH/Phylo BHC/"
path <- setwd(path)

## Import phylogeny
tree = read.tree("/PATH/Datasets/All_Elapids_tree_v2.tre")

## Prepare category schemes and associated datasets
schemes <- c("prey_traditional", "prey_phylo1", "prey_phylo2", 
             "prey_functional1", "prey_functional2", "prey_type")
datafiles <- sprintf("bhc_dataset_%s.csv", schemes)
datafiles <- datafiles[1:6]

## Plot
setwd(path)
datafiles <- sprintf("%s_spp_dietstates_245_max_prob.csv", schemes)
pdf(width = 10, height = 10, file = "FigS3_DPMS_Phylo_States.pdf")
par(mar = c(1,1,1,1),
    oma = c(0,0,0,0),
    mgp=c(0,0,0),
    xpd=TRUE)
plot.phylo(tree, cex = 0.5, type = "fan", label.offset = 2.0)
brewer <- vector("list", 6)
ofst <- 0
shapes <- c(21, 22, 23, 24, 25, 8)

for(CA in 1:length(schemes)){
  infstates <- read.csv(sprintf("BHC_%s/%s_spp_dietstates_245_max_prob.csv", schemes[ CA ],
                                schemes[ CA ]))
  nstates <- length(unique(infstates$V2))
  
  #set colors
  nb.cols <- nstates
  brew <- colorRampPalette(brewer.pal(nstates, "Paired"))(nb.cols)
  brewer[[ CA ]] <- brew
  names(brew) <- unique(infstates$V2)
  match(infstates$V2, names(brew))
  color_vector <- brew[infstates$V2]
  
  #combine states with tip labels
  ape::tiplabels(pch = shapes[CA], bg = color_vector, col = color_vector,
                 cex = 0.5, offset = ofst)
  ofst <- c(ofst) + 0.35
}

dev.off()