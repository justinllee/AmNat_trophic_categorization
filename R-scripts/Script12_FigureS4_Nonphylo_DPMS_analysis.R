### Script objectives:
### Non-Phylogenetic implementation of Dirichlet Proces Multinomial State model 
### of Elapid diet data, including diet breadth analysis

## Prepare and install packages
devtools::install_github('blueraleigh/bhc')
library(ape)
library(devtools)
library(bhc)
library(RColorBrewer)

## Prepare directories (use PATH folder from previous scripts)
wd <- "/PATH/Nonphylo BHC"
setwd("/PATH/Nonphylo BHC")
Elapids <- read.csv("/PATH/Datasets/All_Elapids_squamatabase.csv")
taxalist <- read.csv("/PATH/Datasets/Elapid_full_taxon_list.csv")

# NOTE: How you order the diet dataset matters!

##-----Part 1: Bayesian Hierarchical Clustering------

## Function to cut states into clusters based on BHC residuals
cutree_bhc = function(res, alpha=0.5) {
  n = nrow(res$merge)
  if (res$height[n] > res$height[n-1])
    root = res$merge[n,]
  else
    root = n
  tips = function(node) {
    ans = c()
    foo = function(node) {
      if (node > 0) {
        lf = res$merge[node, 1]
        rt = res$merge[node, 2]
        foo(lf)
        foo(rt)
      } else {
        ans <<- c(ans, node)
      }
    }
    foo(node)
    -ans
  }
  bar = function(node) {
    if (node > 0) {
      if (exp(res$height[node]) < alpha) {
        lf = res$merge[node, 1]
        rt = res$merge[node, 2]
        bar(lf)
        bar(rt)
      } else {
        g <<- g + 1
        ans[tips(node)] <<- g
      }
    } else {
      g <<- g + 1
      ans[tips(node)] <<- g
    }
  }
  ans = integer(n+1)
  g = 0
  for (node in root)
    bar(node)
  structure(ans, names=res$labels)
}

## Prepare category schemes
schemes <- c("prey_traditional", "prey_phylo1", "prey_phylo2", 
             "prey_functional1", "prey_functional2", "prey_type")

## Develop matrix of BHC states
data <- matrix(data = 0, nrow = length(unique(Elapids$predator)), ncol = length(schemes))
rownames(data) <- unique(Elapids$predator)
colnames(data) <- schemes

data1 <- matrix(data = 0, nrow = length(unique(Elapids$predator)), ncol = length(schemes))
rownames(data1) <- unique(Elapids$predator)
colnames(data1) <- schemes

## Perform Bayesian Hierarchical Clustering
for(n in 1:length(schemes)){
  sch <- schemes[n]
  
  #reads in BHC diet datasets
  snakediet <- read.csv(sprintf("bhc_dataset_%s.csv", sch))
  x = unclass(xtabs(snakediet$count ~ snakediet$species + snakediet$food))
  res = bhc.multinomial(x, alpha=1, beta=rep(1, ncol(x)))
  saveRDS(res, sprintf("bhc_%s.rds", sch))
  
  #use cutree_bhc() to estimate states
  bhc_states1 <- as.data.frame(cutree_bhc(res))
  row.names(bhc_states1[bhc_states1 == 1])
  colnames(bhc_states1) <- "states"
  
  #Input estimated states into matrix
  data[,sch] <- bhc_states1$states
}

## Export summary of diet states
write.csv(data, file = "nonphylo_diet_states_summary.csv")

##-----Part 2: Calculate Diet Breadth------

## Calculate diet proportions and diet breadth for each category
states <- as.data.frame(data)
breadth <- matrix(data = 0, nrow = length(unique(Elapids$predator)), ncol = length(schemes))
rownames(breadth) <- unique(Elapids$predator)
colnames(breadth) <- schemes

for(m in 1:length(schemes)){
  sch <- schemes[m]
  snakediet <- read.csv(sprintf("bhc_dataset_%s.csv", sch))
  x = xtabs(snakediet$count ~ snakediet$species + snakediet$food)
  fit = readRDS(sprintf("bhc_%s.rds", sch))
  diet.states = t(sapply(split(rownames(x), cutree_bhc(fit, 0.5)), 
                         function(d)
                         {
                           colSums(x[d, , drop=FALSE])+1
                         }
  ))
  diet.states = sweep(diet.states, 1, rowSums(diet.states), "/")
  write.csv(diet.states, file = sprintf("state_proportions_%s.csv", sch))
  y = structure(diet.states[cutree_bhc(fit, 0.5),], dimnames=list(rownames(x),colnames(diet.states)))
  diet.breadth = apply(y, 1, function(p) 1 / sum(p*p))
  breadth[,sch] <- diet.breadth
  write.csv(diet.breadth, file = sprintf("diet_breadth_%s.csv", sch))
}

breadth <- as.data.frame(breadth)
index3 <- match(rownames(breadth), taxalist$taxon)
breadth$radiation <- taxalist$geog_radiation[index3]
breadth[1:6] <- log(breadth[1:6])

## ANOVA of diet breadth scores
bd <- breadth
anovasummary <- vector("list", 6)
tukeyHSD <- vector("list", 6)
for(an in 1:length(schemes)){
  sch <- schemes[an]
  anovaresult <- aov(as.matrix(log(breadth[sch])) ~ breadth$radiation)
  anovasummary[an] <- summary(anovaresult)
  tukeyHSD[an] <- TukeyHSD(anovaresult, conf.level = 0.95)
}
names(tukeyHSD) <- schemes
capture.output(tukeyHSD, file="TukeyHSD_diet_breadth_nonphylo.txt")

## Rename columns and export diet breadth scores
breadth <- breadth[c("radiation", schemes)]
scheme.lab <- c("[T1]", "[T2]", "[T3]", "[E1]", "[E2]", "[E3]")
colnames(breadth)[2:7] <- scheme.lab
write.csv(breadth, file = "diet_breadth_nonphylo_summary.csv", row.names = F)
# 'diet_breadth_nonphylo_summary.csv' will be used to calculate diet breadth

##-----Part 3: Plot diet niche states ------

## Import Elapid tree
tree <- read.tree("/PATH/Datasets/All_Elapids_tree_v2.tre")

## Prune tree to removed unsampled taxa
diff <- setdiff(tree$tip.label, rownames(states))
tree <- ape::drop.tip(tree, diff)
bad <- setdiff(rownames(states), tree$tip.label)
states <- states[rownames(states) %in% tree$tip.label,]
ord <- match(tree$tip.label, rownames(states))
states <- states[ord,]

## Combine dietary states with tip labels
pdf(width = 10, height = 10, file = "FigS4_DPMS_Nonphylo_States.pdf")
par(mar = c(1,1,1,1),
    oma = c(0,0,0,0),
    mgp=c(0,0,0),
    xpd=TRUE)
plot.phylo(tree, cex = 0.5, type = "fan", label.offset = 2.0)
brewer <- vector("list", 6)
ofst <- 0
shapes <- c(21, 22, 23, 24, 25, 8)

for(CA in 1:length(schemes)){
  sch <- schemes[CA]
  infstates <- states[,sch]
  names(infstates) <- rownames(states)
  nstates <- length(unique(infstates))
  
  #set colors
  nb.cols <- nstates
  brew <- colorRampPalette(brewer.pal(nstates, "Paired"))(nb.cols)
  brewer[[ CA ]] <- brew
  names(brew) <- unique(infstates)
  col_index <- match(infstates, names(brew))
  color_vector <- brew[col_index]
  
  #combine states with tip labels
  ape::tiplabels(pch = shapes[CA], bg = color_vector, col = color_vector,
                 cex = 0.5, offset = ofst)
  ofst <- c(ofst) + 0.35
}
dev.off()

## End script
