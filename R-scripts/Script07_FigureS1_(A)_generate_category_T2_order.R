### Script objectives:
### Generate phylogenetic [T2] category scheme based on order-rank
### data from Timetree of life (Kumar et al. 2022) and existing elapid diet data

## Prepare packages
library(ape)
library(phytools)

## Prepare directories
wd <- "/PATH/phylo category"
setwd(wd)

## Prepare data
data <- "/PATH/Datasets/All_Elapids_squamatabase.csv"
tree1 <- read.tree("multicellular_animals_order.nwk")

## Clean-up diet data
sqbase <- read.csv(data)
sqbase$prey <- chartr(" ", "_", sqbase$prey)
sqbase$prey_taxon <- chartr(" ", "_", sqbase$prey_taxon)

## Subset taxon column (this provides a rank-based classification scheme of each item)
prey <- unique(sqbase$prey_taxon)
orders <- tree1$tip.label

##-----Part 1: Clean-up and generate categories using Timetree of Life ------
string <- "\\W[0-9]+" #regex to capture hyphen + clade number for paraphyletic orders
only1 <- "\\W[1]$" #regex to capture the FIRST mention of a paraphyletic clade
bad <- NULL
good <- NULL
for(f in 1:length(orders)){
  tip <- orders[f]
  tmp <- grepl(string, tip)
  except <- grepl(only1, tip)
  if(tmp == TRUE & except == FALSE){
    bad <- c(bad, tip)
  } else {
    good <- c(good, tip)
  }
}

tree1 <- drop.tip(tree1, bad)
tree1$tip.label <- sub(only1 , "" , tree1$tip.label) #cleans up tip labels
clean_orders <- unique(tree1$tip.label)

## Match the tree tips with diet data at family-rank:
phcat <- NULL
prey <- strsplit(prey, split = ";") #creates big list of prey based on ranks
for(i in 1:length(prey)){
  tmp1 <- prey[[ i ]]
  if(identical(intersect(tmp1, clean_orders), character(0)) == FALSE){
    phcat[ i ] <- intersect(tmp1, clean_orders)
  } else {
    phcat[ i ] <- "unknown"
  }
}
phylo <- unique(phcat)

## Check food items that are missing from tree
missing <- prey[ phcat == 'unknown' ]
# In this case, only a few food items are missing: Chilopoda (Scolopendromorpha).
# Disregard the Bivalves (which are probably secondary prey ingestion)
phylo <- c(phylo, "Scolopendromorpha")

## Now subset tips on timetree of life
missing_tips <- setdiff(tree1$tip.label, phylo)
tree1 <- drop.tip(tree1, missing_tips)

## Add back major squamate suborders (dates based on Timetree of life)
library(phytools)
squam <- which(tree1$tip.label == "Squamata")
tree2 <- phytools::bind.tip(tree = tree1, tip.label = "Gekkota", 
                            edge.length = 186.7, where = squam, 
                            position = 186.7)
tree2 <- phytools::bind.tip(tree = tree2, tip.label = "Scincomorpha", 
                            edge.length = 173.5, where = squam, 
                            position = 173.5)
tree2 <- phytools::bind.tip(tree = tree2, tip.label = "Lacertoidea", 
                            edge.length = 166.9, where = squam, 
                            position = 166.9)
tree2 <- phytools::bind.tip(tree = tree2, tip.label = "Amphisbaenia", 
                            edge.length = 119.4, 
                            where = which(tree2$tip.label == "Lacertoidea"), 
                            position = 119.4)
tree2 <- phytools::bind.tip(tree = tree2, tip.label = "Iguania", 
                            edge.length = 161.0, 
                            where = squam, position = 161.0)
tree2 <- phytools::bind.tip(tree = tree2, tip.label = "Anguiomorpha", 
                            edge.length = 157.2, 
                            where = which(tree2$tip.label == "Iguania"), 
                            position = 157.2)
tree2$tip.label[squam] <- "Serpentes"

## Subset timetree to exclude internal nodes
nTip <- length(tree2$tip.label) #56 tips
nEdge <- length(tree2$edge.length) #104 nodes
tip_lengths <- tree2$edge.length[tree2$edge[,2] <= nTip]

##-----Part 2: Filter tips are older than the K-Pg mass extinction event ------
old <- tree2$tip.label[tip_lengths >= 66.043]
new <- tree2$tip.label[tip_lengths <= 66.043]
write.tree(tree3, file = "phylo_category_order.tre")

## Plot to check old and new tips
plot.new()
par(mar=c(0,1,0,0))
plot(tree2, cex = 0.35)
par(mar=c(0,1,0,0))
tiplabels(tip = which(tip_lengths <= 66.043),  pch = 20, col = "red", cex = 0.5)
dev.off()

## Drop and rename the tips that emerged after the K-Pg extinction
todrop <- c("Bucerotiformes", "Coraciiformes", "Peramelemorphia")
tree3 <- drop.tip(tree2, todrop)
tree3$tip.label[which(tree3$tip.label == "Dasyuromorphia")] <- "Marsupialia"
tree3$tip.label[which(tree3$tip.label == "Piciformes")] <- "Picocoraciae"

## Plot final version of tree
quartz(height = 5, width = 5)
plot.new()
pdf(file = "phylo_category_order_tree2.pdf", height = 7.5, width = 5.0)
par(mar=c(0,0,0,0))
plot(tree3, cex = 0.35)
dev.off()

## File 'phylo_category_order_tree2.pdf' was combined into Figure S1.
## All order-level prey items were added to 'All_Elapids_squamatabase.csv'.

## End script
