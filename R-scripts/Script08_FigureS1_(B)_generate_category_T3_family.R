### Script objectives:
### Generate phylogenetic [T3] category scheme based on family-rank
### data from Timetree of life (Kumar et al. 2022) and existing elapid diet data

## Prepare packages and directory
library(ape)
library(phytools)

## Prepare directories
wd <- "/PATH/"
setwd(wd)

## Prepare data
data <- "/PATH/Datasets/All_Elapids_squamatabase.csv"
tree1 <- read.tree("multicellular_animals_family.nwk")

## Clean-up diet data
sqbase <- read.csv(data)
sqbase$prey <- chartr(" ", "_", sqbase$prey)
sqbase$prey_taxon <- chartr(" ", "_", sqbase$prey_taxon)

## Subset taxon column (this provides a rank-based classification scheme of each item)
prey <- unique(sqbase$prey_taxon)
families <- tree1$tip.label

##-----Part 1: Clean-up and generate categories using Timetree of Life ------
string <- "\\W[0-9]+" #regex to capture hyphen + clade number for paraphyletic orders
only1 <- "\\W[1]$" #regex to capture the FIRST mention of a paraphyletic clade (e.g., Citharinidae-1)
bad <- NULL
good <- NULL
for(f in 1:length(families)){
  tip <- families[f]
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
clean_families <- unique(tree1$tip.label)

## Match the tree tips with diet data at family-rank:
phcat <- NULL
prey <- strsplit(prey, split = ";") #creates list of food items based on rank
for(i in 1:length(prey)){
  tmp1 <- prey[[ i ]]
  if(identical(intersect(tmp1, clean_families), character(0)) == FALSE){
    phcat[ i ] <- intersect(tmp1, clean_families)
  } else {
    phcat[ i ] <- "unknown"
  }
}
phylo <- unique(phcat)

## Check food items that are missing from tree
missing <- prey[ phcat == 'unknown' ] 

# We see that the following food items in the dataset higher than family-rank
# are present: Decapoda, Orthoptera, Blattodea, Scorpiones.
# We can add them below:
add <- c("Cancridae", "Gryllidae", "Blattidae", "Buthidae") #use to represent orders
phylo <- c(phylo, add)

## Now subset tips on timetree of life
missing_tips <- setdiff(tree1$tip.label, phylo)
tree1 <- drop.tip(tree1, missing_tips)
change <- match(add, tree1$tip.label)
order <- c("Decapoda", "Orthoptera", "Blattodea", "Scorpiones")
tree1$tip.label[change] <- order
plot(tree1, cex = 0.5)

## Subset timetree to exclude internal nodes
Ntip <- length(tree1$tip.label) #138 tips
Nnode <- tree1$edge.length #274 nodes
tip_lengths <- tree1$edge.length[tree1$edge[,2] <= Ntip]

##-----Part 2: Filter tips are older than the K-Pg mass extinction event ------
old <- tree1$tip.label[tip_lengths >= 66.043]
new <- tree1$tip.label[tip_lengths <= 66.043]
write.tree(tree1, file = "phylo_category_family.tre")

## Plot to check old and new tips
plot.new()
par(mar=c(0,1,0,0))
plot(tree1, cex = 0.35)
par(mar=c(0,1,0,0))
tiplabels(tip = which(tip_lengths <= 66.043),  pch = 20, col = "red", cex = 0.5)
dev.off()

## Find tips in tree to keep
tokeep <- c("Leiognathidae", "Lutjanidae", "Ariidae", "Dasyuridae", "Herpestidae", 
            "Pteropodidae", "Nesomyidae", "Pygopodidae", "Pythonidae", "Colubridae", 
            "Columbidae", "Picidae", "Turdidae")
todrop <- setdiff(new, tokeep)
tree2 <- drop.tip(tree1, todrop)
kept <- NULL; found <- NULL
for(find in 1:length(tokeep)){
  kept <- tokeep[find]
  found[find] <- which(tree2$tip.label == kept) 
}
names(found) <- tokeep

## Give tips new names for final tree
newnames <- c("Eupercaria-is-1", "Eupercaria-is-2", "Siluroidei", "Marsupialia", "Carnivora",
              "Chiroptera", "Muroidea", "Pygopodoidea", "Constrictores", "Colubriformes", 
              "Columbimorphae", "Picodynastornithes", "Passeriformes")
for(find in 1:length(tokeep)){
  replace <- found[find]
  newname <- newnames[find]
  tree2$tip.label[replace] <- newname
}
write.tree(tree2, file = "phylo_category_family.tre")

## Plot final version of tree
plot.new()
pdf(file = "phylo_category_family_tree.pdf", height = 7.5, width = 5.0)
par(mar=c(0,0,0,0))
plot(tree2, cex = 0.35)
dev.off()

## File 'phylo_category_family_tree.pdf' was combined into Figure S1.
## All family-level prey items were added to 'All_Elapids_squamatabase.csv'.

## End script
