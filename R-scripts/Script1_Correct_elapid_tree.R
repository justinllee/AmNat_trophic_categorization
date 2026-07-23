### Script objectives:
### This script makes enhancements based on previous relationships in Elapidae.
### The base phylogeny is from Title et al. (2024).
### Refer to materials and methods section for further explanations on changes.

## Prepare packages
library(ape)
library(phytools)

## Prepare directories
path <- "/PATH/"
path <- setwd(path)

## Prepare data
setwd("/PATH/Datasets")
tree <- read.tree("best_ultrametric_fulltree_ddBD_revision.tre")
diet <- read.csv("All_Elapids_squamatabase_resubmission1.csv")

## Subset Elapids from Title et al. (2024) tree
span <- c("Calliophis_intestinalis", "Hydrophis_platurus")
mrca <- getMRCA(tree, tip = span)
etree <- ape::extract.clade(tree, node = mrca)

## Remove poorly placed species from tree
bad <- c("Hydrophis_hardwickii", "Pseudonaja_guttata", "Ogmodon_vitianus")
good <- ape::drop.tip(etree, bad)

## Change Parasuta to Suta
suta <- good$tip.label[63:70]
change <- c("Suta_nigriceps", "Echiopsis_curta", 
            "Suta_spectabilis", "Suta_gouldii", "Suta_suta",
            "Suta_fasciata", "Suta_punctata", "Suta_monachus")
good$tip.label[63:70] <- change

## Change Simoselaps bimaculatus to Narophis
narophis <- "Narophis_bimaculatus"
good$tip.label[109] <- narophis

## Case 1: Add Micrurus distans to the coral snake clade
good$edge
distans_edge <- good$edge.length[471] + good$edge.length[472]
good <- phytools::bind.tip(tree = good,
                           tip.label = "Micrurus distans",
                           edge.length = distans_edge + 0.75,
                           where = 460, 
                           position = 0.45)

## Case 2: Fix position of Pseudonaja guttata
pseudonaja <- c("Pseudonaja_modesta" , "Pseudonaja_nuchalis")
pguttata <- getMRCA(good, tip = pseudonaja)

good$edge #tip for Pseudonaja modesta is 98
p_edge <- good$edge.length[206]
good <- phytools::bind.tip(tree = good,
                           tip.label = "Pseudonaja_guttata",
                           edge.length = p_edge + 0.25,
                           where = pguttata,
                           position = 0.25)

## Case 3: Add Salomonelaps and Loveridgelaps to clade with Micropechis
hydro <- c("Micropechis_ikaheca" , "Hydrophis_platurus")
hymrca <- getMRCA(good, tip = hydro)

s_edge <- good$edge.length[284]
good1 <- phytools::bind.tip(tree = good,
                            tip.label = "Salomonelaps_par",
                            edge.length = s_edge + 0.1,
                            where = hymrca, position = 0.1)
quartz(width = 5, height = 10)
par(mar=c(0,1,0,0))
plot(good1, cex = 0.3)

which(good1$tip.label == "Salomonelaps_par")
good1 <- phytools::bind.tip(tree = good1,
                            tip.label = "Loveridgelaps_elapoides",
                            edge.length = 7.5, 
                            where = 141,
                            position = 7.5)

## Case 4: Fix position of Ogmodon vitianus
ogm <- c("Loveridgelaps_elapoides" , "Salomonelaps_par")
ogmrca <-getMRCA(good1, tip = ogm)
o_edge <- good$edge.length[391]
good1 <- phytools::bind.tip(tree = good1,
                            tip.label = "Ogmodon_vitianus", 
                            edge.length = 7.5 + good$edge.length[391] * 0.5, 
                            where = ogmrca,
                            position = good$edge.length[391] * 0.5)

## Case 5: Drop Sinomicrurus hatori based on Smart et al. (2021)
good1 <- drop.tip(good1, tip = "Sinomicrurus_hatori")

## Summarize changes in new tree
write.tree(good1, file = "All_Elapids_tree_v2.tre")

## End script
