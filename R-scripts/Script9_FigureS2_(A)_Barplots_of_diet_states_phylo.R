### Script objectives:
### Visualize dietary niche states estimated from phylogenetic implementation
### of Dirichlet Process Multinomial State (DPMS) model

## Prepare directories
wd <- "/PATH/"
setwd(wd)

## Prepare packages
library(classInt)
library(ape)
library(scales)

##-----Part 1: Generate barplots for Phylogenetic Results ------

## Traditional
trad_states <- read.csv("/PATH/BHC_prey_traditional/prey_traditional_inferred_states.csv")
trad_states <- trad_states[1:9,] #Omit bottommost row
trad_states <- t(trad_states)
trad_states <- apply(trad_states, 2, function(x){x*100/sum(x,na.rm=T)}) #converts to percent

state_order <- c("amphibians", "frogs", "caecilians", "lizards", "reptile_eggs", 
                 "amphisbaenians", "snakes", "birds", "bird_eggs", "mammals", 
                 "fish", "invertebrates", "other")
names(state_order) <- c("indianred2", "indianred", "indianred4", "deepskyblue", 
                        "turquoise", "cyan1", "seagreen3", "darkolivegreen2",
                        "darkolivegreen", "royalblue3", "mediumpurple", 
                        "goldenrod1", "gray")
trad_states <- trad_states[state_order,]

## Phylogenetic (order)
phylo1 <- read.tree(file = "/PATH/phylo category/phylo_category_order.tre")
fish <- phylo1$tip.label[1:18]
amphibian <- phylo1$tip.label[19:21]
mammals <- phylo1$tip.label[22:28]
serpentes <- phylo1$tip.label[29]
squamates <- phylo1$tip.label[30:35]
turtles <- phylo1$tip.label[42]
birds <- phylo1$tip.label[36:41]
velvet <- phylo1$tip.label[43]
arthropods <- phylo1$tip.label[44:51]
worms <- phylo1$tip.label[52]
mollusk <- phylo1$tip.label[53]

actinopterygii <- colorRampPalette(c("plum1", "mediumpurple"))
actinopterygii <- actinopterygii(19)
names(actinopterygii) <- c(fish, "Actinopterygii")

amphibia <- c("indianred4", "indianred2", "indianred")
names(amphibia) <- amphibian

mammalia <- colorRampPalette(c("royalblue2", "royalblue4"))
mammalia <- mammalia(10)
names(mammalia) <- c(mammals,"Dasyuromorphia", "Peramelemorphia", "Mammalia")

squamata <- colorRampPalette(c("cyan1", "deepskyblue"))
squamata <- squamata(8)
names(squamata) <- c(squamates, "Squamata", "Reptilia")

serpentes <- "seagreen1"
names(serpentes) <- "Serpentes"

testudines <- "darkgreen"
names(testudines) <- turtles

aves <- colorRampPalette(c("darkolivegreen2", "springgreen"))
aves <- aves(9)
names(aves) <- c(birds, "Coraciiformes", "Piciformes", "Aves")

onychophora <- "lightgoldenrod2"
names(onychophora) <- "Euonychophora"
  
arthropoda <- colorRampPalette(c("goldenrod4", "goldenrod1"))
arthropoda <- arthropoda(10)
names(arthropoda) <- c(arthropods, "Insecta", "Chilopoda")

annelid <- "gold1"
names(annelid) <- worms

mollusca <- colorRampPalette(c("slategray1", "slategray4"))
mollusca <- mollusca(3)
names(mollusca) <- c(mollusk, "Mollusca", "Cephalopoda")

other <- c("gray", "gray20")
names(other) <- c("unknown", "Chordata")

phylo1colors <- c(actinopterygii, amphibia, mammalia, 
                  squamata, serpentes, testudines, aves, onychophora, 
                  arthropoda, annelid, mollusca, other)

phylo1_states <- read.csv("/PATH/BHC_prey_phylo1/prey_phylo1_inferred_states.csv")
phylo1_states <- phylo1_states[1:6,]
phylo1_states <- t(phylo1_states)
phylo1_states <- apply(phylo1_states, 2, function(x){x*100/sum(x,na.rm=T)}) #converts to percent

phylo1colors <- phylo1colors[names(phylo1colors) %in% rownames(phylo1_states)]
phylo1_states <- phylo1_states[names(phylo1colors),]

## Phylogenetic (family)
phylo2 <- read.tree(file = "/PATH/phylo category/phylo_category_family.tre")
fish <- phylo2$tip.label[1:51]
fish[1] <- "Eupercaria.is.1"
fish[3] <- "Eupercaria.is.2"
amphibian <- phylo2$tip.label[52:70]
mammals <- phylo2$tip.label[71:79]
snakes <- phylo2$tip.label[88:94]
squamates <- phylo2$tip.label[c(80:87, 95:106)]
turtles <- phylo2$tip.label[113]
birds <- phylo2$tip.label[107:112]
velvet <- phylo2$tip.label[114]
arthropods <- phylo2$tip.label[115:121]
worms <- phylo2$tip.label[122]
mollusk <- phylo2$tip.label[123]

actinopterygii <- colorRampPalette(c("plum1", "mediumpurple"))
actinopterygii <- actinopterygii(54)
names(actinopterygii) <- c(fish, "Actinopterygii", "Anguilliformes", "Perciformes")

amphibia <- colorRampPalette(c("indianred4", "indianred"))
amphibia <- amphibia(22)
names(amphibia) <- c(amphibian, "Anura", "Gymnophiona")

mammalia <- colorRampPalette(c("royalblue2", "royalblue4"))
mammalia <- mammalia(11)
names(mammalia) <- c(mammals, "Mammalia", "Rodentia")

squamata <- colorRampPalette(c("cyan1", "deepskyblue"))
squamata <- squamata(23)
names(squamata) <- c(squamates, "Squamata", "Reptilia", "Amphisbaenia")

serpentes <- colorRampPalette(c("seagreen1", "seagreen4"))
serpentes <- serpentes(8)
names(serpentes) <- c(snakes, "Serpentes")

testudines <- "darkgreen"
names(testudines) <- turtles

aves <- colorRampPalette(c("darkolivegreen2", "darkolivegreen4"))
aves <- aves(7)
names(aves) <- c(birds, "Aves")

onychophora <- "lightgoldenrod2"
names(onychophora) <- "Euonychophora"

arthropoda <- colorRampPalette(c("goldenrod4", "goldenrod1"))
arthropoda <- arthropoda(9)
names(arthropoda) <- c(arthropods, "Insecta", "Chilopoda")

annelid <- colorRampPalette(c("gold1", "gold"))
annelid <- annelid(2)
names(annelid) <- c(worms, "Haplotaxida")

mollusca <- colorRampPalette(c("slategray1", "slategray4"))
mollusca <- mollusca(3)
names(mollusca) <- c(mollusk, "Mollusca", "Cephalopoda")

other <- c("gray", "gray20")
names(other) <- c("unknown", "Chordata")

phylo2colors <- c(actinopterygii, amphibia, mammalia, 
                  squamata, serpentes, testudines, aves, 
                  onychophora, arthropoda, annelid, mollusca, other)

phylo2_states <- read.csv("/PATH/BHC_prey_phylo2/prey_phylo2_inferred_states.csv")
phylo2_states <- phylo2_states[1:6,]
phylo2_states <- t(phylo2_states)
phylo2_states <- apply(phylo2_states, 2, function(x){x*100/sum(x,na.rm=T)}) #converts to percent

phylo2colors <- phylo2colors[names(phylo2colors) %in% rownames(phylo2_states)]
phylo2_states <- phylo2_states[names(phylo2colors),]

## Functional
functional_states <- read.csv("/PATH/BHC_prey_functional1/prey_functional1_inferred_states.csv")
functional_states <- functional_states[1:7,]
functional_states <- t(functional_states)
functional_states <- apply(functional_states, 2, function(x){x*100/sum(x,na.rm=T)}) #converts to percent

functional_order <- c("slippery", "scaly", "armor.durophagy", "elongate.scaly", 
                      "amniote.eggs", "fish.eggs", "elongate.slippery", 
                      "hairy.bulky", "winged.feathered", "soft.bodied", "chitinous")
names(functional_order) <- c("indianred2", "deepskyblue", "turquoise", 
                              "seagreen3", "darkolivegreen2", "royalblue3",
                             "indianred4", "mediumpurple",
                             "darkolivegreen", "sienna3", "goldenrod1")
functional_states <- functional_states[functional_order,]

## Habitat
habitat_states <- read.csv("/PATH/BHC_prey_functional2/prey_functional2_inferred_states.csv")
habitat_states <- habitat_states[1:9,]
habitat_states <- t(habitat_states)
habitat_states <- apply(habitat_states, 2, function(x){x*100/sum(x,na.rm=T)}) #converts to percent

habitat_order <- c("fossorial", "terrestrial", "aquatic", "marine", "arboreal")
names(habitat_order) <- c("indianred2", "deepskyblue", "sienna3",
                          "seagreen3", "mediumpurple")
habitat_states <- habitat_states[habitat_order,]

## MBT Prey Type
mbt_states <- read.csv("/PATH/BHC_prey_type/prey_type_inferred_states.csv")
mbt_states <- mbt_states[1:8,]
mbt_states <- t(mbt_states)
mbt_states <- apply(mbt_states, 2, function(x){x*100/sum(x,na.rm=T)}) #converts to percent

mbt_order <- c("X1", "X2", "X3", "X4")
names(mbt_order) <- c("indianred2", "deepskyblue", "sienna3", "seagreen3")
mbt_states <- mbt_states[mbt_order,]


##-----Part 2: Combine barplots into one figure ------
## Adjust rownames of states for Legend
rownames(trad_states) <- c("Amphibians", "Frogs", "Caecilians", "Lizards", 
                           "Reptile Eggs", "Amphisbaenians", "Snakes", "Birds", 
                           "Bird Eggs", "Mammals", "Fish", "Invertebrates", "Other")
rownames(functional_states) <- c("Slippery", "Scaly", "Armor Durophagy", 
                                 "Elongate Scaly", "Amniote Eggs", "Fish Eggs", 
                                 "Elongate Slippery", "Hairy Bulky", 
                                 "Winged Feathered", "Soft Bodied", "Chitinous")
rownames(habitat_states) <- c("Fossorial", "Terrestrial", "Aquatic", 
                              "Marine", "Arboreal" )
rownames(mbt_states) <- c("Type 1", "Type 2", "Type 3", "Type 4")

## Plot
plot.new()
pdf(height = 9, width = 12, file = "Phylo_BHC_dietstates_figure.pdf")
par(mfrow = c(2,3), 
    oma = c(1,1,1,1), 
    mar = c(3,6,6,6),
    cex.axis = 0.75,
    cex.lab = 0.75,
    mgp=c(1.5,0.5,0),
    xpd=TRUE)
barplot(as.matrix(trad_states), col = names(state_order), 
        border = NA, cex.axis = 1.25, cex.names = 1.0, cex.lab = 1.25,
        ylab = "Proportion of Diet (%)", mgp=c(2.25,0.5,0), 
        ylim = c(0,100),
        axes = TRUE, las = 1,
        legend.text = rownames(trad_states),
        args.legend = list(bty = "n", inset = -0.35, 
                           x = "right", ncol = 1, cex = 0.7, border = NA))
title(expression(paste("Traditional ", bold("[T1]"))), cex = 1.5, line = 0.75)

barplot(as.matrix(phylo1_states), col = phylo1colors, 
        border = NA, cex.axis = 1.25, cex.names = 1.0, cex.lab = 1.25,
        ylab = "Proportion of Diet (%)", mgp=c(2.25,0.5,0), 
        ylim = c(0,100),
        axes = TRUE, las = 1,)
title(expression(paste("Order ", bold("[T2]"))), cex = 1.5, line = 0.75)

barplot(as.matrix(phylo2_states), col = phylo2colors, 
        border = NA, cex.axis = 1.25, cex.names = 1.0, cex.lab = 1.25,
        ylab = "Proportion of Diet (%)", mgp=c(2,0.5,0), 
        ylim = c(0,100),
        axes = TRUE, las = 1)
title(expression(paste("Family ", bold("[T3]"))), cex = 1.5, line = 0.75)

barplot(as.matrix(functional_states), col = names(functional_order), 
        border = NA, cex.axis = 1.25, cex.names = 1.0, cex.lab = 1.25,
        ylab = "Proportion of Diet (%)", mgp=c(2.25,0.5,0), 
        ylim = c(0,100),
        axes = TRUE, las = 1, 
        legend.text = rownames(functional_states),
        args.legend = list(bty = "n", inset = -0.35, x = "right", 
                           ncol = 1, cex = 0.75, border = NA))
title(expression(paste("Functional ", bold("[E1]"))), cex = 1.5, line = 0.75)

barplot(as.matrix(habitat_states), col = names(habitat_order), 
        border = NA, cex.axis = 1.25, cex.names = 1.0, cex.lab = 1.25,
        ylab = "Proportion of Diet (%)", mgp=c(2.25,0.5,0), 
        ylim = c(0,100),
        axes = TRUE, las = 1,
        legend.text = rownames(habitat_states),
        args.legend = list(bty = "n", inset = -0.25, x = "right", 
                           ncol = 1, cex = 0.8, border = NA))
title(expression(paste("Habitat ", bold("[E2]"))), cex = 1.5, line = 0.75)

barplot(as.matrix(mbt_states), col = names(mbt_order), border = NA, 
        cex.axis = 1.25, cex.names = 1.0, cex.lab = 1.25,
        ylab = "Proportion of Diet (%)", mgp=c(2.25,0.5,0), 
        ylim = c(0,100),
        axes = TRUE, las = 1,
        legend.text = rownames(mbt_states),
        args.legend = list(bty = "n", inset = -0.2, x = "right", 
                           ncol = 1, cex = 0.9, border = NA))
title(expression(paste("MBT Prey Type ", bold("[E3]"))), cex = 1.5,  line = 0.75)
dev.off()

## End script
