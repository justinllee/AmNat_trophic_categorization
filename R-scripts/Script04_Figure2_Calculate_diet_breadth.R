### Script objectives:
### Estimate dietary breadth scores based on phylogenetic run
### of Dirichlet Process Multinomial State (DPMS) model.

## Prepare and install packages
library(macroevolution)
library(ape)
library(RColorBrewer)
library(phytools)

## Prepare directories, different method used here so that loop works
path <- "/PATH/"
path <- setwd(path)

## Import phylogeny
phy <- macroevolution::read.newick("/PATH/Datasets/All_Elapids_tree_v2.tre")
tree = read.tree("/PATH/Datasets/All_Elapids_tree_v2.tre")

## Prepare category schemes and associated datasets
schemes <- c("prey_traditional", "prey_phylo1", "prey_phylo2", 
             "prey_functional1", "prey_functional2", "prey_type")
setwd(path)

##-----Part 1: Calculate Diet Breadth------
breadth <- read.csv("/PATH/Datasets/Elapid_taxon_list_intree.csv")

for(MI in 1:length(schemes)){
  setwd(path) #resets working directory after each categorization scheme is finished
  folder <- sprintf("BHC_%s", schemes[MI])
  setwd(folder)
  
  ## Calculate diet breath scores
  db <- data.matrix(read.csv(sprintf("diet-proportions-%s.csv", schemes[MI]), row.names=1))
  diet.breadth <- apply(db, 1, function(p) 1 / sum(p*p))
  diet.breadth <- as.data.frame(diet.breadth)
  taxa <- row.names(diet.breadth)
  taxa <- read.csv("/PATH/Datasets/Elapid_taxon_list_intree.csv")
  diet.breadth <- cbind(taxa, diet.breadth$diet.breadth)
  labels <- c("species", "radiation", "diet_breadth")
  colnames(diet.breadth) <- labels
  write.csv(diet.breadth, sprintf("diet-breadth-%s.csv", schemes[MI]))
  
  # Add scores to summary file
  breadth[,2+MI] <- diet.breadth[,"diet_breadth"]
}

## Log transform the diet breadth scores
breadth[3:8] <- log(breadth[3:8])

## Rename columns and export diet breadth scores
setwd(path)
scheme.lab <- c("[T1]", "[T2]", "[T3]", "[E1]", "[E2]", "[E3]")
colnames(breadth)[3:8] <- scheme.lab
write.csv(breadth, file = "diet_breadth_phylo_summary.csv", row.names = F)

## ANOVA of diet breadth scores
bd <- breadth
anovasummary <- vector("list", 6)
tukeyHSD <- vector("list", 6)
for(an in 1:length(scheme.lab)){
  sch <- scheme.lab[an]
  anovaresult <- aov(as.matrix(bd[sch]) ~ bd$geog_radiation)
  anovasummary[an] <- summary(anovaresult)
  tukeyHSD[an] <- TukeyHSD(anovaresult, conf.level = 0.95)
}
names(tukeyHSD) <- scheme.lab
capture.output(tukeyHSD, file="TukeyHSD_diet_breadth_phylo.txt")


##-----Part 2: Boxplots of diet niche scores ------
# breadth <- read.csv('diet_breadth_phylo_summary.csv')
# colnames(breadth) <- c("taxon", "geog_radiation", 
# "[T1]", "[T2]", "[T3]", "[E1]", "[E2]", "[E3]")

brew1 <- c("#E78AC3", "#FFD92F", "#FC8D62", "#8DA0CB","#66C2A5", "#A6D854")

## Plot
plot.new()
pdf(file = "Fig2_Diet_Breadth_Phylo.pdf", width = 9.75, height = 6.5)
par(mfrow = c(2,3), 
    oma = c(2,2,2,2), 
    mar = c(3,3,3,2),
    mgp=c(2,0.5,0),
    xpd=TRUE)

boxplot(breadth[,"[T1]"] ~ geog_radiation, data = breadth, frame = TRUE, outline=FALSE,
        xlab = "Radiation", xaxt = "n",
        ylab = "log(Diet Breadth)", 
        ylim = c(1.75, 2.6), las = 1,
        cex = 0.75, cex.axis = 1, cex.lab = 1.25,
        col = NULL,
        border = c(brew1))
axis(side = 1, at = c(1:6), cex.axis = 1, label = c("A", "B", "C", "D", "E", "F"))
stripchart(breadth[,"[T1]"] ~ geog_radiation, vertical = TRUE, data = breadth, 
           method = "jitter", add = TRUE, 
           pch = 20, 
           cex = 1.0, 
           col = c(brew1))
title("[T1]", line = 0.5)

boxplot(breadth[,"[T2]"] ~ geog_radiation, data = breadth, frame = TRUE, outline=FALSE,
        xlab = "Radiation", xaxt = "n",
        ylab = "log(Diet Breadth)", 
        ylim = c(1.5, 4.5), las = 1,
        cex = 0.75, cex.axis = 1, cex.lab = 1.25,
        col = NULL,
        border = c(brew1))
axis(side = 1, at = c(1:6), cex.axis = 1, label = c("A", "B", "C", "D", "E", "F"))
stripchart(breadth[,"[T2]"] ~ geog_radiation, vertical = TRUE, data = breadth, 
           method = "jitter", add = TRUE, 
           pch = 20, 
           cex = 1.0, 
           col = c(brew1))
title("[T2]", line = 0.5)

boxplot(breadth[,"[T3]"] ~ geog_radiation, data = breadth, frame = TRUE, outline=FALSE,
        xlab = "Radiation", xaxt = "n",
        ylab = "log(Diet Breadth)", 
        ylim = c(2.0, 5.0), las = 1,
        cex = 0.75, cex.axis = 1, cex.lab = 1.25,
        col = NULL,
        border = c(brew1))
axis(side = 1, at = c(1:6), cex.axis = 1, label = c("A", "B", "C", "D", "E", "F"))
stripchart(breadth[,"[T3]"] ~ geog_radiation, vertical = TRUE, data = breadth, 
           method = "jitter", add = TRUE, 
           pch = 20, 
           cex = 1.0, 
           col = c(brew1))
title("[T3]", line = 0.5)

boxplot(breadth[,"[E1]"] ~ geog_radiation, data = breadth, frame = TRUE, outline=FALSE,
        xlab = "Radiation", xaxt = "n",
        ylab = "log(Diet Breadth)", 
        ylim = c(0.8, 2.3), las = 1,
        cex = 0.75, cex.axis = 1, cex.lab = 1.25,
        col = NULL,
        border = c(brew1))
axis(side = 1, at = c(1:6), cex.axis = 1, label = c("A", "B", "C", "D", "E", "F"))
stripchart(breadth[,"[E1]"] ~ geog_radiation, vertical = TRUE, data = breadth, 
           method = "jitter", add = TRUE, 
           pch = 20, 
           cex = 1.0, 
           col = c(brew1))
title("[E1]", line = 0.5)

boxplot(breadth[,"[E2]"] ~ geog_radiation, data = breadth, frame = TRUE, outline=FALSE,
        xlab = "Radiation", xaxt = "n",
        ylab = "log(Diet Breadth)", 
        ylim = c(0.6, 1.8), las = 1,
        cex = 0.75, cex.axis = 1, cex.lab = 1.2,
        col = NULL,
        border = c(brew1))
axis(side = 1, at = c(1:6), cex.axis = 1, label = c("A", "B", "C", "D", "E", "F"))
stripchart(breadth[,"[E2]"] ~ geog_radiation, vertical = TRUE, data = breadth, 
           method = "jitter", add = TRUE, 
           pch = 20, 
           cex = 1.0, 
           col = c(brew1))
title("[E2]", line = 0.5)

boxplot(breadth[,"[E3]"] ~ geog_radiation, data = breadth, frame = TRUE, outline=FALSE,
        xlab = "Radiation", xaxt = "n",
        ylab = "log(Diet Breadth)", 
        ylim = c(0.5, 1.8), las = 1,
        cex = 0.75, cex.axis = 1, cex.lab = 1.25,
        col = NULL,
        border = c(brew1))
axis(side = 1, at = c(1:6), cex.axis = 1, label = c("A", "B", "C", "D", "E", "F"))
stripchart(breadth[,"[E3]"] ~ geog_radiation, vertical = TRUE, data = breadth, 
           method = "jitter", add = TRUE, 
           pch = 20, 
           cex = 1.0, 
           col = c(brew1))
title("[E3]", line = 0.5)
dev.off()

## End script
