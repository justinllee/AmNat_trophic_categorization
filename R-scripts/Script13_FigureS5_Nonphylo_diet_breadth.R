### Script objectives:
### Non-Phylogenetic Estimations of Diet Breadth

## Prepare directories (use same PATH as previous scripts)
setwd("/PATH/Nonphylo BHC")
setwd(wd)

## Prepare data
Elapids <- read.csv("/PATH/Datasets/All_Elapids_squamatabase.csv")
taxalist <- read.csv("/PATH/Datasets/Elapid_full_taxon_list.csv")
breadth <- read.csv('diet_breadth_nonphylo_summary.csv')
colnames(breadth) <- c("radiation", "[T1]", "[T2]", "[T3]", "[E1]", "[E2]", "[E3]")

##----- Figure S4 boxplots ------

## Create color scheme
brew1 <- c("#E78AC3", "#FFD92F", "#FC8D62", "#8DA0CB","#66C2A5", "#A6D854")

## Plot
plot.new()
pdf(file = "FigS5_Diet_Breadth_Nonphylo.pdf", width = 9.75, height = 6.5)
par(mfrow = c(2,3), 
    oma = c(2,2,2,2), 
    mar = c(3,3,3,1),
    cex.axis = 1, las = 1,
    cex.lab = 1.25, 
    mgp=c(2,0.5,0),
    xpd=TRUE)

boxplot(breadth[,"[T1]"] ~ radiation, data = breadth, frame = TRUE, outline=FALSE,
        xlab = "Radiation", xaxt = "n",
        ylab = "log(Diet Breadth)", ylim = c(0, 2.5),
        cex = 0.5,
        col = NULL,
        border = c(brew1))
axis(side = 1, at = c(1:6), cex.axis = 1, label = c("A", "B", "C", "D", "E", "F"))
stripchart(breadth[,"[T1]"] ~ radiation, vertical = TRUE, data = breadth, 
           method = "jitter", add = TRUE, 
           pch = 20, 
           cex = 1.0, 
           col = c(brew1))
title("[T1]", line = 0.5)

boxplot(breadth[,"[T2]"] ~ radiation, data = breadth, frame = TRUE, outline=FALSE,
        xlab = "Radiation", xaxt = "n",
        ylab = "log(Diet Breadth)", ylim = c(0.5, 4.0),
        cex = 0.5,
        col = NULL,
        border = c(brew1))
axis(side = 1, at = c(1:6), cex.axis = 1, label = c("A", "B", "C", "D", "E", "F"))
stripchart(breadth[,"[T2]"] ~ radiation, vertical = TRUE, data = breadth, 
           method = "jitter", add = TRUE, 
           pch = 20, 
           cex = 1.0, 
           col = c(brew1))
title("[T2]", line = 0.5)

boxplot(breadth[,"[T3]"] ~ radiation, data = breadth, frame = TRUE, outline=FALSE,
        xlab = "Radiation", xaxt = "n",
        ylab = "log(Diet Breadth)", ylim = c(1, 5.0),
        cex = 0.5,
        col = NULL,
        border = c(brew1))
axis(side = 1, at = c(1:6), cex.axis = 1, label = c("A", "B", "C", "D", "E", "F"))
stripchart(breadth[,"[T3]"] ~ radiation, vertical = TRUE, data = breadth, 
           method = "jitter", add = TRUE, 
           pch = 20, 
           cex = 1.0, 
           col = c(brew1))
title("[T3]", line = 0.5)

boxplot(breadth[,"[E1]"] ~ radiation, data = breadth, frame = TRUE, outline=FALSE,
        xlab = "Radiation", xaxt = "n",
        ylab = "log(Diet Breadth)", ylim = c(0, 2.0),
        cex = 0.5,
        col = NULL,
        border = c(brew1))
axis(side = 1, at = c(1:6), cex.axis = 1, label = c("A", "B", "C", "D", "E", "F"))
stripchart(breadth[,"[E1]"] ~ radiation, vertical = TRUE, data = breadth, 
           method = "jitter", add = TRUE, 
           pch = 20, 
           cex = 1.0, 
           col = c(brew1))
title("[E1]", line = 0.5)

boxplot(breadth[,"[E2]"] ~ radiation, data = breadth, frame = TRUE, outline=FALSE,
        xlab = "Radiation", xaxt = "n",
        ylab = "log(Diet Breadth)", ylim = c(0, 1.2),
        cex = 0.5,
        col = NULL,
        border = c(brew1))
axis(side = 1, at = c(1:6), cex.axis = 1, label = c("A", "B", "C", "D", "E", "F"))
stripchart(breadth[,"[E2]"] ~ radiation, vertical = TRUE, data = breadth, 
           method = "jitter", add = TRUE, 
           pch = 20, 
           cex = 1.0, 
           col = c(brew1))
title("[E2]", line = 0.5)

boxplot(breadth[,"[E3]"] ~ radiation, data = breadth, frame = TRUE, outline=FALSE,
        xlab = "Radiation", xaxt = "n",
        ylab = "log(Diet Breadth)", ylim = c(0, 1.5),
        cex = 0.5,
        col = NULL,
        border = c(brew1))
axis(side = 1, at = c(1:6), cex.axis = 1, label = c("A", "B", "C", "D", "E", "F"))
stripchart(breadth[,"[E3]"] ~ radiation, vertical = TRUE, data = breadth, 
           method = "jitter", add = TRUE, 
           pch = 20, 
           cex = 1.0, 
           col = c(brew1))
title("[E3]", line = 0.5)
dev.off()

## End script
