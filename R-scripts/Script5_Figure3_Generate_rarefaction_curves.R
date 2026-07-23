### Script objectives:
### Generate rarefaction curves to compare diets of Australo-Papuan hydrophiines
### and New World coralsnakes

## Prepare directories
wd <- "/PATH/"
setwd(wd)

## Prepare data
data <- "/PATH/Datasets/All_Elapids_squamatabase_resubmission1.csv"
data <- read.csv(data) # reads data in, but no radiations matched yet
good <- c("predator", "prey", "prey_traditional", "prey_phylo1", "prey_phylo2",
          "prey_functional1", "prey_functional2", "prey_type", "radiation")
cat <- c("prey", "traditional", "phylo1", "phylo2", "functional1", "functional2", "type")
data <- data[good]

#subset based on Coral Snakes and Hydrophiines
micrurus <- data[data$radiation == "nw-corals",]
australo <- data[data$radiation == "australo_papuan",]

##-----Part 1: Generate rarefaction curves for coral snakes ------
usp <- unique(micrurus$predator) #creates vector of unique species
Nsims <- 1000
nsamples <- seq(from = 1, to = 50, by = 1) #goes from 1 to 50
tmp <- matrix(data = NA, nrow = Nsims, ncol = length(cat), dimnames = list(1:Nsims, cat) )
m <- matrix(data = NA, nrow = length(nsamples), ncol = length(cat), dimnames = list(nsamples, cat) )
msd <- matrix(data = NA, nrow = length(nsamples), ncol = length(cat), dimnames = list(nsamples, cat) )

## for() loop that performs rarefaction analysis
for(h in 1:length(nsamples)){
  K <- nsamples[h] #reads in number of samples
  for(i in 1:Nsims){
  sub_sp <- sample(usp, K, replace = FALSE) #unique species, no. samples; replace = false
  for (ii in 1:length(sub_sp)){
  curr_sp <- sub_sp[ ii ] #chooses unique species sampled
  sub_df <- micrurus[ micrurus$predator == curr_sp , ] #creates data frame of diet records
  if (ii == 1){
    FINAL <- sub_df[ sample(1:nrow(sub_df), 1, replace = FALSE),  2:8] #stores prey category info
  } else {
    FINAL <- rbind(FINAL, sub_df[ sample(1:nrow(sub_df), 1),  2:8])
  }
  }
  for (iii in 1:length(cat)){
    si <- prop.table(table(FINAL[,iii])) #creates diet proportion table
    invsim <- 1 / sum(si * si) #inverse Simpson's diversity index
    tmp[i, iii] <- invsim
  }
  rm(FINAL) #removes final so you start from scratch next simulation
  }
  for(iv in 1:length(cat)){
    avg <- mean(tmp[,iv])
    stdev <- sd(tmp[,iv])
    m[h,iv] <- avg
    msd[h,iv] <- stdev
  }
}
curves <- as.data.frame(m)
sd_curves <- as.data.frame(msd)
write.csv(curves, file = "Micrurus_mean_simpson_results.csv")
write.csv(sd_curves, file = "Micrurus_sd_simpson_results.csv")

##-----Part 2: Same loop applied to Australian Hydrophiines ------
asp <- unique(australo$predator) #creates vector of unique species
atmp <- matrix(data = NA, nrow = Nsims, ncol = length(cat), dimnames = list(1:Nsims, cat) )
ap <- matrix(data = NA, nrow = length(nsamples), ncol = length(cat), dimnames = list(nsamples, cat) )
asd <- matrix(data = NA, nrow = length(nsamples), ncol = length(cat), dimnames = list(nsamples, cat) )

#Loop
for(h in 1:length(nsamples)){
  K <- nsamples[h] #reads in number of samples
  for(i in 1:Nsims){
    sub_sp <- sample(asp, K, replace = FALSE) #unique species, no. samples; replace = false
    for (ii in 1:length(sub_sp)){
      curr_sp <- sub_sp[ ii ] #chooses unique species sampled
      sub_df <- australo[ australo$predator == curr_sp , ] #creates data frame of diet records
      if (ii == 1){
        FINAL <- sub_df[ sample(1:nrow(sub_df), 1, replace = FALSE),  2:8] #stores prey category info
      } else {
        FINAL <- rbind(FINAL, sub_df[ sample(1:nrow(sub_df), 1),  2:8])
      }
    }
    for (iii in 1:length(cat)){
      si <- prop.table(table(FINAL[,iii])) #creates diet proportion table
      invsim <- 1 / sum(si * si) #inverse Simpson's diversity index
      atmp[i, iii] <- invsim
    }
    rm(FINAL) #removes final so you start from scratch next simulation
  }
  for(iv in 1:length(cat)){
    aavg <- mean(atmp[,iv])
    astdev <- sd(atmp[,iv])
    ap[h,iv] <- aavg
    asd[h,iv] <- astdev
  }
}
acurves <- as.data.frame(ap)
asd_curves <- as.data.frame(asd)
write.csv(acurves, file = "Aus_mean_simpson_results.csv")
write.csv(asd_curves, file = "Aus_sd_simpson_results.csv")

##-----Part 3: Plot rarefaction curves from each category ------

## Read in the data
curves <- read.csv("Micrurus_mean_simpson_results.csv")
sd_curves <- read.csv("Micrurus_sd_simpson_results.csv")
acurves <- read.csv("Aus_mean_simpson_results.csv")
asd_curves <- read.csv("Aus_sd_simpson_results.csv")

## Generate color scheme
cor <- "#66C2A5"
aus <- "#FC8D62"

## Plot
plot.new()
pdf(file = "Fig3_Rarefaction_Curves.pdf", width = 9.75, height = 6.5)
par(mfrow = c(2,3), 
    oma = c(2,2,2,2), 
    mar = c(3,3,3,2),
    cex.axis = 0.75,
    cex.lab = 0.75,
    mgp=c(1.75,0.5,0),
    xpd=TRUE)
ticks = c(0, 10, 20, 30, 40, 50)
breaks = c(0, 10, 20, 30, 40, 50)

plot(curves$traditional, pch = 21, bg = cor, col = F,
     xlab = 'Samples', xaxt = "n", cex.lab = 1.25, cex.axis = 1.25,
     ylab = 'Mean No. Sampled Prey', ylim = c(0,4), las = 1)
arrows(c(1:50), curves$traditional-sd_curves$traditional, 
       c(1:50), curves$traditional+sd_curves$traditional, 
       length=0.0, angle=90, code=1, col = cor)
points(acurves$traditional, pch = 21, bg = aus, col = F)
arrows(c(1:50), acurves$traditional-asd_curves$traditional, 
       c(1:50), acurves$traditional+asd_curves$traditional, 
       length=0.0, angle=90, code=1, col = aus)
title(main = '[T1]', line = 0.75)
axis(1, at=breaks, labels=ticks, cex.axis = 1.25)
legend("bottomright", legend=c("Coral Snakes", "Hydrophiines"), col=c(cor, aus),
       pch=20, cex= 1, ncol = 1, bty = "n")

plot(curves$phylo1, pch = 21, bg = cor, col = F,
     xlab = 'Samples', xaxt = "n", cex.lab = 1.25, cex.axis = 1.25,
     ylab = 'Mean No. Sampled Prey', ylim = c(0,5.5), las = 1)
arrows(c(1:50), curves$phylo1-sd_curves$phylo1, 
       c(1:50), curves$phylo1+sd_curves$phylo1, 
       length=0.0, angle=90, code=1, col = cor)
points(acurves$phylo1, pch = 21, bg = aus, col = F)
arrows(c(1:50), acurves$phylo1-asd_curves$phylo1, 
       c(1:50), acurves$phylo1+asd_curves$phylo1, 
       length=0.0, angle=90, code=1, col = aus)
title(main = '[T2]', line = 0.75)
axis(1, at=breaks, labels=ticks, cex.axis = 1.25)
legend("bottomright", legend=c("Coral Snakes", "Hydrophiines"), col=c(cor, aus),
       pch=20, cex= 1, ncol = 1, bty = "n")

plot(curves$phylo2, pch = 21, bg = cor, col = F,
     xlab = 'Samples', xaxt = "n", cex.lab = 1.25, cex.axis = 1.25,
     ylab = 'Mean No. Sampled Prey', ylim = c(0,6.5), las = 1)
arrows(c(1:50), curves$phylo2-sd_curves$phylo2, 
       c(1:50), curves$phylo2+sd_curves$phylo2, 
       length=0.0, angle=90, code=1, col = cor)
points(acurves$phylo2, pch = 21, bg = aus, col = F)
arrows(c(1:50), acurves$phylo2-asd_curves$phylo2, 
       c(1:50), acurves$phylo2+asd_curves$phylo2, 
       length=0.0, angle=90, code=1, col = aus)
title(main = '[T3]', line = 0.75)
axis(1, at=breaks, labels=ticks, cex.axis = 1.25)
legend("bottomright", legend=c("Coral Snakes", "Hydrophiines"), col=c(cor, aus),
       pch=20, cex= 1, ncol = 1, bty = "n")

plot(curves$functional1, pch = 21, bg = cor, col = F,
     xlab = 'Samples', xaxt = "n", cex.lab = 1.25, cex.axis = 1.25,
     ylab = 'Mean No. Sampled Prey', ylim = c(0,4), las = 1)
arrows(c(1:50), curves$functional1-sd_curves$functional1, 
       c(1:50), curves$functional1+sd_curves$functional1, 
       length=0.0, angle=90, code=1, col = cor)
points(acurves$functional1, pch = 21, bg = aus, col = F)
arrows(c(1:50), acurves$functional1-asd_curves$functional1, 
       c(1:50), acurves$functional1+asd_curves$functional1, 
       length=0.0, angle=90, code=1, col = aus)
title(main = '[E1]', line = 0.75)
axis(1, at=breaks, labels=ticks, cex.axis = 1.25)
legend("bottomright", legend=c("Coral Snakes", "Hydrophiines"), col=c(cor, aus),
       pch=20, cex= 1, ncol = 1, bty = "n")

plot(curves$functional2, pch = 21, bg = cor, col = F,
     xlab = 'Samples', xaxt = "n", cex.lab = 1.25, cex.axis = 1.25,
     ylab = 'Mean No. Sampled Prey', ylim = c(0,4), las = 1)
arrows(c(1:50), curves$functional2-sd_curves$functional2, 
       c(1:50), curves$functional2+sd_curves$functional2, 
       length=0.0, angle=90, code=1, col = cor)
points(acurves$functional2, pch = 21, bg = aus, col = F)
arrows(c(1:50), acurves$functional2-asd_curves$functional2, 
       c(1:50), acurves$functional2+asd_curves$functional2, 
       length=0.0, angle=90, code=1, col = aus)
title(main = '[E2]', line = 0.75)
axis(1, at=breaks, labels=ticks, cex.axis = 1.25)
legend("bottomright", legend=c("Coral Snakes", "Hydrophiines"), col=c(cor, aus),
       pch=20, cex= 1, ncol = 1, bty = "n")

plot(curves$type, pch = 21, bg = cor, col = F,
     xlab = 'Samples', xaxt = "n", cex.lab = 1.25, cex.axis = 1.25,
     ylab = 'Mean No. Sampled Prey', ylim = c(0,4), las = 1)
arrows(c(1:50), curves$type-sd_curves$type, 
       c(1:50), curves$type+sd_curves$type, 
       length=0.0, angle=90, code=1, col = cor)
points(acurves$type, pch = 21, bg = aus, col = F)
arrows(c(1:50), acurves$type-asd_curves$type, 
       c(1:50), acurves$type+asd_curves$type, 
       length=0.0, angle=90, code=1, col = aus)
title(main = '[E3]', line = 0.75)
axis(1, at=breaks, labels=ticks, cex.axis = 1.25)
legend("bottomright", legend=c("Coral Snakes", "Hydrophiines"), col=c(cor, aus),
       pch=20, cex=1, ncol = 1, bty = "n")
dev.off()

## End script
