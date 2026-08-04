### Script objectives:
### Prepare Elapid diet datasets for Dirichlet Process Multinomial State (DPMS)
### model.

## Prepare packages
library(ape)

## Prepare directories
path <- "/PATH/"
path <- setwd(path)

## Prepare files
setwd("/PATH/")
Elapids <- read.csv("/PATH/Datasets/All_Elapids_squamatabase.csv")
tree <- read.tree("/PATH/Datasets/All_Elapids_tree_v2.tre") 

## Adjust prey/predator counts for missing data
Elapids[Elapids == ""] <- NA
Elapids$predator_count[is.na(Elapids$predator_count)] <- 1
Elapids$prey_count[is.na(Elapids$prey_count)] <- 1 

## Change taxon names to underscores
Elapids$predator <- chartr(" ", "_", Elapids$predator)
Elapids$predator_taxon <- chartr(" ", "_", Elapids$predator_taxon)
Elapids$prey <- chartr(" ", "_", Elapids$prey)
Elapids$prey_taxon <- chartr(" ", "_", Elapids$prey_taxon)

## Organize prey category schemes
schemes <- c("prey_traditional", "prey_phylo1", "prey_phylo2", 
             "prey_functional1", "prey_functional2", "prey_type")

## Count prey categories in each predator species
for(j in 1:length(schemes)){
  sch <- schemes[j]
  Elapids <- Elapids[!is.na(Elapids[sch]),]
}

for(j in 1:length(schemes)){
  sch <- schemes[j]
  categories <- unique(Elapids[, sch])
  categories <- categories[!is.na(categories)]
  species <- unique(Elapids$predator)
  diet_matrix <- matrix(data = 0, nrow = length(species), ncol = length(categories))
  rownames(diet_matrix) <- species
  colnames(diet_matrix) <- categories
  
  for(a in 1:length(species)){
    taxon <- Elapids[Elapids$predator == species[a], ]
    for( jj in 1:length(taxon$prey)){
      ii <- taxon$predator[jj]
      kk <- taxon[jj,sch]
      uv <- taxon$predator_voucher[jj]
      vv <- match(taxon$predator_voucher, uv)
      ww <- which(vv == 1)
      tmp <- taxon[ww, ]
      maxtot <- max(taxon$predator_count[jj], taxon$prey_count[jj])
      if( any(is.na(tmp$predator_voucher)) == TRUE ){
        diet_matrix[ii, kk] <- diet_matrix[ii, kk] + maxtot
      } else if( any(is.na(tmp$predator_voucher)) == FALSE ){
        count_row <- nrow(tmp)
        if(count_row == 1){
          diet_matrix[ii, kk] <- diet_matrix[ii, kk] + maxtot
        } else if(jj == ww[1]){
          diet_matrix[ii, kk] <- diet_matrix[ii, kk] + maxtot
        } else if(jj != ww[1]){
          prey <- taxon$prey[jj]
          if(prey == tmp$prey[1]){
            diet_matrix[ii, kk] <- diet_matrix[ii, kk] + 0
          } else {
            diet_matrix[ii, kk] <- diet_matrix[ii, kk] + maxtot
          }
        }
      }
    }
  }
  
  #append prey category names to matrix
  diet_matrix <- as.data.frame(diet_matrix)
  datafile <- sprintf("diet_matrix_%s.csv", sch)
  write.csv(diet_matrix, file = datafile, row.names = TRUE)
  rm(diet_matrix)
  
  #Prepare dataset for BHC clusters 
  bhc_dataset <- data.frame("species"= numeric(0), "food"= numeric(0), "count" = numeric(0))
  
  for(a in 1:length(species)){
    taxon <- Elapids[Elapids$predator == species[a], ]
    for(b in 1:length(categories)){
      count <- sum(taxon[, sch] == categories[b])
      if(count != 0){
        data <- c(species[a], categories[b], count)
        bhc_dataset[nrow(bhc_dataset)+1, ] <- data
        bhc_dataset
      }
    }
  }
  
  bhcfile <- sprintf("bhc_dataset_%s.csv", sch)
  write.csv(bhc_dataset, file = bhcfile, row.names = FALSE)
  rm(bhc_dataset)  
}

## End script
