library(viridisLite)
library(ggplot2)

########################### Plot Data with Low mu H s ##########################
FC_LmuHs.res <- readRDS("FusionChange_LowMuHighS_CYS.rds")
num_sims <- 25
col.vec <- viridis(num_sims)
num_gens <- nrow(FC_LmuHs.res)/(4*num_sims)
num_indv <- ncol(FC_LmuHs.res)
FC_LmuHs.count <- matrix(nrow = num_sims, ncol = num_gens)
FC_LmuHs.prop <- matrix(nrow = num_sims, ncol = num_gens)

X.count <- 0
Y.count <- 0
rows <- 1:4

for(sim in 1:nrow(FC_LmuHs.count)){
  print(paste(c("Sim: ", sim), collapse = ""))
  for(gen in 1:ncol(FC_LmuHs.count)){
    # Get count of all Y chr with fusions
    FC_LmuHs.count[sim, gen] <- sum(FC_LmuHs.res[rows[4],] > 0)
    # Get count of all Y chromosomes
    Y.count <- sum(FC_LmuHs.res[rows[2],])
    # Get proportion of Y chr with fusions
    FC_LmuHs.prop[sim, gen] <- FC_LmuHs.count[sim, gen]/Y.count
    
    # Get rows for next set
    rows <- rows + 4
  }
}

write.csv(FC_LmuHs.prop, file = "PropFusions_FC_LmuHs.csv", row.names = F)

plot(FC_LmuHs.prop[1,], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "Chiasmatic Y Low mu High s",
     col = col.vec[1])

for(i in 2:nrow(FC_LmuHs.prop)){
  lines(y = FC_LmuHs.prop[i, ], x = 1:ncol(FC_LmuHs.prop), col = col.vec[i])
}


########################### Plot Data with H mu Low s ##########################
FC_HmuLs.res <- readRDS("FusionChange_HighMuLowS_CYS.rds")
num_sims <- 25
col.vec <- viridis(num_sims)
num_gens <- nrow(FC_HmuLs.res)/(4*num_sims)
num_indv <- ncol(FC_HmuLs.res)
FC_HmuLs.count <- matrix(nrow = num_sims, ncol = num_gens)
FC_HmuLs.prop <- matrix(nrow = num_sims, ncol = num_gens)

X.count <- 0
Y.count <- 0
rows <- 1:4

for(sim in 1:nrow(FC_HmuLs.count)){
  print(paste(c("Sim: ", sim), collapse = ""))
  for(gen in 1:ncol(FC_HmuLs.count)){
    # Get count of all Y chr with fusions
    FC_HmuLs.count[sim, gen] <- sum(FC_HmuLs.res[rows[4],] > 0)
    # Get count of all Y chromosomes
    Y.count <- sum(FC_HmuLs.res[rows[2],])
    # Get proportion of Y chr with fusions
    FC_HmuLs.prop[sim, gen] <- FC_HmuLs.count[sim, gen]/Y.count
    
    # Get rows for next set
    rows <- rows + 4
  }
}

write.csv(FC_HmuLs.prop, file = "PropFusions_FC_HmuLs.csv", row.names = F)

plot(FC_HmuLs.prop[1,], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "Chiasmatic Y High mu Low s",
     col = col.vec[1])

for(i in 2:nrow(FC_HmuLs.prop)){
  lines(y = FC_HmuLs.prop[i, ], x = 1:ncol(FC_HmuLs.prop), col = col.vec[i])
}

