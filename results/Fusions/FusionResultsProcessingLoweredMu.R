library(viridisLite)
library(ggplot2)

########################### Plot Data with Lowered mu ##########################
loweredmu.res <- readRDS("LowEREDMuHighS_CYS.rds")
num_sims <- 25
col.vec <- viridis(num_sims)
num_gens <- nrow(loweredmu.res)/(4*num_sims)
num_indv <- ncol(loweredmu.res)
loweredmu.count <- matrix(nrow = num_sims, ncol = num_gens)
loweredmu.prop <- matrix(nrow = num_sims, ncol = num_gens)

X.count <- 0
Y.count <- 0
rows <- 1:4

for(sim in 1:nrow(loweredmu.count)){
  print(paste(c("Sim: ", sim), collapse = ""))
  for(gen in 1:ncol(loweredmu.count)){
    # Get count of all Y chr with fusions
    loweredmu.count[sim, gen] <- sum(loweredmu.res[rows[4],] > 0)
    # Get count of all Y chromosomes
    Y.count <- sum(loweredmu.res[rows[2],])
    # Get proportion of Y chr with fusions
    loweredmu.prop[sim, gen] <- loweredmu.count[sim, gen]/Y.count
    
    # Get rows for next set
    rows <- rows + 4
  }
}

write.csv(loweredmu.prop, file = "PropFusions_LoweredMu.csv", row.names = F)

plot(loweredmu.prop[1,], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "Chiasmatic Y Low mu High s",
     col = col.vec[1])

for(i in 2:nrow(loweredmu.prop)){
  lines(y = loweredmu.prop[i, ], x = 1:ncol(loweredmu.prop), col = col.vec[i])
}
