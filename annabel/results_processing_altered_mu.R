# The goal of this script is to etermine whether changing the way mutations are
# introduced fixes the problem where none of the fusions were fixing 

library(viridisLite)
library(ggplot2)

setwd("/Users/knigh/Documents/GitHub/achiasmy-sims/annabel/")

# Read in data from running simulations where...
# 1. Males are chiasmatic
# 2. Fusions can only occur between Y chromosomes and small (S) autosomes
CYS_HmuLs.res <- readRDS("HmuLs_CYS_06-10-22.rds")
CYS_LmuHs.res <- readRDS("LmuHs_CYS_06-10-22.rds")

############################## General Parameters ##############################
# Number of simulations must be known a priori
num_sims <- 25
col.vec <- viridis(num_sims)
# Number of generations is number of rows in a matrix divided by 4*num_sims
num_gens <- nrow(CYS_LmuHs.res)/(4*num_sims)
# Number of individuals is number of columns in a matrix
num_indv <- ncol(CYS_LmuHs.res)

# Create matrices to store the counts and proportions of fusions of each 
# simulation (row) at each generation (column)
CYS_HmuLs.count <- CYS_HmuLs.prop <- CYS_LmuHs.count <- CYS_LmuHs.prop <- 
  matrix(nrow = num_sims, ncol = num_gens)

Y.count <- 0
rows <- 1:4

for(sim in 1:num_sims){
  print(paste(c("Sim: ", sim), collapse = ""))
  for(gen in 1:num_gens){
    # For High mu Low s simulation...
    # Get count of all Y chr with fusions
    CYS_HmuLs.count[sim, gen] <- sum(CYS_HmuLs.res[rows[4],] > 0)
    # Get count of all Y chromosomes
    Y.count <- sum(CYS_HmuLs.res[rows[2],])
    # Get proportion of Y chr with fusions
    CYS_HmuLs.prop[sim, gen] <- CYS_HmuLs.count[sim, gen]/Y.count
    
    # Repeat for Low mu, High s simulation
    CYS_LmuHs.count[sim, gen] <- sum(CYS_LmuHs.res[rows[4],] > 0)
    Y.count <- sum(CYS_LmuHs.res[rows[2],])
    CYS_LmuHs.prop[sim, gen] <- CYS_LmuHs.count[sim, gen]/Y.count
    
    # Get rows for next set
    rows <- rows + 4
  }
}

# Plot fixation patterns Over generations for each
plot(CYS_HmuLs.prop[1,], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "Chiasmatic Y High mu Low s",
     col = col.vec[1])

for(i in 2:num_sims){
  lines(y = CYS_HmuLs.prop[i, ], x = 1:num_gens, col = col.vec[i])
}


plot(CYS_LmuHs.prop[1,], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "Chiasmatic Y Low mu High s",
     col = col.vec[1])

for(i in 2:num_sims){
  lines(y = CYS_LmuHs.prop[i, ], x = 1:num_gens, col = col.vec[i])
}
