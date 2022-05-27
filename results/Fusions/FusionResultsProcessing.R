library(viridisLite)
library(ggplot2)

###################### Identify Appropriate mu and s ###########################
setwd("/Users/knigh/Documents/GitHub/achiasmy-sims/results/Fusions/")

# Data from running simulation on...
# A-_----:  Heterogametic sex is achiasmatic
# C-_----:  Heterogametic sex is chiasmatic
# -X_----:  Fusions occur on X
# -Y_----:  Fusions occur on Y
# --_Lmu--: Low mutation rate
# --_Hmu--: High mutation rate
# --_--Ls:  Low selection coefficient on SAL
# --_--Hs:  High selection coefficient on SAL
CX_LmuLs.res <- readRDS("XChiasmaticResults_LowMuLowS.rds")
CX_HmuLs.res <- readRDS("XChiasmaticResults_HighMuLowS.rds")
CX_HmuHs.res <- readRDS("XChiasmaticResults_HighMuHighS.rds")
CX_LmuHs.res <- readRDS("XChiasmaticResults_LowMuHighS.rds")

CY_LmuLs.res <- readRDS("YChiasmaticResults_LowMuLowS.rds")
CY_HmuLs.res <- readRDS("YChiasmaticResults_HighMuLowS.rds")
CY_HmuHs.res <- readRDS("YChiasmaticResults_HighMuHighS.rds")
CY_LmuHs.res <- readRDS("YChiasmaticResults_LowMuHighS.rds")

AX_LmuLs.res <- readRDS("XAchiasmaticResults_LowMuLowS.rds")
AX_HmuLs.res <- readRDS("XAchiasmaticResults_HighMuLowS.rds")
AX_HmuHs.res <- readRDS("XAchiasmaticResults_HighMuHighS.rds")
AX_LmuHs.res <- readRDS("XAchiasmaticResults_LowMuHighS.rds")

AY_LmuLs.res <- readRDS("YAchiasmaticResults_LowMuLowS.rds")
AY_HmuLs.res <- readRDS("YAchiasmaticResults_HighMuLowS.rds")
AY_HmuHs.res <- readRDS("YAchiasmaticResults_HighMuHighS.rds")
AY_LmuHs.res <- readRDS("YAchiasmaticResults_LowMuHighS.rds")
############################## General Parameters ##############################
# Number of simulations must be known a priori
num_sims <- 25
col.vec <- viridis(num_sims)
# Number of generations is number of rows in a matrix divided by 4*num_sims
num_gens <- nrow(CX_LmuLs.res)/(4*num_sims)
# Number of individuals is number of columns in a matrix
num_indv <- ncol(CX_LmuLs.res)

# Get total number of fusions at each generation for each simulation
# Matrices:
# Columns = Generations
# Rows = Simulations
# Cell = Total count OR proportion of fusions (all types)
CX_LmuLs.count <- matrix(nrow = num_sims, ncol = num_gens)
CX_LmuLs.prop <- matrix(nrow = num_sims, ncol = num_gens)

CX_HmuLs.count <- matrix(nrow = num_sims, ncol = num_gens)
CX_HmuLs.prop <- matrix(nrow = num_sims, ncol = num_gens)

CX_HmuHs.count <- matrix(nrow = num_sims, ncol = num_gens)
CX_HmuHs.prop <- matrix(nrow = num_sims, ncol = num_gens)

CX_LmuHs.count <- matrix(nrow = num_sims, ncol = num_gens)
CX_LmuHs.prop <- matrix(nrow = num_sims, ncol = num_gens)

AX_LmuLs.count <- matrix(nrow = num_sims, ncol = num_gens)
AX_LmuLs.prop <- matrix(nrow = num_sims, ncol = num_gens)

AX_HmuLs.count <- matrix(nrow = num_sims, ncol = num_gens)
AX_HmuLs.prop <- matrix(nrow = num_sims, ncol = num_gens)

AX_HmuHs.count <- matrix(nrow = num_sims, ncol = num_gens)
AX_HmuHs.prop <- matrix(nrow = num_sims, ncol = num_gens)

AX_LmuHs.count <- matrix(nrow = num_sims, ncol = num_gens)
AX_LmuHs.prop <- matrix(nrow = num_sims, ncol = num_gens)

CY_LmuLs.count <- matrix(nrow = num_sims, ncol = num_gens)
CY_LmuLs.prop <- matrix(nrow = num_sims, ncol = num_gens)

CY_HmuLs.count <- matrix(nrow = num_sims, ncol = num_gens)
CY_HmuLs.prop <- matrix(nrow = num_sims, ncol = num_gens)

CY_HmuHs.count <- matrix(nrow = num_sims, ncol = num_gens)
CY_HmuHs.prop <- matrix(nrow = num_sims, ncol = num_gens)

CY_LmuHs.count <- matrix(nrow = num_sims, ncol = num_gens)
CY_LmuHs.prop <- matrix(nrow = num_sims, ncol = num_gens)

AY_LmuLs.count <- matrix(nrow = num_sims, ncol = num_gens)
AY_LmuLs.prop <- matrix(nrow = num_sims, ncol = num_gens)

AY_HmuLs.count <- matrix(nrow = num_sims, ncol = num_gens)
AY_HmuLs.prop <- matrix(nrow = num_sims, ncol = num_gens)

AY_HmuHs.count <- matrix(nrow = num_sims, ncol = num_gens)
AY_HmuHs.prop <- matrix(nrow = num_sims, ncol = num_gens)

AY_LmuHs.count <- matrix(nrow = num_sims, ncol = num_gens)
AY_LmuHs.prop <- matrix(nrow = num_sims, ncol = num_gens)

X.count <- 0
Y.count <- 0
rows <- 1:4

for(sim in 1:nrow(CX_LmuLs.count)){
  print(paste(c("Sim: ", sim), collapse = ""))
  for(gen in 1:ncol(CX_LmuLs.count)){
    # X
    # Add count of fusions at current gen + sim to appropriate cell
    CX_LmuLs.count[sim, gen] <- sum(CX_LmuLs.res[rows[3],] > 0) + sum(CX_LmuLs.res[rows[4],] > 0)
    # Get number of X chromosomes at current gen + sim
    X.count <- num_indv + sum(CX_LmuLs.res[rows[2],] == 0)
    # Get proportion of X chromosomes with fusions and add to appropriate cell
    CX_LmuLs.prop[sim, gen] <- CX_LmuLs.count[sim, gen]/X.count
    
    # Repeat for all other X chromosome simulations
    CX_HmuLs.count[sim, gen] <- sum(CX_HmuLs.res[rows[3],] > 0) + sum(CX_HmuLs.res[rows[4],] > 0)
    X.count <- num_indv + sum(CX_HmuLs.res[rows[2],] == 0)
    CX_HmuLs.prop[sim, gen] <- CX_HmuLs.count[sim, gen]/X.count
    
    CX_HmuHs.count[sim, gen] <- sum(CX_HmuHs.res[rows[3],] > 0) + sum(CX_HmuHs.res[rows[4],] > 0)
    X.count <- num_indv + sum(CX_HmuHs.res[rows[2],] == 0)
    CX_HmuHs.prop[sim, gen] <- CX_HmuHs.count[sim, gen]/X.count
    
    CX_LmuHs.count[sim, gen] <- sum(CX_LmuHs.res[rows[3],] > 0) + sum(CX_LmuHs.res[rows[4],] > 0)
    X.count <- num_indv + sum(CX_LmuHs.res[rows[2],] == 0)
    CX_LmuHs.prop[sim, gen] <- CX_LmuHs.count[sim, gen]/X.count
    
    AX_LmuLs.count[sim, gen] <- sum(AX_LmuLs.res[rows[3],] > 0) + sum(AX_LmuLs.res[rows[4],] > 0)
    X.count <- num_indv + sum(AX_LmuLs.res[rows[2],] == 0)
    AX_LmuLs.prop[sim, gen] <- AX_LmuLs.count[sim, gen]/X.count
    
    AX_HmuLs.count[sim, gen] <- sum(AX_HmuLs.res[rows[3],] > 0) + sum(AX_HmuLs.res[rows[4],] > 0)
    X.count <- num_indv + sum(AX_HmuLs.res[rows[2],] == 0)
    AX_HmuLs.prop[sim, gen] <- AX_HmuLs.count[sim, gen]/X.count
    
    AX_HmuHs.count[sim, gen] <- sum(AX_HmuHs.res[rows[3],] > 0) + sum(AX_HmuHs.res[rows[4],] > 0)
    X.count <- num_indv + sum(AX_HmuHs.res[rows[2],] == 0)
    AX_HmuHs.prop[sim, gen] <- AX_HmuHs.count[sim, gen]/X.count
    
    AX_LmuHs.count[sim, gen] <- sum(AX_LmuHs.res[rows[3],] > 0) + sum(AX_LmuHs.res[rows[4],] > 0)
    X.count <- num_indv + sum(AX_LmuHs.res[rows[2],] == 0)
    AX_LmuHs.prop[sim, gen] <- AX_LmuHs.count[sim, gen]/X.count
    
    # Y
    # Get count of all Y chr with fusions
    CY_LmuLs.count[sim, gen] <- sum(CY_LmuLs.res[rows[4],] > 0)
    # Get count of all Y chromosomes
    Y.count <- sum(CY_LmuLs.res[rows[2],])
    # Get proportion of Y chr with fusions
    CY_LmuLs.prop[sim, gen] <- CY_LmuLs.count[sim, gen]/Y.count
    
    CY_HmuLs.count[sim, gen] <- sum(CY_HmuLs.res[rows[4],] > 0)
    Y.count <- sum(CY_HmuLs.res[rows[2],])
    CY_HmuLs.prop[sim, gen] <- CY_HmuLs.count[sim, gen]/Y.count
    
    CY_HmuHs.count[sim, gen] <- sum(CY_HmuHs.res[rows[4],] > 0)
    Y.count <- sum(CY_HmuHs.res[rows[2],])
    CY_HmuHs.prop[sim, gen] <- CY_HmuHs.count[sim, gen]/Y.count
    
    CY_LmuHs.count[sim, gen] <- sum(CY_LmuHs.res[rows[4],] > 0)
    Y.count <- sum(CY_LmuHs.res[rows[2],])
    CY_LmuHs.prop[sim, gen] <- CY_LmuHs.count[sim, gen]/Y.count
    
    AY_LmuLs.count[sim, gen] <- sum(AY_LmuLs.res[rows[4],] > 0)
    Y.count <- sum(AY_LmuLs.res[rows[2],])
    AY_LmuLs.prop[sim, gen] <- AY_LmuLs.count[sim, gen]/Y.count
    
    AY_HmuLs.count[sim, gen] <- sum(AY_HmuLs.res[rows[4],] > 0)
    Y.count <- sum(AY_HmuLs.res[rows[2],])
    AY_HmuLs.prop[sim, gen] <- AY_HmuLs.count[sim, gen]/Y.count
    
    AY_HmuHs.count[sim, gen] <- sum(AY_HmuHs.res[rows[4],] > 0)
    Y.count <- sum(AY_HmuHs.res[rows[2],])
    AY_HmuHs.prop[sim, gen] <- AY_HmuHs.count[sim, gen]/Y.count
    
    AY_LmuHs.count[sim, gen] <- sum(AY_LmuHs.res[rows[4],] > 0)
    Y.count <- sum(AY_LmuHs.res[rows[2],])
    AY_LmuHs.prop[sim, gen] <- AY_LmuHs.count[sim, gen]/Y.count
    
    # Get rows for next set
    rows <- rows + 4
  }
}

# Write proportion of all fusions under each condition to a file
write.csv(CX_LmuLs.prop, file = "PropFusionsCXLmuLs.csv", row.names = F)
write.csv(CX_HmuLs.prop, file = "PropFusionsCXHmuLs.csv", row.names = F)
write.csv(CX_LmuHs.prop, file = "PropFusionsCXLmuHs.csv", row.names = F)
write.csv(CX_HmuHs.prop, file = "PropFusionsCXHmuHs.csv", row.names = F)

write.csv(CY_LmuLs.prop, file = "PropFusionsCYLmuLs.csv", row.names = F)
write.csv(CY_HmuLs.prop, file = "PropFusionsCYHmuLs.csv", row.names = F)
write.csv(CY_LmuHs.prop, file = "PropFusionsCYLmuHs.csv", row.names = F)
write.csv(CY_HmuHs.prop, file = "PropFusionsCYHmuHs.csv", row.names = F)

write.csv(AX_LmuLs.prop, file = "PropFusionsAXLmuLs.csv", row.names = F)
write.csv(AX_HmuLs.prop, file = "PropFusionsAXHmuLs.csv", row.names = F)
write.csv(AX_LmuHs.prop, file = "PropFusionsAXLmuHs.csv", row.names = F)
write.csv(AX_HmuHs.prop, file = "PropFusionsAXHmuHs.csv", row.names = F)

write.csv(AY_LmuLs.prop, file = "PropFusionsAYLmuLs.csv", row.names = F)
write.csv(AY_HmuLs.prop, file = "PropFusionsAYHmuLs.csv", row.names = F)
write.csv(AY_LmuHs.prop, file = "PropFusionsAYLmuHs.csv", row.names = F)
write.csv(AY_HmuHs.prop, file = "PropFusionsAYHmuHs.csv", row.names = F)

# Ensure sims with Y fusions never have fusions to X
# (This was a bug fixed between Mar15 and Mar 26. Proportion of Y fusions is #
# fusion to X/Y over total # Y chr, so if there are fusion to X total prop at
# fixation will be > 1)
sum(c(CY_LmuLs.prop > 1, CY_HmuLs.prop > 1, CY_LmuHs.prop > 1, CY_HmuHs.prop > 1,
      AY_LmuLs.prop > 1, AY_HmuLs.prop > 1, AY_LmuHs.prop > 1, AY_HmuHs.prop > 1))

# Plot fixation patterns Over generations for each
plot(CX_LmuLs.prop[1,], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "Chiasmatic X Low mu Low s",
     col = col.vec[1])

for(i in 2:nrow(CX_LmuLs.prop)){
  lines(y = CX_LmuLs.prop[i, ], x = 1:ncol(CX_LmuLs.prop), col = col.vec[i])
}

plot(CX_HmuHs.prop[1,], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions",
     main = "Chiasmatic X High mu High s",
     col = col.vec[1])

for(i in 2:nrow(CX_HmuHs.prop)){
  lines(y = CX_HmuHs.prop[i, ], x = 1:ncol(CX_HmuHs.prop), col = col.vec[i])
}

plot(CX_HmuLs.prop[1,], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions",
     main = "Chiasmatic X High mu Low s",
     col = col.vec[1])

for(i in 2:nrow(CX_HmuLs.prop)){
  lines(y = CX_HmuLs.prop[i, ], x = 1:ncol(CX_HmuLs.prop), col = col.vec[i])
}

plot(CX_LmuHs.prop[1,], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions",
     main = "Chiasmatic X Low mu High s",
     col = col.vec[1])

for(i in 2:nrow(CX_LmuHs.prop)){
  lines(y = CX_LmuHs.prop[i, ], x = 1:ncol(CX_LmuHs.prop), col = col.vec[i])
}

plot(AX_LmuLs.prop[1,], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions",
     main = "Achiasmatic X Low mu Low s",
     col = col.vec[1])

for(i in 2:nrow(AX_LmuLs.prop)){
  lines(y = AX_LmuLs.prop[i, ], x = 1:ncol(AX_LmuLs.prop), col = col.vec[i])
}

plot(AX_HmuHs.prop[1,], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions",
     main = "Achiasmatic X High mu High s",
     col = col.vec[1])

for(i in 2:nrow(AX_HmuHs.prop)){
  lines(y = AX_HmuHs.prop[i, ], x = 1:ncol(AX_HmuHs.prop), col = col.vec[i])
}

plot(AX_HmuLs.prop[1,], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions",
     main = "Achiasmatic X High mu Low s",
     col = col.vec[1])

for(i in 2:nrow(AX_HmuLs.prop)){
  lines(y = AX_HmuLs.prop[i, ], x = 1:ncol(AX_HmuLs.prop), col = col.vec[i])
}

plot(AX_LmuHs.prop[1,], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions",
     main = "Achiasmatic X Low mu High s",
     col = col.vec[1])

for(i in 2:nrow(AX_LmuHs.prop)){
  lines(y = AX_LmuHs.prop[i, ], x = 1:ncol(AX_LmuHs.prop), col = col.vec[i])
}


plot(CY_LmuLs.prop[1,], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions",
     main = "Chiasmatic Y Low mu Low s",
     col = col.vec[1])

for(i in 2:nrow(CY_LmuLs.prop)){
  lines(y = CY_LmuLs.prop[i, ], x = 1:ncol(CY_LmuLs.prop), col = col.vec[i])
}

plot(CY_HmuHs.prop[1,], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions",
     main = "Chiasmatic Y High mu High s",
     col = col.vec[1])

for(i in 2:nrow(CY_HmuHs.prop)){
  lines(y = CY_HmuHs.prop[i, ], x = 1:ncol(CY_HmuHs.prop), col = col.vec[i])
}

plot(CY_HmuLs.prop[1,], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions",
     main = "Chiasmatic Y High mu Low s",
     col = col.vec[1])

for(i in 2:nrow(CY_HmuLs.prop)){
  lines(y = CY_HmuLs.prop[i, ], x = 1:ncol(CY_HmuLs.prop), col = col.vec[i])
}

plot(CY_LmuHs.prop[1,], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions",
     main = "Chiasmatic Y Low mu High s",
     col = col.vec[1])

for(i in 2:nrow(CY_LmuHs.prop)){
  lines(y = CY_LmuHs.prop[i, ], x = 1:ncol(CY_LmuHs.prop), col = col.vec[i])
}

plot(AY_LmuLs.prop[1,], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions",
     main = "Achiasmatic Y Low mu Low s",
     col = col.vec[1])

for(i in 2:nrow(AY_LmuLs.prop)){
  lines(y = AY_LmuLs.prop[i, ], x = 1:ncol(AY_LmuLs.prop), col = col.vec[i])
}

plot(AY_HmuHs.prop[1,], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions",
     main = "Achiasmatic Y High mu High s",
     col = col.vec[1])

for(i in 2:nrow(AY_HmuHs.prop)){
  lines(y = AY_HmuHs.prop[i, ], x = 1:ncol(AY_HmuHs.prop), col = col.vec[i])
}

plot(AY_HmuLs.prop[1,], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions",
     main = "Achiasmatic Y High mu Low s",
     col = col.vec[1])

for(i in 2:nrow(AY_HmuLs.prop)){
  lines(y = AY_HmuLs.prop[i, ], x = 1:ncol(AY_HmuLs.prop), col = col.vec[i])
}

plot(AY_LmuHs.prop[1,], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions",
     main = "Achiasmatic Y Low mu High s",
     col = col.vec[1])

for(i in 2:nrow(AY_LmuHs.prop)){
  lines(y = AY_LmuHs.prop[i, ], x = 1:ncol(AY_LmuHs.prop), col = col.vec[i])
}



############################### Data collection ################################
# All matrices start with the same skeleton:
# col = generation
# row = simulation
num.all.A1.X <- all.A1.X <- num.all.A2.X <- all.A2.X <- 
  num.all.A1.Y <- all.A1.Y <- num.all.A2.Y <- all.A2.Y <-
  num.achias.A1.X <- achias.A1.X <- num.achias.A2.X <- achias.A2.X <- 
  num.achias.A1.Y <- achias.A1.Y <- num.achias.A2.Y <- achias.A2.Y <-
  num.chias.A1.X <- chias.A1.X <- num.chias.A2.X <- chias.A2.X <- 
  num.chias.A1.Y <- chias.A1.Y <- num.chias.A2.Y <- chias.A2.Y <- matrix(nrow = num_sims, ncol = num_gens)

# Add data from each group of 4 rows, as a group of 4 rows represents a gen
cur.gen.rows <- 1:4
cur.sim <- 1
cur.gen <- 1

while(cur.sim <= num_sims){
  print(cur.gen)
  
  # All Sims, X
  # Number of X chromosomes in current generation
  num_X <- ncol(achiasX.res) + sum(achiasX.res[cur.gen.rows[2],] == 0) +
             ncol(chiasX.res) + sum(chiasX.res[cur.gen.rows[2],] == 0) 
  
  # If fusions were permitted only on X chromosomes, then any X/Y loci with 
  # fusions must be X chromosomes
  num.all.A1.X[cur.sim, cur.gen] <- sum(achiasX.res[cur.gen.rows[3],] == 1,
                                        achiasX.res[cur.gen.rows[4],] == 1,
                                        chiasX.res[cur.gen.rows[3],] == 1,
                                        chiasX.res[cur.gen.rows[4],] == 1)
  all.A1.X[cur.sim, cur.gen] <- num.all.A1.X[cur.sim, cur.gen]/num_X
  
  num.all.A2.X[cur.sim, cur.gen] <- sum(achiasX.res[cur.gen.rows[3],] == 2,
                                        achiasX.res[cur.gen.rows[4],] == 2,
                                        chiasX.res[cur.gen.rows[3],] == 2,
                                        chiasX.res[cur.gen.rows[4],] == 2)
  all.A2.X[cur.sim, cur.gen] <- num.all.A2.X[cur.sim, cur.gen]/num_X
  
  # Achiasmy, X
  num_X_achias <- ncol(achiasX.res) + sum(achiasX.res[cur.gen.rows[2],] == 0)
  
  num.achias.A1.X[cur.sim, cur.gen] <- sum(achiasX.res[cur.gen.rows[3],] == 1,
                                        achiasX.res[cur.gen.rows[4],] == 1)
  achias.A1.X[cur.sim, cur.gen] <- num.achias.A1.X[cur.sim, cur.gen]/num_X_achias
  
  num.achias.A2.X[cur.sim, cur.gen] <- sum(achiasX.res[cur.gen.rows[3],] == 2,
                                        achiasX.res[cur.gen.rows[4],] == 2)
  achias.A2.X[cur.sim, cur.gen] <- num.achias.A2.X[cur.sim, cur.gen]/num_X_achias
  
  # Chiasmy, X
  num_X_chias <- ncol(chiasX.res) + sum(chiasX.res[cur.gen.rows[2],] == 0)
  
  num.chias.A1.X[cur.sim, cur.gen] <- sum(chiasX.res[cur.gen.rows[3],] == 1,
                                           chiasX.res[cur.gen.rows[4],] == 1)
  chias.A1.X[cur.sim, cur.gen] <- num.chias.A1.X[cur.sim, cur.gen]/num_X_chias
  
  num.chias.A2.X[cur.sim, cur.gen] <- sum(chiasX.res[cur.gen.rows[3],] == 2,
                                           chiasX.res[cur.gen.rows[4],] == 2)
  chias.A2.X[cur.sim, cur.gen] <- num.chias.A2.X[cur.sim, cur.gen]/num_X_chias
  
  # All Sims, Y
  num_Y <- sum(achiasY.res[cur.gen.rows[2],] == 1, chiasY.res[cur.gen.rows[2],] == 1) 
  
  num.all.A1.Y[cur.sim, cur.gen] <- sum(achiasY.res[cur.gen.rows[3],] == 1,
                                        achiasY.res[cur.gen.rows[4],] == 1,
                                        chiasY.res[cur.gen.rows[3],] == 1,
                                        chiasY.res[cur.gen.rows[4],] == 1)
  all.A1.Y[cur.sim, cur.gen] <- num.all.A1.Y[cur.sim, cur.gen]/num_Y
  
  num.all.A2.Y[cur.sim, cur.gen] <- sum(achiasY.res[cur.gen.rows[3],] == 2,
                                        achiasY.res[cur.gen.rows[4],] == 2,
                                        chiasY.res[cur.gen.rows[3],] == 2,
                                        chiasY.res[cur.gen.rows[4],] == 2)
  all.A2.Y[cur.sim, cur.gen] <- num.all.A2.Y[cur.sim, cur.gen]/num_Y
    
  
  
  # Achiasmy, Y
  num_Y_achias <- sum(achiasY.res[cur.gen.rows[2],] == 1)
  
  num.achias.A1.Y[cur.sim, cur.gen] <- sum(achiasY.res[cur.gen.rows[3],] == 1,
                                           achiasY.res[cur.gen.rows[4],] == 1)
  achias.A1.Y[cur.sim, cur.gen] <- num.achias.A1.Y[cur.sim, cur.gen]/num_Y_achias
  
  num.achias.A2.Y[cur.sim, cur.gen] <- sum(achiasY.res[cur.gen.rows[3],] == 2,
                                           achiasY.res[cur.gen.rows[4],] == 2)
  achias.A2.Y[cur.sim, cur.gen] <- num.achias.A2.Y[cur.sim, cur.gen]/num_Y_achias
  
  # Chiasmy, Y
  num_Y_chias <- sum(chiasY.res[cur.gen.rows[2],] == 1)
  
  num.chias.A1.Y[cur.sim, cur.gen] <- sum(chiasY.res[cur.gen.rows[3],] == 1,
                                          chiasY.res[cur.gen.rows[4],] == 1)
  chias.A1.Y[cur.sim, cur.gen] <- num.chias.A1.Y[cur.sim, cur.gen]/num_Y_chias
  
  num.chias.A2.Y[cur.sim, cur.gen] <- sum(chiasY.res[cur.gen.rows[3],] == 2,
                                          chiasY.res[cur.gen.rows[4],] == 2)
  chias.A2.Y[cur.sim, cur.gen] <- num.chias.A2.Y[cur.sim, cur.gen]/num_Y_chias
  
  
  # If you've reached the last generation in this simulation, move to the next
  if(cur.gen == num_gens){
    cur.sim <- cur.sim + 1
    cur.gen <- 1
    cur.gen.rows <- cur.gen.rows + 4
    next
  }else{
    cur.gen <- cur.gen + 1
    cur.gen.rows <- cur.gen.rows + 4
    next
  }
}

################################## Question 2 ##################################
# Question 2: Do large and small autosomes differ in frequency of fusion reg-
# ardless of achiasmy status of species?

# Finding Q2: 
t.test(num.all.A1.X[,num_gens], num.all.A2.X[,num_gens])
t.test(num.all.A1.Y[,num_gens], num.all.A2.Y[,num_gens])

# Plot 2.1:
# 1 color per simulation
# X = Generations (100)
# Y = Frequency of small-X fusion across both chiasmatic and achiasmatic
#     conditions

plot(all.A1.X[1,], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Small Fusions on X Chr in Both A/C",
     col = col.vec[1])

for(i in 2:nrow(all.A1.X)){
  lines(y = all.A1.X[i, ], x = 1:ncol(all.A1.X), col = col.vec[i])
}

# Plot 2.2:
plot(all.A2.X[1,], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Large Fusions on X Chr in Both A/C",
     col = col.vec[1])

for(i in 2:nrow(all.A2.X)){
  lines(y = all.A2.X[i, ], x = 1:ncol(all.A2.X), col = col.vec[i])
}

# Plot 2.3:
plot(all.A1.Y[1,], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Small Fusions on Y Chr in Both A/C",
     col = col.vec[1])

for(i in 2:nrow(all.A1.Y)){
  lines(y = all.A1.Y[i, ], x = 1:ncol(all.A1.Y), col = col.vec[i])
}

# Plot 2.4:
plot(all.A2.Y[1,], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Large Fusions on Y Chr in Both A/C",
     col = col.vec[1])

for(i in 2:nrow(all.A2.Y)){
  lines(y = all.A2.Y[i, ], x = 1:ncol(all.A2.Y), col = col.vec[i])
}


################################## Question 3 ##################################
# Question 3: Does the proportion of small and large autosomal fusions differ 
# between chiasmatic and achiasmatic species?

################################## Question 4 ##################################
# Question 4: Does the proportion of small vs. large fusions differ between X
# and Y chromosomes?
# Finding: 

################################## Question 5 ##################################
# Question 5: Do achiasmatic species have lower fequency of large autosomal
# fusions compared to chiasmatic species?
# Finding: 

