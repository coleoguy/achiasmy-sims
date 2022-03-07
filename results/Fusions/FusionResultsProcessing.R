setwd("/home/blackmonlab/Documents/Annabel/AchiasmyFusionSim/")

# Read in results for simulation under achiasmatic and chiasmatic conditions
achis.res <- readRDS("AchiasmaticResults.rds")
chis.res <- readRDS("ChiasmaticResults.rds")
library(viridisLite)
library(ggplot2)
col.vec <- viridis(length(achis.res))


# Question 1: Do fusions occur more often on X or Y?
# Prediction: Fusions more common on Y
# Conclusion: No difference in frequencies

# vectors in which to store frequencies of X & Y fusions in last generation of 
# each simulation
all.fus.X <- rep(NA, 1000)
all.fus.Y <- rep(NA, 1000)

# Iterate through each simulation and store the total frequency of X and Y 
# fusions at last generation
for(sim in 1:length(chis.res)){
  # Sum all X vs. Y chr in both simulation classes
  X.denom <- ncol(achis.res[[sim]]$SDR) + ncol(chis.res[[sim]]$SDR) + sum(!achis.res[[sim]]$SDR[nrow(achis.res[[sim]]$SDR),]) + sum(!chis.res[[sim]]$SDR[nrow(chis.res[[sim]]$SDR),])
  Y.denom <- sum(achis.res[[sim]]$SDR[nrow(achis.res[[sim]]$SDR),]) + sum(chis.res[[sim]]$SDR[nrow(chis.res[[sim]]$SDR),])
  
  # Sum all X and Y with fusions in both simulation classes
  X.num <- sum(achis.res[[sim]]$FusionLocus[nrow(achis.res[[sim]]$FusionLocus), (achis.res[[sim]]$SDR[nrow(achis.res[[sim]]$SDR),] == 0)] != 0) +
    sum(chis.res[[sim]]$FusionLocus[nrow(chis.res[[sim]]$FusionLocus), (chis.res[[sim]]$SDR[nrow(chis.res[[sim]]$SDR),] == 0)] != 0) +
    sum(achis.res[[sim]]$FusionLocus[nrow(achis.res[[sim]]$FusionLocus) - 1,] != 0) +
    sum(chis.res[[sim]]$FusionLocus[nrow(chis.res[[sim]]$FusionLocus) - 1,] != 0)
    
  Y.num <- sum(achis.res[[sim]]$FusionLocus[nrow(achis.res[[sim]]$FusionLocus), (achis.res[[sim]]$SDR[nrow(achis.res[[sim]]$SDR),] == 1)] != 0) +
    sum(chis.res[[sim]]$FusionLocus[nrow(chis.res[[sim]]$FusionLocus), (chis.res[[sim]]$SDR[nrow(chis.res[[sim]]$SDR),] == 1)] != 0)
  
  # Get frequencies of both X and Y fusions in both simulations
  all.fus.X[sim] <- X.num/X.denom
  all.fus.Y[sim] <- Y.num/Y.denom
}



################################ Data for Q2-Q5 ################################
# Matrices: Each datapoint is the total frequency of an allele at a generation
# row = sim
# col = gen

# All Sims, All Sex
num.all.A1.sex <- matrix(nrow = length(achis.res),
                         ncol = nrow(achis.res[[1]]$FusionLocus)/2)
all.A1.sex <- matrix(nrow = length(achis.res),
                     ncol = nrow(achis.res[[1]]$FusionLocus)/2)
num.all.A2.sex <- matrix(nrow = length(achis.res),
                         ncol = nrow(achis.res[[1]]$FusionLocus)/2)
all.A2.sex <- matrix(nrow = length(achis.res),
                     ncol = nrow(achis.res[[1]]$FusionLocus)/2)

# Achiasmatic Sims, All Sex
num.achias.A1.sex <- matrix(nrow = length(achis.res),
                         ncol = nrow(achis.res[[1]]$FusionLocus)/2)
achias.A1.sex <- matrix(nrow = length(achis.res),
                        ncol = nrow(achis.res[[1]]$FusionLocus)/2)
num.achias.A2.sex <- matrix(nrow = length(achis.res),
                            ncol = nrow(achis.res[[1]]$FusionLocus)/2)
achias.A2.sex <- matrix(nrow = length(achis.res),
                        ncol = nrow(achis.res[[1]]$FusionLocus)/2)

# Chiasmatic Sims, All Sex
num.chias.A1.sex <- matrix(nrow = length(chis.res),
                            ncol = nrow(chis.res[[1]]$FusionLocus)/2)
chias.A1.sex <- matrix(nrow = length(chis.res),
                        ncol = nrow(chis.res[[1]]$FusionLocus)/2)
num.chias.A2.sex <- matrix(nrow = length(chis.res),
                            ncol = nrow(chis.res[[1]]$FusionLocus)/2)
chias.A2.sex <- matrix(nrow = length(chis.res),
                        ncol = nrow(chis.res[[1]]$FusionLocus)/2)

# All Sims, X
num.all.A1.X <- matrix(nrow = length(achis.res),
                         ncol = nrow(achis.res[[1]]$FusionLocus)/2)
all.A1.X <- matrix(nrow = length(achis.res),
                     ncol = nrow(achis.res[[1]]$FusionLocus)/2)
num.all.A2.X <- matrix(nrow = length(achis.res),
                         ncol = nrow(achis.res[[1]]$FusionLocus)/2)
all.A2.X <- matrix(nrow = length(achis.res),
                     ncol = nrow(achis.res[[1]]$FusionLocus)/2)

# All Sims, Y
num.all.A1.Y <- matrix(nrow = length(achis.res),
                       ncol = nrow(achis.res[[1]]$FusionLocus)/2)
all.A1.Y <- matrix(nrow = length(achis.res),
                   ncol = nrow(achis.res[[1]]$FusionLocus)/2)
num.all.A2.Y <- matrix(nrow = length(achis.res),
                       ncol = nrow(achis.res[[1]]$FusionLocus)/2)
all.A2.Y <- matrix(nrow = length(achis.res),
                   ncol = nrow(achis.res[[1]]$FusionLocus)/2)

# Chiasmatic vs. Achiasmatic A2 on Y
num.achias.A2.Y <- matrix(nrow = length(achis.res),
                         ncol = nrow(achis.res[[1]]$FusionLocus)/2)
achias.A2.Y <- matrix(nrow = length(achis.res),
                     ncol = nrow(achis.res[[1]]$FusionLocus)/2)
num.chias.A2.Y <- matrix(nrow = length(chis.res),
                           ncol = nrow(chis.res[[1]]$FusionLocus)/2)
chias.A2.Y <- matrix(nrow = length(chis.res),
                       ncol = nrow(chis.res[[1]]$FusionLocus)/2)

for(sim in 1:nrow(all.A1.sex)){ 
  gen.pair <- 1:2
  for(gen in 1:(nrow(achis.res[[sim]]$FusionLocus)/2)){
    # All Sims, All Sex
    num.all.A1.sex[sim,gen] <- sum(sum(achis.res[[sim]]$FusionLocus[gen.pair,] == 1) +
                                     sum(chis.res[[sim]]$FusionLocus[gen.pair,] == 1))
    all.A1.sex[sim,gen] <- num.all.A1.sex[sim,gen]/(ncol(achis.res[[sim]]$FusionLocus)*4)
    num.all.A2.sex[sim,gen] <- sum(sum(achis.res[[sim]]$FusionLocus[gen.pair,] == 2) +
                                     sum(chis.res[[sim]]$FusionLocus[gen.pair,] == 2))
    all.A2.sex[sim,gen] <- num.all.A2.sex[sim,gen]/(ncol(achis.res[[sim]]$FusionLocus)*4)
    
    # Achiasmatic Sims, All Sex
    num.achias.A1.sex[sim,gen] <- sum(achis.res[[sim]]$FusionLocus[gen.pair,] == 1)
    achias.A1.sex[sim,gen] <- num.achias.A1.sex[sim,gen]/(ncol(achis.res[[sim]]$FusionLocus)*2)
    num.achias.A2.sex[sim,gen] <- sum(achis.res[[sim]]$FusionLocus[gen.pair,] == 2)
    achias.A2.sex[sim,gen] <- num.achias.A2.sex[sim,gen]/(ncol(achis.res[[sim]]$FusionLocus)*2)
    
    # Chiasmatic Sims, All Sex
    num.chias.A1.sex[sim,gen] <- sum(chis.res[[sim]]$FusionLocus[gen.pair,] == 1)
    chias.A1.sex[sim,gen] <- num.chias.A1.sex[sim,gen]/(ncol(chis.res[[sim]]$FusionLocus)*2)
    num.chias.A2.sex[sim,gen] <- sum(chis.res[[sim]]$FusionLocus[gen.pair,] == 2)
    chias.A2.sex[sim,gen] <- num.chias.A2.sex[sim,gen]/(ncol(chis.res[[sim]]$FusionLocus)*2)
    
    # All Sims, X
    X.denom <- ncol(achis.res[[sim]]$SDR) + ncol(chis.res[[sim]]$SDR) + sum(!achis.res[[sim]]$SDR[nrow(achis.res[[sim]]$SDR),]) + sum(!chis.res[[sim]]$SDR[nrow(chis.res[[sim]]$SDR),])
    num.all.A1.X[sim,gen] <-sum(achis.res[[sim]]$FusionLocus[gen.pair[1],] == 1) +
      sum(chis.res[[sim]]$FusionLocus[gen.pair[1],] == 1) +
      sum(achis.res[[sim]]$FusionLocus[gen.pair[2], (achis.res[[sim]]$SDR[gen.pair[2],] == 0)] == 1) +
      sum(chis.res[[sim]]$FusionLocus[gen.pair[2], (chis.res[[sim]]$SDR[gen.pair[2],] == 0)] == 1)
    all.A1.X[sim,gen] <- num.all.A1.X[sim,gen]/X.denom
    num.all.A2.X[sim,gen] <-sum(achis.res[[sim]]$FusionLocus[gen.pair[1],] == 2) +
      sum(chis.res[[sim]]$FusionLocus[gen.pair[1],] == 2) +
      sum(achis.res[[sim]]$FusionLocus[gen.pair[2], (achis.res[[sim]]$SDR[gen.pair[2],] == 0)] == 2) +
      sum(chis.res[[sim]]$FusionLocus[gen.pair[2], (chis.res[[sim]]$SDR[gen.pair[2],] == 0)] == 2)
    all.A2.X[sim,gen] <- num.all.A2.X[sim,gen]/X.denom
    
    # All Sims, Y
    Y.denom <- sum(achis.res[[sim]]$SDR[nrow(achis.res[[sim]]$SDR),]) + sum(chis.res[[sim]]$SDR[nrow(chis.res[[sim]]$SDR),])
    num.all.A1.Y[sim,gen] <- sum(achis.res[[sim]]$FusionLocus[gen.pair[2], (achis.res[[sim]]$SDR[gen.pair[2],] )] == 1) +
      sum(chis.res[[sim]]$FusionLocus[gen.pair[2], (chis.res[[sim]]$SDR[gen.pair[2],] )] == 1)
    all.A1.Y[sim,gen] <- num.all.A1.Y[sim,gen]/Y.denom
    num.all.A2.Y[sim,gen] <-sum(achis.res[[sim]]$FusionLocus[gen.pair[2], (achis.res[[sim]]$SDR[gen.pair[2],] )] == 2) +
      sum(chis.res[[sim]]$FusionLocus[gen.pair[2], (chis.res[[sim]]$SDR[gen.pair[2],] )] == 2)
    all.A2.Y[sim,gen] <- num.all.A2.Y[sim,gen]/Y.denom
    
    # Chiasmatic vs. Achiasmatic A2 on Y
    Y.denom.achias <- sum(achis.res[[sim]]$SDR[nrow(achis.res[[sim]]$SDR),])
    num.achias.A2.Y[sim,gen] <-sum(achis.res[[sim]]$FusionLocus[gen.pair[2], (achis.res[[sim]]$SDR[gen.pair[2],] )] == 2)
    achias.A2.Y[sim,gen] <- num.achias.A2.Y[sim,gen]/Y.denom.achias
    
    Y.denom.chias <- sum(chis.res[[sim]]$SDR[nrow(chis.res[[sim]]$SDR),])
    num.chias.A2.Y[sim,gen] <- sum(chis.res[[sim]]$FusionLocus[gen.pair[2], (chis.res[[sim]]$SDR[gen.pair[2],] )] == 2)
    chias.A2.Y[sim,gen] <- num.chias.A2.Y[sim,gen]/Y.denom.chias
    
    gen.pair <- gen.pair + 2
  }
}

################################## Question 2 ##################################
# Question 2: Do large and small autosomes differ in frequency of fusion reg-
# ardless of achiasmy status of species?
# Finding 2: Small but significantly greater proportion of small fusions
t.test(num.all.A1.sex[,100], num.all.A2.sex[,100])

# Plot 2.1:
# 1 color per simulation
# X = Generations (100)
# Y = Frequency of Autosome1-Sex fusion across both chiasmatic and achiasmatic
#     condition

plot(all.A1.sex[1,], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency A1 Fusion across All Sex Chr in Both A/C",
     col = col.vec[1])

for(i in 2:nrow(all.A1.sex)){
  lines(y = all.A1.sex[i, ], x = 1:ncol(all.A1.sex), col = col.vec[i])
}

# Plot 2.2:
# 1 color per simulation
# X = Generations (100)
# Y = Frequency of Autosome2-Sex fusion across both chiasmatic and achiasmatic
#     conditions

plot(all.A2.sex[1,], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency A2 Fusion across All Sex Chr in Both A/C",
     col = col.vec[1])

for(i in 2:nrow(all.A2.sex)){
  lines(y = all.A2.sex[i, ], x = 1:ncol(all.A2.sex), col = col.vec[i])
}


################################## Question 3 ##################################
# Question 3: Does the proportion of small and large autosomal fusions differ 
# between chiasmatic and achiasmatic species?
# Finding: Proportion of small:large NOT significantly greater in achiasmatic 
# than in chiasmatic
t.test(num.achias.A1.sex[,100]/num.achias.A2.sex[,100], num.chias.A1.sex[,100]/num.chias.A2.sex[,100])

# Finding: Total number of small fusions NOT significantly greater in achiasmatic
t.test(num.achias.A1.sex[,100], num.chias.A1.sex[,100])

# Plot 3.1
# 1 color per simulation
# Y = Frequency of small autosomal fusion in achiasmatic simulation
# X = Generation
plot(achias.A1.sex[1,], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency A2 Fusion Across All Sex Chr in Achiasmatic",
     col = col.vec[1])

for(i in 2:nrow(achias.A1.sex)){
  lines(y = achias.A1.sex[i, ], x = 1:ncol(achias.A1.sex), col = col.vec[i])
}

# Plot 3.2
# 1 color per simulation
# Y = Frequency of large autosomal fusion in achiasmatic simulation
# X = Generation
plot(achias.A2.sex[1,], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency A2 Fusion Across All Sex Chr in Achiasmatic",
     col = col.vec[1])

for(i in 2:nrow(achias.A2.sex)){
  lines(y = achias.A2.sex[i, ], x = 1:ncol(achias.A2.sex), col = col.vec[i])
}

# Plot 3.3
# 1 color per simulation
# Y = Frequency of small autosomal fusion in chiasmatic simulation
# X = Generation
plot(chias.A1.sex[1,], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency A2 Fusion Across All Sex Chr in chiasmatic",
     col = col.vec[1])

for(i in 2:nrow(chias.A1.sex)){
  lines(y = chias.A1.sex[i, ], x = 1:ncol(chias.A1.sex), col = col.vec[i])
}

# Plot 3.4
# 1 color per simulation
# Y = Frequency of large autosomal fusion in chiasmatic simulation
# X = Generation
plot(chias.A2.sex[1,], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency A2 Fusion Across All Sex Chr in chiasmatic",
     col = col.vec[1])

for(i in 2:nrow(chias.A2.sex)){
  lines(y = chias.A2.sex[i, ], x = 1:ncol(chias.A2.sex), col = col.vec[i])
}

################################## Question 4 ##################################
# Question 3: Does the proportion of small vs. large fusions differ between X
# and Y chromosomes?
# Finding: Proportion of small to large fusions is significantly higher on Y 
# than on X chromosomes

# Replace 0s with 1s so 0 is never in denominator
ed.all.A2.X <- all.A2.X
ed.num.all.A2.X <- num.all.A2.X
ed.all.A2.Y <- all.A2.Y
ed.num.all.A2.Y <- num.all.A2.Y
for(i in 1:nrow(num.all.A2.Y)){
  if(num.all.A2.Y[i,100] == 0){
    ed.num.all.A2.Y[i,100] <- 1
    ed.all.A2.Y[,100] <- 10^-100
  }
  if(num.all.A2.X[i,100] == 0){
    ed.num.all.A2.X[i,100] <- 1
    ed.all.A2.X[,100] <- 10^-100
  }
}

Y.A1_A2 <- num.all.A1.Y[,100]/ed.num.all.A2.Y[,100]
X.A1_A2 <- num.all.A1.X[,100]/ed.num.all.A2.X[,100]
Y.Y.t <- t.test(Y.A1_A2, X.A1_A2)

# Question 5: Do achiasmatic species have lower fequency of large autosomal
# fusions compared to chiasmatic species?
# Finding: Achiasmy does NOT impact the number of large fusions to Y chromosomes
t.test(num.achias.A2.Y, num.chias.A2.Y)

