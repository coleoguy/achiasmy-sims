############ Question 1 is Deprecated (Added "fus.type" parameter) #############
# Question 1: Do fusions occur more often on X or Y?
# Prediction: Fusions more common on Y
# Conclusion from Experiment 1: No difference in frequencies

# vectors in which to store frequencies of X & Y fusions in last generation of 
# each simulation
all.fus.X <- rep(NA, 1000)
all.fus.Y <- rep(NA, 1000)

# Iterate through each simulation and store the total frequency of X and Y 
# fusions at last generation
for(sim in 1:length(chiasX.res)){
  # Sum all X vs. Y chr in both simulation classes
  X.denom <- ncol(achiasX.res[[sim]]$SDR) + ncol(chiasX.res[[sim]]$SDR) + sum(!achiasX.res[[sim]]$SDR[nrow(achiasX.res[[sim]]$SDR),]) + sum(!chiasX.res[[sim]]$SDR[nrow(chiasX.res[[sim]]$SDR),])
  Y.denom <- sum(achiasX.res[[sim]]$SDR[nrow(achiasX.res[[sim]]$SDR),]) + sum(chiasX.res[[sim]]$SDR[nrow(chiasX.res[[sim]]$SDR),])
  
  # Sum all X and Y with fusions in both simulation classes
  X.num <- sum(achiasX.res[[sim]]$FusionLocus[nrow(achiasX.res[[sim]]$FusionLocus), (achiasX.res[[sim]]$SDR[nrow(achiasX.res[[sim]]$SDR),] == 0)] != 0) +
    sum(chiasX.res[[sim]]$FusionLocus[nrow(chiasX.res[[sim]]$FusionLocus), (chiasX.res[[sim]]$SDR[nrow(chiasX.res[[sim]]$SDR),] == 0)] != 0) +
    sum(achiasX.res[[sim]]$FusionLocus[nrow(achiasX.res[[sim]]$FusionLocus) - 1,] != 0) +
    sum(chiasX.res[[sim]]$FusionLocus[nrow(chiasX.res[[sim]]$FusionLocus) - 1,] != 0)
  
  Y.num <- sum(achiasX.res[[sim]]$FusionLocus[nrow(achiasX.res[[sim]]$FusionLocus), (achiasX.res[[sim]]$SDR[nrow(achiasX.res[[sim]]$SDR),] == 1)] != 0) +
    sum(chiasX.res[[sim]]$FusionLocus[nrow(chiasX.res[[sim]]$FusionLocus), (chiasX.res[[sim]]$SDR[nrow(chiasX.res[[sim]]$SDR),] == 1)] != 0)
  
  # Get frequencies of both X and Y fusions in both simulations
  all.fus.X[sim] <- X.num/X.denom
  all.fus.Y[sim] <- Y.num/Y.denom
}



##################### Deprecated Data Creation for Q2-Q5 #######################
# Matrices: Each datapoint is the total frequency of an allele at a generation
# row = sim
# col = gen

# All Sims, All Sex
num.all.A1.sex <- matrix(nrow = length(achiasX.res),
                         ncol = nrow(achiasX.res[[1]]$FusionLocus)/2)
all.A1.sex <- matrix(nrow = length(achiasX.res),
                     ncol = nrow(achiasX.res[[1]]$FusionLocus)/2)
num.all.A2.sex <- matrix(nrow = length(achiasX.res),
                         ncol = nrow(achiasX.res[[1]]$FusionLocus)/2)
all.A2.sex <- matrix(nrow = length(achiasX.res),
                     ncol = nrow(achiasX.res[[1]]$FusionLocus)/2)

# Achiasmatic Sims, All Sex
num.achias.A1.sex <- matrix(nrow = length(achiasX.res),
                            ncol = nrow(achiasX.res[[1]]$FusionLocus)/2)
achias.A1.sex <- matrix(nrow = length(achiasX.res),
                        ncol = nrow(achiasX.res[[1]]$FusionLocus)/2)
num.achias.A2.sex <- matrix(nrow = length(achiasX.res),
                            ncol = nrow(achiasX.res[[1]]$FusionLocus)/2)
achias.A2.sex <- matrix(nrow = length(achiasX.res),
                        ncol = nrow(achiasX.res[[1]]$FusionLocus)/2)

# Chiasmatic Sims, All Sex
num.chias.A1.sex <- matrix(nrow = length(chiasX.res),
                           ncol = nrow(chiasX.res[[1]]$FusionLocus)/2)
chias.A1.sex <- matrix(nrow = length(chiasX.res),
                       ncol = nrow(chiasX.res[[1]]$FusionLocus)/2)
num.chias.A2.sex <- matrix(nrow = length(chiasX.res),
                           ncol = nrow(chiasX.res[[1]]$FusionLocus)/2)
chias.A2.sex <- matrix(nrow = length(chiasX.res),
                       ncol = nrow(chiasX.res[[1]]$FusionLocus)/2)

# All Sims, X
num.all.A1.X <- matrix(nrow = length(achiasX.res),
                       ncol = nrow(achiasX.res[[1]]$FusionLocus)/2)
all.A1.X <- matrix(nrow = length(achiasX.res),
                   ncol = nrow(achiasX.res[[1]]$FusionLocus)/2)
num.all.A2.X <- matrix(nrow = length(achiasX.res),
                       ncol = nrow(achiasX.res[[1]]$FusionLocus)/2)
all.A2.X <- matrix(nrow = length(achiasX.res),
                   ncol = nrow(achiasX.res[[1]]$FusionLocus)/2)

# All Sims, Y
num.all.A1.Y <- matrix(nrow = length(achiasX.res),
                       ncol = nrow(achiasX.res[[1]]$FusionLocus)/2)
all.A1.Y <- matrix(nrow = length(achiasX.res),
                   ncol = nrow(achiasX.res[[1]]$FusionLocus)/2)
num.all.A2.Y <- matrix(nrow = length(achiasX.res),
                       ncol = nrow(achiasX.res[[1]]$FusionLocus)/2)
all.A2.Y <- matrix(nrow = length(achiasX.res),
                   ncol = nrow(achiasX.res[[1]]$FusionLocus)/2)

# Chiasmatic vs. Achiasmatic A2 on Y
num.achias.A2.Y <- matrix(nrow = length(achiasX.res),
                          ncol = nrow(achiasX.res[[1]]$FusionLocus)/2)
achias.A2.Y <- matrix(nrow = length(achiasX.res),
                      ncol = nrow(achiasX.res[[1]]$FusionLocus)/2)
num.chias.A2.Y <- matrix(nrow = length(chiasX.res),
                         ncol = nrow(chiasX.res[[1]]$FusionLocus)/2)
chias.A2.Y <- matrix(nrow = length(chiasX.res),
                     ncol = nrow(chiasX.res[[1]]$FusionLocus)/2)

for(sim in 1:nrow(all.A1.sex)){ 
  gen.pair <- 1:2
  for(gen in 1:(nrow(achiasX.res[[sim]]$FusionLocus)/2)){
    # All Sims, All Sex
    num.all.A1.sex[sim,gen] <- sum(sum(achiasX.res[[sim]]$FusionLocus[gen.pair,] == 1) +
                                     sum(chiasX.res[[sim]]$FusionLocus[gen.pair,] == 1))
    all.A1.sex[sim,gen] <- num.all.A1.sex[sim,gen]/(ncol(achiasX.res[[sim]]$FusionLocus)*4)
    num.all.A2.sex[sim,gen] <- sum(sum(achiasX.res[[sim]]$FusionLocus[gen.pair,] == 2) +
                                     sum(chiasX.res[[sim]]$FusionLocus[gen.pair,] == 2))
    all.A2.sex[sim,gen] <- num.all.A2.sex[sim,gen]/(ncol(achiasX.res[[sim]]$FusionLocus)*4)
    
    # Achiasmatic Sims, All Sex
    num.achias.A1.sex[sim,gen] <- sum(achiasX.res[[sim]]$FusionLocus[gen.pair,] == 1)
    achias.A1.sex[sim,gen] <- num.achias.A1.sex[sim,gen]/(ncol(achiasX.res[[sim]]$FusionLocus)*2)
    num.achias.A2.sex[sim,gen] <- sum(achiasX.res[[sim]]$FusionLocus[gen.pair,] == 2)
    achias.A2.sex[sim,gen] <- num.achias.A2.sex[sim,gen]/(ncol(achiasX.res[[sim]]$FusionLocus)*2)
    
    # Chiasmatic Sims, All Sex
    num.chias.A1.sex[sim,gen] <- sum(chiasX.res[[sim]]$FusionLocus[gen.pair,] == 1)
    chias.A1.sex[sim,gen] <- num.chias.A1.sex[sim,gen]/(ncol(chiasX.res[[sim]]$FusionLocus)*2)
    num.chias.A2.sex[sim,gen] <- sum(chiasX.res[[sim]]$FusionLocus[gen.pair,] == 2)
    chias.A2.sex[sim,gen] <- num.chias.A2.sex[sim,gen]/(ncol(chiasX.res[[sim]]$FusionLocus)*2)
    
    # All Sims, X
    X.denom <- ncol(achiasX.res[[sim]]$SDR) + ncol(chiasX.res[[sim]]$SDR) + sum(!achiasX.res[[sim]]$SDR[nrow(achiasX.res[[sim]]$SDR),]) + sum(!chiasX.res[[sim]]$SDR[nrow(chiasX.res[[sim]]$SDR),])
    num.all.A1.X[sim,gen] <-sum(achiasX.res[[sim]]$FusionLocus[gen.pair[1],] == 1) +
      sum(chiasX.res[[sim]]$FusionLocus[gen.pair[1],] == 1) +
      sum(achiasX.res[[sim]]$FusionLocus[gen.pair[2], (achiasX.res[[sim]]$SDR[gen.pair[2],] == 0)] == 1) +
      sum(chiasX.res[[sim]]$FusionLocus[gen.pair[2], (chiasX.res[[sim]]$SDR[gen.pair[2],] == 0)] == 1)
    all.A1.X[sim,gen] <- num.all.A1.X[sim,gen]/X.denom
    num.all.A2.X[sim,gen] <-sum(achiasX.res[[sim]]$FusionLocus[gen.pair[1],] == 2) +
      sum(chiasX.res[[sim]]$FusionLocus[gen.pair[1],] == 2) +
      sum(achiasX.res[[sim]]$FusionLocus[gen.pair[2], (achiasX.res[[sim]]$SDR[gen.pair[2],] == 0)] == 2) +
      sum(chiasX.res[[sim]]$FusionLocus[gen.pair[2], (chiasX.res[[sim]]$SDR[gen.pair[2],] == 0)] == 2)
    all.A2.X[sim,gen] <- num.all.A2.X[sim,gen]/X.denom
    
    # All Sims, Y
    Y.denom <- sum(achiasX.res[[sim]]$SDR[nrow(achiasX.res[[sim]]$SDR),]) + sum(chiasX.res[[sim]]$SDR[nrow(chiasX.res[[sim]]$SDR),])
    num.all.A1.Y[sim,gen] <- sum(achiasX.res[[sim]]$FusionLocus[gen.pair[2], (achiasX.res[[sim]]$SDR[gen.pair[2],] )] == 1) +
      sum(chiasX.res[[sim]]$FusionLocus[gen.pair[2], (chiasX.res[[sim]]$SDR[gen.pair[2],] )] == 1)
    all.A1.Y[sim,gen] <- num.all.A1.Y[sim,gen]/Y.denom
    num.all.A2.Y[sim,gen] <-sum(achiasX.res[[sim]]$FusionLocus[gen.pair[2], (achiasX.res[[sim]]$SDR[gen.pair[2],] )] == 2) +
      sum(chiasX.res[[sim]]$FusionLocus[gen.pair[2], (chiasX.res[[sim]]$SDR[gen.pair[2],] )] == 2)
    all.A2.Y[sim,gen] <- num.all.A2.Y[sim,gen]/Y.denom
    
    # Chiasmatic vs. Achiasmatic A2 on Y
    Y.denom.achias <- sum(achiasX.res[[sim]]$SDR[nrow(achiasX.res[[sim]]$SDR),])
    num.achias.A2.Y[sim,gen] <-sum(achiasX.res[[sim]]$FusionLocus[gen.pair[2], (achiasX.res[[sim]]$SDR[gen.pair[2],] )] == 2)
    achias.A2.Y[sim,gen] <- num.achias.A2.Y[sim,gen]/Y.denom.achias
    
    Y.denom.chias <- sum(chiasX.res[[sim]]$SDR[nrow(chiasX.res[[sim]]$SDR),])
    num.chias.A2.Y[sim,gen] <- sum(chiasX.res[[sim]]$FusionLocus[gen.pair[2], (chiasX.res[[sim]]$SDR[gen.pair[2],] )] == 2)
    chias.A2.Y[sim,gen] <- num.chias.A2.Y[sim,gen]/Y.denom.chias
    
    gen.pair <- gen.pair + 2
  }
}


library(viridisLite)
library(ggplot2)

setwd("/home/blackmonlab/Documents/Annabel/AchiasmyFusionSim/")

# Read in results for simulation under achiasmatic and chiasmatic conditions
achiasX.res <- readRDS("XAchiasmaticResults.rds")
chiasX.res <- readRDS("XChiasmaticResults.rds")
achiasY.res <- readRDS("YAchiasmaticResults.rds")
chiasY.res <- readRDS("YChiasmaticResults.rds")
col.vec <- viridis(length(achiasX.res))

############################## General Parameters ##############################
# Number of simulations must be known a priori
num_sims <- 2
# Number of generations is number of rows in a matrix divided by 4*num_sims
num_gens <- nrow(chiasY.res)/(4*num_sims)
# Number of individuals is number of columns in a matrix
num_indv <- ncol(chiasY.res)

################################## Question 2 ##################################
# Question 2: Do large and small autosomes differ in frequency of fusion reg-
# ardless of achiasmy status of species?

# Finding Q2 (Experiment 1): Small but significantly greater proportion of small fusions
t.test(num.all.A1.sex[,100], num.all.A2.sex[,100])

# TODO Finding Q2 (Experiment 2):


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

