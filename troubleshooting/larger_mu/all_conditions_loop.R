# July 2022
# Annabel Perry
# This script loops through the simulation with all possible combinations of
# meiotic & fusion parameters as well as all selection coefficients from 0:1

source("../../scripts/functions.R")

#################### Create Cluster for Parallel Processing #################### 
sim_cluster <- makeCluster(spec = 100, type = "SOCK")
registerDoSNOW(cl = sim_cluster)

################################# Define Parameters ############################ 
num_sims <- 1000
gen_no <- 1000
pop_size <- 1000     
# Run at s of 0.1, because at this selection coefficient SOME Y-S simulations
# will fix but not all will fix
s <- 0.1
# Empirically-derived mutation rate for Drosophila
mus <- c(5e-8, 5e-7, 5e-6)

# To start, run only on CYS and AYS
meiotic_types <- c("C", "A")
fus.type <- "Y"
fus.large <- F
size <- "S"
# Run first sim chiasmatic and second achiasmatic
chiasm <- T

# Run with each possible combination of mu and chiasmy
for(i in 1:length(mus)){
  mu <- mus[i]
  for(cond in 1:length(meiotic_types)){ 
    time.start <- Sys.time()
    meiotic_type <- meiotic_types[cond]
    
    # Prep filename based on current parameters 
    filename <- paste(c(meiotic_type, fus.type, size, "_mu=", mu,".rds"), 
                      collapse = "")
    print(filename)
    
    # Create a list where each entry is a vector of frequencies for each 
    # generation of a simulation
    fusion_frequencies <- vector(mode = "list", length = num_sims)
    # Run sims in parallel
    dfe <- GetDFE()
    fusion_frequencies <- foreach(sim = 1:num_sims) %dopar% {
      Evolve(pop_size, gen_no, s, chiasm, fus.type,  mu, fus.large, dfe)
    } 
    saveRDS(fusion_frequencies, filename)
    time.elapsed <- Sys.time() - time.start
    print(time.elapsed)
    # Run first sim chiasmatic and second achiasmatic
    chiasm <- F
  }
}

stopCluster(sim_cluster)
