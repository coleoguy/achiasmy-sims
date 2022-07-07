# July 2022
# Annabel Perry
# This script loops through the simulation with all possible combinations of
# meiotic & fusion parameters as well as all selection coefficients from 0:1

source("functions.R")

#################### Create Cluster for Parallel Processing #################### 
nCores <- detectCores() - 1
sim_cluster <- makeCluster(spec = nCores, type = "SOCK")
registerDoSNOW(cl = sim_cluster)

################################# Define Parameters ############################ 
num_sims <- 2
gen_no <- 5
pop_size <- 100
# Selection coefficients ranging from 0 (negative control) to 1 (most extreme)
s_coeffs <- (0:10)/10
# Empirically-derived mutation rate for Drosophila
mu <- 5e-9
# There are 8 possible combinations of parameters because there are 2 values 
# each for "chiasm", "fus.type", and "fus.size"
num_conditions <- 8

# Run with each possible combination of parameters. 
for(i in 1:length(s_coeffs)){
  s <- s_coeffs[i]
  for(cond in 1:num_conditions){
    # Run first 4 sims chiasmatic and final four achiasmatic
    if(cond <= 4){
      chiasm <- T
      meiotic_type <- "C"
    }else{
      chiasm <- F
      meiotic_type <- "A"
    }
    
    # Every 2 runs, switch fusing from Y to X
    if(cond %in% c(1, 2, 5, 6)){
      fus.type <- "Y"
    }else{
      fus.type <- "X"
    }
    
    # Every other run, switch from small to large
    if(cond%%2){
      fus.large <- F
      size <- "S"
    }else{
      fus.large <- T
      size <- "L"
    }
    
    # Prep filename based on current parameters
    filename <- paste(c(meiotic_type, fus.type, size, "_s=", s, "_", 
                        as.character(Sys.Date()), ".rds"), collapse = "")
    print(filename)
    
    # Create a list where each entry is a vector of frequencies for each 
    # generation of a simulation
    fusion_frequencies <- vector(mode = "list", length = num_sims)
    # Run sims in parallel
    fusion_frequencies <- foreach(sim = 1:num_sims) %dopar% {
      print(paste(c("Simulation: ", sim), collapse = ""))
      Evolve(pop_size, gen_no, s, chiasm, fus.type, mu, fus.large)
    }
    
    saveRDS(fusion_frequencies, filename)
  }
}
