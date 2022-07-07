# July 2022
# Annabel Perry
# This script loops through the simulation with all possible combinations of
# meiotic & fusion parameters

source("functions.R")

# Create a cluster to run each simulation simultaneously
sim_cluster <- makeCluster(spec = 100, type = "SOCK")
registerDoSNOW(cl = sim_cluster)

num_sims <- 2
gen_no <- 1000
pop_size <- 1000
s <- 0
mu <- 5e-9

# Run with each possible combination of parameters. There are 8 possible 
# combinations of parameters because there are 2 values each for "chiasm",
# "fus.type", and "fus.size"
num_conditions <- 8
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
  filename <- paste(c(meiotic_type, fus.type, size, "_", 
                      as.character(Sys.Date()), ".rds"), collapse = "")
  print(filename)
  
  # Create matrix where...
  # row = generation
  # column = simulation
  # cell = frequency of fusions
  fusion_frequencies <- matrix
  
  # Run sims in parallel and store the resulting fusion frequencies after each
  # generation as a column
  foreach(sim = 1:num_sims) %dopar% {
    col <- Evolve(pop_size, gen_no, s, chiasm, fus.type, mu, fus.large)
    fusion_frequencies[, sim] <- as.vector(col)
  }
  
  saveRDS(fusion_frequencies, filename)
}