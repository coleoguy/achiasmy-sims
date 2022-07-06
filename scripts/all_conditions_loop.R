# July 2022
# Annabel Perry
# This script loops through the simulation with all possible combinations of
# meiotic & fusion parameters

source("functions.R")

num_sims <- 10
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
  }else{
    chiasm <- F
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
  }else{
    fus.large <- T
  }
  # Prep filename based on current parameters
  if(chiasm){
    meiotic_type <- "C"
  }else{
    meiotic_type <- "A"
  }
  if(fus.large){
    size <- "L"
  }else{
    size <- "S"
  }
  filename <- paste(c(meiotic_type, fus.type, size, "_", 
                      as.character(Sys.Date()), ".rds"), collapse = "")
  print(filename)
  
  # Run sims
  res <- Evolve(num_sims, pop_size, gen_no, s, chiasm, fus.type, mu, fus.large)
  
  saveRDS(res, filename)
}