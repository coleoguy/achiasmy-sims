# July 2022
# Annabel Perry
# This script plots the proportion of fixations across all generations for all
# simulations run under a single set of parameters. The purpose is to check data
library(viridis)

# set in parameters from the run of simulations
num_sims <- 1000
gen_no <- 1000
s_coeffs <- (0:5)/5

# There are 8 possible combinations of parameters because there are 2 values 
# each for "chiasm", "fus.type", and "fus.size"
num_conditions <- 8

# One line for each simulation
col.vec <- viridis(num_sims)

setwd("../results/")
filenames <- matrix(c(
  "CYS_s=1_2022-07-19.rds", "AXL_s=0.2_2022-07-16.rds", 
  "AXL_s=0.4_2022-07-17.rds", "AXL_s=0.6_2022-07-18.rds",
  "AXL_s=0.8_2022-07-19.rds", "AXL_s=0_2022-07-15.rds",
  "AXL_s=1_2022-07-20.rds", "AXS_s=0.2_2022-07-16.rds",
  "AXS_s=0.4_2022-07-17.rds", "AXS_s=0.6_2022-07-18.rds",
  "AXS_s=0.8_2022-07-19.rds", "AXS_s=0_2022-07-15.rds",
  "AXS_s=1_2022-07-20.rds", "AYL_s=0.2_2022-07-16.rds",
  "AYL_s=0.4_2022-07-17.rds", "AYL_s=0.6_2022-07-18.rds",
  "AYL_s=0.8_2022-07-19.rds", "AYL_s=0_2022-07-15.rds",
  "AYL_s=1_2022-07-20.rds", "AYS_s=0.2_2022-07-16.rds",
  "AYS_s=0.4_2022-07-17.rds", "AYS_s=0.6_2022-07-18.rds",
  "AYS_s=0.8_2022-07-19.rds", "AYS_s=0_2022-07-15.rds",
  "AYS_s=1_2022-07-20.rds", "CXL_s=0.2_2022-07-15.rds",
  "CXL_s=0.4_2022-07-17.rds", "CXL_s=0.6_2022-07-18.rds",
  "CXL_s=0.8_2022-07-19.rds", "CXL_s=0_2022-07-14.rds", 
  "CXL_s=1_2022-07-20.rds", "CXS_s=0.2_2022-07-15.rds",
  "CXS_s=0.4_2022-07-16.rds", "CXS_s=0.6_2022-07-17.rds",
  "CXS_s=0.8_2022-07-19.rds", "CXS_s=0_2022-07-14.rds",
  "CXS_s=1_2022-07-20.rds", "CYL_s=0.2_2022-07-15.rds",
  "CYL_s=0.4_2022-07-16.rds", "CYL_s=0.6_2022-07-17.rds",
  "CYL_s=0.8_2022-07-18.rds", "CYL_s=0_2022-07-14.rds",
  "CYL_s=1_2022-07-19.rds", "CYS_s=0.2_2022-07-15.rds",
  "CYS_s=0.4_2022-07-16.rds", "CYS_s=0.6_2022-07-17.rds",
  "CYS_s=0.8_2022-07-18.rds", "CYS_s=0_2022-07-14.rds"))
# Plot each possible combination of parameters. 
for(i in 1:length(s_coeffs)){
  s <- s_coeffs[i]
  for(cond in 1:num_conditions){
    # Run first 4 sims chiasmatic and final four achiasmatic
    if(cond <= 4){
      meiotic_type <- c("C", "Chiasmatic")
    }else{
      meiotic_type <- c("A", "Achiasmatic")
    }
    
    # Every 2 runs, switch fusing from Y to X
    if(cond %in% c(1, 2, 5, 6)){
      fus.type <- "Y"
    }else{
      fus.type <- "X"
    }
    
    # Every other run, switch from small to large
    if(cond%%2){
      size <- c("S", "Small")
    }else{
      size <- c("L", "Large")
    }
    
    # Retrieve the name of the file which matches the current parameters
    cur_file <- filenames[apply(filenames,  MARGIN = 1, FUN = grepl,
                                pattern = paste(c(meiotic_type[1], fus.type, 
                                                  size[1], "_s=", s, "_"),
                                                collapse = ""))]
    print(cur_file)
    # read in the results of a run of simulations as a list where each entry is 
    # a simulation and each element of the entry is the frequency of fusions at 
    # the end of a single generation
    res <- readRDS(cur_file)
    
    # Plot frequency of fusions across all generations for each simulation
    plot(res[[1]], type = "l", ylim=c(0,1), 
         xlab = "Generation", ylab = "Frequency Fusions", 
         main = paste(c(meiotic_type[2], " ", fus.type, "-", size[2], ", s = ", s),
                      collapse = ""),
         col = col.vec[1])
  
    for(i in 2:length(res)){
      lines(y = res[[i]], x = 1:num_gen, col = col.vec[i])
    }
  }
}


