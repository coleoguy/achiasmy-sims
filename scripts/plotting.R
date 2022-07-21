# July 2022
# Annabel Perry
# This script plots all simulations from a single set of parameters and 
# tabulates the proportion of simulations fixed for fusions at the final 
# generation
library(viridis)

# set in parameters from the run of simulations
num_sims <- 1000
num_gen <- 1000
s_coeffs <- (0:10)/10

# There are 8 possible combinations of parameters because there are 2 values 
# each for "chiasm", "fus.type", and "fus.size"
num_conditions <- 8

# One line for each simulation
col.vec <- viridis(num_sims)
# Create dataframe which stores, for each combination of parameters, the 
# frequency of simulations which FIXED
Fixations <- data.frame(
  Meiotic_Type = rep(c("Chiasmatic", "Achiasmatic"), each = 20),
  Sex_Chr = rep(rep(c("Y", "X"), each = 10), times = 2),
  Autosome_Size = rep(c("Small", "Large"), times = 20),
  s = rep(rep(s_coeffs, each = 2), times = 20),
  Frequency_Fixations = rep(NA, times = 40)
)

# Plot and tabulate each possible combination of parameters. 
for(i in 1:length(s_coeffs)){
  s <- s_coeffs[i]
  for(cond in 5:num_conditions){
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
    
    # Retrieve filename based on current parameters
    filename <- paste(c("../results/", 
                        meiotic_type[1], fus.type, size[1], "_s=", s, 
                        "_", 
                        #TODO: add date code,
                        ".rds"), collapse = "")
    print(filename)
    # read in the results of a run of simulations as a list where each entry is 
    # a simulation and each element of the entry is the frequency of fusions at 
    # the end of a single generation
    res <- readRDS(filename)
    
    # Add the frequency of simulations which FIXED for fusions to the 
    # corresponding row of the dataframe
    fixation_count <- 0
    for(i in 1:length(res)){
      if(res[[i]][gen_no] == 1){
        fixation_count <- fixation_count + 1
      }
    }
    Fixations$Frequency_Fixations[Fixations$Meiotic_Type == meiotic_type[2] &
                                  Fixations$Sex_Chr == fus.type &
                                  Fixations$Autosome_Size == size[2] &
                                  Fixations$s == s] <- fixation_count/num_sims

    plot(res[[1]], type = "l", ylim=c(0,1), 
         xlab = "Generation", ylab = "Frequency Fusions", 
         main = paste(c(meiotic_type[2], " ", fus.type, "-", size[2]),
                      collapse = ""),
         col = col.vec[1])
  
    for(i in 2:length(res)){
      lines(y = res[[i]], x = 1:num_gen, col = col.vec[i])
    }
  }
}


