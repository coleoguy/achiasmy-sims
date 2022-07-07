# July 2022
# Annabel Perry
# This script plots all simulations from a single set of parameters
library(viridis)

# read in the results of a run of simulations as a list where each entry is a
# simulation and each element of the entry is the frequency of fusions at the 
# end of a single generation
res <- readRDS("CYS_s=1_2022-07-07.rds")

num_gen <- length(res[[1]])
num_sims <- length(res)
col.vec <- viridis(num_sims)

plot(res[[1]], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "Chiasmatic Y-Small: mu = 5e^(-9)",
     col = col.vec[1])

for(i in 2:length(res)){
  lines(y = res[[i]], x = 1:num_gen, col = col.vec[i])
}
