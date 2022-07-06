source("functions.R")

num_sims <- 1
gen_no <- 1000
s <- 0.3
pop_size <- 1000
chiasm <- T
fus.type <- "Y"
fus.large <- F
mu <- 10^-14

start.time <- proc.time()

results <- Evolve(num_sims, pop_size, gen_no, s, chiasm, fus.type, mu, fus.large)

total.time <- proc.time() - start.time

job_hours <- 16*10*(total.time[[1]])/(60*60)
