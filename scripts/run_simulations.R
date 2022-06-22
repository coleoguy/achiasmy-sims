source("functions.R")

s <- 0.3
pop_size <- 1000
chiasm <- T
fus.type <- "Y"
fus.large <- F
mu <- 10^-9

results <- Evolve(num_sims, pop_size, gen_no, s, chiasm, fus.type, mu, fus.large)
