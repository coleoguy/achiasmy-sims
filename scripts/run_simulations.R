source("functions.R")

num_sims <- 2
gen_no <- 2
s <- 0
pop_size <- 10
chiasm <- T
fus.type <- "Y"
fus.large <- F
mu <- 0

results <- Evolve(num_sims, pop_size, gen_no, s, chiasm, fus.type, mu, fus.large)
