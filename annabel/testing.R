source("functions.R")

num_sims <- 2
pop_size <- 1000
gen_no <- 5
s <- 0.3
chiasm <- T
fus.type <- "Y"
mu <- 0.5
fus.large <- F


results <- Evolve(num_sims, pop_size, gen_no, s, chiasm, fus.type, mu, fus.large)