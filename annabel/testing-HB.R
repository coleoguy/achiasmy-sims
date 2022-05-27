source("functions.R")

# Run simulation multiple times, storing output of each run as a new element in
# list
runs <- 5
pop_size <- 100
gen_no <- 300
mu <- 1^-6
fus.type <- "Y"
fus.large <- F
s <- .3
# Get probabilities of 0, 1, 2, and 3 mutations in a generation


pop <- GetPop(pop_size)

# Get a distribution of fitness effects
dfe <- GetDFE()
mon <- c()
for(gen in 1:gen_no){
  print(paste(c("Generation: ", gen), collapse = ""))
  # Mutate the starting population
  pop <- ActofGod(pop, dfe, fus.type, mu, fus.large)
  # Assess the fitness of the mutated population
  fits <- GetFit(pop, s)
  mon[gen] <- mean(fits)
  # Add loci and fitnesses to corresponding generation in output
  # Select parents based on the fitness
  parents <- Maury(pop, fits, length(pop))
  # Select gametes from these parents
  gametes <- MakeGametes(pop, parents, chiasm=T)
  # Breed the parents
  pop <- MiracleOfLife(gametes)
  # Move to next two generation's rows
}
