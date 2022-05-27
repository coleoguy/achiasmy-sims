source("functions.R")

# Run simulation multiple times, storing output of each run as a new element in
# list
runs <- 5
pop_size <- 100
gen_no <- 1000
mu <- 0.00000001
fus.type <- "Y"
fus.large <- F
s <- 1
# Get probabilities of 0, 1, 2, and 3 mutations in a generation


pop <- GetPop(pop_size)

# Get a distribution of fitness effects
dfe <- GetDFE()

# Initialize list of matrices
# For each matrix...
# Rows = Generations
# Columns = Individuals
results <- matrix(nrow = 4*gen_no, ncol = pop_size,
                  dimnames = list(
                    rep(c("X-SDR", "X/Y-SDR", "X-FuseLocus", "X/Y-FuseLocus"), 
                        gen_no), 1:pop_size))

# Create iterator for the FOUR rows corresponding to the current generation
gen_rows <- 1:4

# For each generation in "gen_no"...
for(gen in 1:gen_no){
  print(paste(c("Generation: ", gen), collapse = ""))
  
  # Mutate the starting population
  pop <- ActofGod(pop, dfe, fus.type, mu.table, fus.large)
  
  # Assess the fitness of the mutated population
  fits <- GetFit(pop, s)
  # Add loci and fitnesses to corresponding generation in output
  for(i in 1:pop_size){
    results[gen_rows, i] <- c(pop[[i]][,13], pop[[i]][,25])
  }
  
  # Select parents based on the fitness
  parents <- Maury(pop, fits, length(pop))
  
  # Select gametes from these parents
  gametes <- MakeGametes(pop, parents, chiasm=T)
  
  # Breed the parents
  pop <- MiracleOfLife(gametes)
  
  # Move to next two generation's rows
  gen_rows <- gen_rows + 4
}
