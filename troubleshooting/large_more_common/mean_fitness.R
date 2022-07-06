# 06-24-2022
# Find mean fitness of large and small chromosomes after several rounds of 
# mutation (in absence of fusions and sexual antagonism) to see if large 
# chromosomes on average have higher fitness than small chromosomes such that 
# fusions to small chromosomes have greater risk of being weeded out prior to 
# fixation
source("../../scripts/functions.R")

# Modified version of GetFit() function which outputs fitnesses of GENERAL 
# FITNESS loci on both large and small chromosomes
GetFit_Sizes <- function(pop){
  # Initialize matrix of fitness values for the large and small chromosomes in
  # the population
  fits <- matrix(nrow = 2, ncol = length(pop), dimnames = list(c("Large","Small")))
  
  # Iterate through each individual in the population
  for(i in 1:length(pop)){
    # Retrieve averaged fitness of all GENERAL fitness loci on the large and on
    # the small chromosome
    fits[1, i] <- prod(colSums(pop[[i]][,c(51:54,56:100)])/2)
    fits[2, i] <- prod(colSums(pop[[i]][,c(26:27,29:50)])/2)
  }
  # Return vector of fitness values of each individual
  return(fits)
}

# Modified version of ActofGod() where no fusions are possible
ActofGod <- function(pop, dfe, mu){
  # Initialize a vector to store sites which receive mutations
  mut_sites <- rep(0, ncol(pop[[1]]))
  # Iterate through each individual
  for(i in 1:length(pop)){
    # Decide which sites in this individual, if any, get mutations
    mut_sites[-c(13,25,28,55)] <- sample(0:1, size = (length(mut_sites) - 4), prob = c(1-mu, mu), replace = T)
    # If the individual gets any mutations...
    if(sum(mut_sites)){
      # At each site to be mutated, replace the value of a random homolog with
      # a randomly-selected value from dfe
      
      # In case Heath asks why you are LOOPING through all mutated sites:
      # If you were to instead do this...
      # pop[[i]][sample(1:2, size = 1), mut_sites[locus]]
      # ...ALL mutations would be added to the SAME homolog (1 or 2)
      for(locus in which(mut_sites == 1)){
        pop[[i]][sample(1:2, size = 1), locus] <-
          sample(dfe, size = 1)
      }
    }
  }
  
  # Return the population 
  return(pop)
}

# Modified version of MoL() where no fusions are possible
MiracleOfLife <- function(gametes){
  # New generation should be same size as number of couples
  newgen <- vector(mode = "list", length = length(gametes$eggs))
  
  # Scramble the order of the eggs and sperm to pick a mom and dad at random
  e <- sample(1:length(gametes$eggs))
  s <- sample(1:length(gametes$sperm))
  
  for(i in 1:length(newgen)){
    # Combine 1 gamete from each randomly selected parent
    egg_genome <- strsplit(gametes$eggs[e[i]], ",")[[1]]
    sperm_genome <- strsplit(gametes$sperm[s[i]], ",")[[1]]
    newgen[[i]] <- matrix(as.numeric(c(egg_genome, sperm_genome)), 
                          nrow = 2, byrow = T)
    
  }
  
  return(newgen)
}

# Run population for several rounds of mutation with no fusions, then use 
# modified GetFit function to calculate fitnesses of large and small chromosomes
gen_no <- 100
s <- 0
pop_size <- 1000
chiasm <- T
mu <- 0.5

# Get a distribution of fitness effects
dfe <- GetDFE()

# Get a new population
pop <- GetPop(pop_size)

# For each generation in "gen_no"...
for(gen in 1:gen_no){
  print(paste(c("Generation: ", gen), collapse = ""))
      
  # Mutate the starting population with no fusions
  pop <- ActofGod(pop, dfe, mu)
      
  # Assess the fitness of the mutated population
  if(gen == gen_no){
    # After all generations have been run, store the fitness of large and small chr
    large_small <- GetFit_Sizes(pop)
    break
  }else{
    fits <- GetFit(pop, s)
  }
  # Select parents based on the fitness
  parents <- Maury(pop, fits, length(pop))
  # Select gametes from these parents
  gametes <- MakeGametes(pop, parents, chiasm = T)
  # Breed the parents
  pop <- MiracleOfLife(gametes)
  
}

par(mfrow = c(1, 2), cex = 0.75)
hist(large_small[1,], main = "Fitness Values of Large Chromosomes")
abline(v = mean(large_small[1,]), col = "red")
hist(large_small[2,], main = "Fitness Values of Small Chromosomes")
abline(v = mean(large_small[2,]), col = "red")

for(i in 1:length(pop)){
  print(which(colSums(pop[[i]][]) != 2))
}
