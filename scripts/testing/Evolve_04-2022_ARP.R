# Issue: Preliminary results of the achiasmy simulation show that fixations
# failed to occur under conditions of LOW mutation rate and HIGH s

# Hypothesis: As written, benefit of being SAME SEX homozygote is SO HIGH that
# NOT heterozygotes at SAL are surviving such that ALL individuals are 
# homozygous for same sex SAL and fusions thus have no benefit.

setwd("C:/Users/knigh/Documents/GitHub/achiasmy-sims/scripts")
source("sizes.R")

# Chunk of code from Evolve():

# Get a new population
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

fus.large <- T
gen_no <- 1000
chiasm <- T
fus.type <- "Y"
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
  gametes <- MakeGametes(pop, parents, chiasm)
  
  # Breed the parents
  pop <- MiracleOfLife(gametes)
  
  # Move to next two generation's rows
  gen_rows <- gen_rows + 4
  
  # Continue to next round of mutation
  next
}

final.pop <- pop
# Get number of individuals with different fitness values at final generation
hist(fits, breaks = 50)

# Genotype everyone at final generation
# MHomo_28_# = Homozygous for male-benefit at 28
# MHomo_55_# = Homozygous for male-benefit at 55
# MHomo_B_# = Homozygous for male-benefit at BOTH
MHomo_28_M <- c()
MHomo_55_M <- c()
MHomo_B_M <- c()

FHomo_28_M <- c()
FHomo_55_M <- c()
FHomo_B_M <- c()

Het_28_M <- c()
Het_55_M <- c()
Het_B_M <- c()

MHomo_28_F <- c()
MHomo_55_F <- c()
MHomo_B_F <- c()

FHomo_28_F <- c()
FHomo_55_F <- c()
FHomo_B_F <- c()

Het_28_F <- c()
Het_55_F <- c()
Het_B_F <- c()
for(i in 1:length(pop)){
  if(sum(pop[[i]][,28] == 2)){
    if(pop[[i]][2,13]){
      MHomo_28_M <- append(MHomo_28_M, i)
    }else{
      MHomo_28_F <- append(MHomo_28_F, i)
    }
  }else if(sum(pop[[i]][,28] == 1)){
    if(pop[[i]][2,13]){
      Het_28_M <- append(Het_28_M, i)
    }else{
      Het_28_F <- append(Het_28_F, i)
    }
  }else{
    if(pop[[i]][2,13]){
      FHomo_28_M <- append(FHomo_28_M, i)
    }else{
      FHomo_28_F <- append(FHomo_28_F, i)
    }
  }
  if(sum(pop[[i]][,55] == 2)){
    if(pop[[i]][2,13]){
      MHomo_55_M <- append(MHomo_55_M, i)
    }else{
      MHomo_55_F <- append(MHomo_55_F, i)
    }
  }else if(sum(pop[[i]][,55] == 1)){
    if(pop[[i]][2,13]){
      Het_55_M <- append(Het_55_M, i)
    }else{
      Het_55_F <- append(Het_55_F, i)
    }
  }else{
    if(pop[[i]][2,13]){
      FHomo_55_M <- append(FHomo_55_M, i)
    }else{
      FHomo_55_F <- append(FHomo_55_F, i)
    }
  }
  if(sum(pop[[i]][,28] == 2) & sum(pop[[i]][,55] == 2)){
    if(pop[[i]][2,13]){
      MHomo_B_M <- append(MHomo_B_M, i)
    }else{
      MHomo_B_F <- append(MHomo_B_F, i)
    }
  }else if(sum(pop[[i]][,28] == 1) & sum(pop[[i]][,55] == 1)){
    if(pop[[i]][2,13]){
      Het_B_M <- append(Het_B_M, i)
    }else{
      Het_B_F <- append(Het_B_F, i)
    }
  }else{
    if(pop[[i]][2,13]){
      FHomo_B_M <- append(FHomo_B_M, i)
    }else{
      FHomo_B_F <- append(FHomo_B_F, i)
    }
  }
}
