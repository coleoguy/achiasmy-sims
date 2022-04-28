############################### 04-20 ##########################################
# Issue: Preliminary results of the achiasmy simulation show that fixations
# failed to occur under conditions of LOW mutation rate and HIGH s

# Hypothesis: As written, benefit of being SAME SEX homozygote is SO HIGH that
# NO heterozygotes at SAL are surviving such that ALL individuals are 
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



############################# Fitness at Fusion ################################
# Hypothesis 1.0: Individuals with fusions start out being more unfit than 
# individuals who do not get fusions
setwd("C:/Users/knigh/Documents/GitHub/achiasmy-sims/scripts")
source("sizes.R")

# Parameters:
gen_no <- 1000
pop_size <- 100
mu <- 0.00000001
s <- 1
fus.large <- F
chiasm <- T
fus.type <- "Y"
# Get probabilities of 0, 1, 2, and 3 mutations in a generation
mu.table <- as.data.frame(
  table(rbinom(n=100000, size=4*1587000, prob = mu))
)[1:4,2]/sum(as.data.frame(table(rbinom(n=100000, size=4*1587000, prob = mu)))[1:4,2])

# Initialize a dataframe with the following columns:
# Generation
# Fusion Status (Newly Fused, Not Newly Fused)
# Fitness
H1.0_results <- data.frame(matrix(ncol = 4, nrow = 0))
names(H1.0_results) <- c("Generation","Fusion Status", "Fitness",
                         "SAL Genotype of Fused Small Autosome")

# Run through ActofGod, but each generation where a fusion is introduced, store
# the generation number and the fusion status of all individuals in that 
# generation
ActofGod <- function(pop, dfe, fus.type, mu.table, fus.large, cur.gen){
  # Pick the individuals who get mutations, as well as the number of mutations
  # each mutated individual gets.
  hit <- sample(0:3, size = length(pop), replace = T, prob = mu.table)
  
  # Get indices of males and females to later sample from
  F.indx <- c()
  M.indx <- c()
  for(i in 1:length(pop)){
    if(pop[[i]][2, 13]){
      M.indx <- append(M.indx, i)
    }else{
      F.indx <- append(F.indx, i)
    }
  }
  # Each generation has 10% chance of having 1 individual with a fusion
  fuse.bool <- sample(0:1, size = 1, prob = c(0.9,0.1))
  
  # If this is one of the generations with a fusion, randomly sample an 
  # applicable chromosome (X or Y) to which the fusion should be introduced
  if(fuse.bool){
    # Add generation to temp vector to later add fitness of all fused and
    # unfused individuals
    temp <- data.frame(Generation = rep(cur.gen, length(pop)),
                       `Fusion Status` = rep(NA, length(pop)),
                       Fitness = rep(NA, length(pop)),
                       `SAL Genotype of Fused Small Autosome` = rep(NA, length(pop)))
    names(temp) <- c("Generation","Fusion Status", "Fitness",
                     "SAL Genotype of Fused Small Autosome")
    
    if(fus.type == "X"){
      # Choose whether a male (1) or female (0) gets the fusion. Females are
      # 2x as likely to get X fusions because they have 2 X chr
      fuse.sex <- sample(0:1, size = 1, prob = c(2, 1))
      # Sample a random individual of the select sex
      if(fuse.sex){
        fuse.indv <- sample(M.indx, size = 1)
        # If male was selected to be fused, sample only the X chr
        fuse.chr <- 1
      }else{
        fuse.indv <- sample(F.indx, size = 1)
        # If female was selected to be fused, sample 1 of the X chr
        fuse.chr <- sample(1:2, size = 1)
      }
      
      # Check if the sampled X is already fused
      if(!pop[[fuse.indv]][fuse.chr, 25]){
        temp$`Fusion Status`[fuse.indv] <- "Newly Fused"
        temp$`Fusion Status`[-fuse.indv] <- "Not Newly Fused"
        # If not, introduce large or small fusion based on value of "fus.large"
        if(fus.large){
          pop[[fuse.indv]][fuse.chr, 25] <- 2
        }else{
          pop[[fuse.indv]][fuse.chr, 25] <- 1
        }
      }else{
        temp$`Fusion Status`[fuse.indv] <- "Not Newly Fused"
      }
      # If only Y is permitted to have fusions, select a male, check if he already
      # has a fusion, then introduce fusion to Y
    }else if(fus.type == "Y"){
      fuse.indv <- sample(M.indx, size = 1)
      if(!pop[[fuse.indv]][2, 25]){
        temp$`Fusion Status`[fuse.indv] <- "Newly Fused"
        temp$`Fusion Status`[-fuse.indv] <- "Not Newly Fused"
        # introduce large or small fusion based on value of "fus.large"
        if(fus.large){
          pop[[fuse.indv]][2, 25] <- 2
        }else{
          pop[[fuse.indv]][2, 25] <- 1
        }
      # If he does already have a fusion, mark him as such
      }else{
        temp$`Fusion Status`[fuse.indv] <- "Not Newly Fused"
      }
    }
  }
  
  # Iterate through each individual
  for(i in 1:length(hit)){
    # If this is a mutated individual, sample the indicated number of sites to
    # be mutated from the non-SDL, non-SAL
    if(hit[i]){
      mut_sites <- sample((1:100)[-c(13,25,28,55)], size = hit[i], replace = F)
      # At each site to be mutated, replace the value of a random homolog with
      # a randomly-selected value from dfe
      pop[[i]][sample(1:2, size = 1, replace = F), mut_sites] <-
        sample(dfe, size = hit[i])
    }
  }
  if(fuse.bool){
    return(list(pop, temp))
  }else{
    return(list(pop, data.frame()))
  }
}

#### Run Evolve for 1 simulation on 1000 generations
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

# For each generation in "gen_no"...
for(gen in 1:gen_no){
  print(paste(c("Generation: ", gen), collapse = ""))
  # Mutate the starting population
  AoGres <- ActofGod(pop, dfe, fus.type, mu.table, fus.large, 
                     cur.gen = gen)
  pop <- AoGres[[1]]
  
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
  # If a fusion occurred in this generation, record the fitnesses and the SAL28
  # genotype of the small autosome which segregated with the Y
  if(nrow(AoGres[[2]])){
    H1.0_results <- rbind(H1.0_results, AoGres[[2]])
    H1.0_results$Fitness[H1.0_results$Generation == gen] <- fits
    for(i in 1:length(gametes$sperm)){
      cur.sperm <- as.numeric(strsplit(gametes$sperm[i], ",")[[1]])
      if((cur.sperm[25] == 1) & (cur.sperm[13] == 1)){
        # Collect genotype of fused chr at SAL28
        if(cur.sperm[28] == 1){
          fusedSAL28 <- "Male-Benefit"
        }else{
          fusedSAL28 <- "Female-Benefit"
        }
        H1.0_results$`SAL Genotype of Fused Small Autosome`[H1.0_results$Generation == gen & H1.0_results$`Fusion Status` == "Newly Fused"] <- fusedSAL28
      }
    } 
  }

  # Breed the parents
  pop <- MiracleOfLife(gametes)
  
  # Move to next two generation's rows
  gen_rows <- gen_rows + 4
  
  # Continue to next round of mutation
  next
}

final.pop <- pop

# Create a stripchart where...
# X = Generation
# Y = Fitness
# Color = Newly fused or not newly fused
library(ggplot2)
# Basic stripchart
ggplot(H1.0_results, 
       aes(x=as.factor(Generation), y=Fitness, color = `Fusion Status`)) + 
  geom_jitter(position=position_jitter(0.2), alpha = 0.5) + 
  theme_bw()

# Issue: Appears fused individuals have lower fitness
FusedFitH1.0 <- H1.0_results[H1.0_results$`Fusion Status` == "Newly Fused",]




