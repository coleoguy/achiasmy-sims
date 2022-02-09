# The purpose of this code is to compare the impact of
# achiasmatic meiosis on fusions of large or small autosomes
# with sex chromosomes we will build a canonically inspired
# model where the sex chromosome has an SDL and the autosomes have 
# SA loci

#  sex chromosome         autosome 1         autosome 2
#  1-12(PAR) SDL 14-25    26-35 SA 37-50     51-60 SA 62-100
#

# SDL 0=X 1=Y
# SA 0=female benefit; 1=male benefit
# All other loci are general fitness loci with values from 0-1
# reflecting their fitness. Fitness will be multiplicative

# genotype at SA    male      female
#  00                1-s      1
#  01                1-.5s    1-.5s
#  11                1        1-s
#

# for all models females will recombine on all chromosomes
# for achiasmatic model males will have no recombination on
# any chromosomes but females will recombine as normal
# males that are chiasmatic will recombine only in the PAR 
# region on the sex chromosome.



# variables
# s: the selection coefficient for SA .5
# dfe: vector of selection coefficients for general fitness loci
# gen: number of generations to run
# iter: number of trials to use
# 

# 1 make random starting genomes
# Input the number of genomes, output a list where each entry is a genome
GetPop <- function(N){
  pop <- vector(mode = "list", length = N)
  for(i in 1:N){
    # Set all general fitness loci equal to 1
    pop[[i]] <- matrix(rep(1, 200), nrow = 2, ncol = 100)
    # Set each SDL to 1 or 0 with prob of 0.5, because 50% chance of being male
    # or female, but make it impossible to have 2 Y chromosomes
    pop[[i]][,13] <- sample(c(0,0,0,1), size = 2, replace = F)
    # TODO Should SAL have equal probability of being male or female benefit?
    pop[[i]][,36] <- sample(0:1, size = 2, replace = T, prob = rep(0.5,2))
    pop[[i]][,61] <- sample(0:1, size = 2, replace = T, prob = rep(0.5,2))
  }
  return(pop)
}

# 2 assess fitness
# pop = list of genomes from GetPop()
# s36 = selection coefficient on sexually-antagonistic locus at 36
# s61 = selection coefficient on sexually-antagonistic locus at 61
GetFit <- function(pop, s36, s61){
  # Initialize a vector of fitness values for all individuals in the population
  pop_fits <- vector(mode = "numeric", length = length(pop))
  
  # Iterate through each individual in the population
  for(i in 1:length(pop)){
    # Create a vector of averaged fitness values of ALL homologous loci, 
    # excluding the SDL and SAL
    locus_fits <- colSums(pop[[i]][,-c(13,36,61)])/2
    
    # If the individual is a male...
    if(1 %in% pop[[i]][,13]){
      # If SAL36 genotype is 00, append 1-s36 to vector
      if(sum(pop[[i]][,36]) == 0){
        locus_fits <- append(locus_fits, 1 - s36)
      # If SAL36 genotype is 01, append 1 - 0.5*s36 to vector
      }else if(sum(pop[[i]][,36]) == 1){
        locus_fits <- append(locus_fits, 1 - 0.5*s36)
      }
      # If SAL36 genotype is 11, do nothing (appending 1 would not change product)
      
      # Repeat for SAL61
      if(sum(pop[[i]][,61]) == 0){
        locus_fits <- append(locus_fits, 1 - s61)
      }else if(sum(pop[[i]][,61]) == 1){
        locus_fits <- append(locus_fits, 1 - 0.5*s61)
      }
    # If the individual is a female...
    }else{
      # If SAL36 genotype is 11, append 1-s36 to vector
      if(sum(pop[[i]][,36]) == 2){
        locus_fits <- append(locus_fits, 1 - s36)
        # If SAL36 genotype is 01, append 1 - 0.5*s36 to vector
      }else if(sum(pop[[i]][,36]) == 1){
        locus_fits <- append(locus_fits, 1 - 0.5*s36)
      }
      # If SAL36 genotype is 00, do nothing (appending 1 would not change product)
      
      # Repeat for SAL61
      if(sum(pop[[i]][,61]) == 2){
        locus_fits <- append(locus_fits, 1 - s61)
      }else if(sum(pop[[i]][,61]) == 1){
        locus_fits <- append(locus_fits, 1 - 0.5*s61)
      }
      
    }
    # Take product of fitness vector and add to vector of population fitness
    # values
    pop_fits[[i]] <- prod(locus_fits)
  }
  # Return vector of fitness values of each individual
  return(pop_fits)
}

# 3 pick parents

# 4 make gametes

# 5 lay down mutations
dfe <- rgamma(5000, shape = .28, scale=113)
# values less than -1 are not meaningful
dfe[dfe>1] <- 1
dfe <- 1-dfe
hist(dfe)
# 6 make next gen

# 7 return to step 2


