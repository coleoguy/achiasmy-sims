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
    pop[[i]][1,13] <- 1
    pop[[i]][2,13] <- sample(0:1, size = 1, replace = F, prob = rep(0.5,2))
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
GetFit <- function(pop, s){
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
        locus_fits <- append(locus_fits, 1 - s)
      # If SAL36 genotype is 01, append 1 - 0.5*s36 to vector
      }else if(sum(pop[[i]][,36]) == 1){
        locus_fits <- append(locus_fits, 1 - 0.5*s)
      }
      # If SAL36 genotype is 11, do nothing (appending 1 would not change product)
      
      # Repeat for SAL61
      if(sum(pop[[i]][,61]) == 0){
        locus_fits <- append(locus_fits, 1 - s)
      }else if(sum(pop[[i]][,61]) == 1){
        locus_fits <- append(locus_fits, 1 - 0.5*s)
      }
    # If the individual is a female...
    }else{
      # If SAL36 genotype is 11, append 1-s36 to vector
      if(sum(pop[[i]][,36]) == 2){
        locus_fits <- append(locus_fits, 1 - s)
        # If SAL36 genotype is 01, append 1 - 0.5*s36 to vector
      }else if(sum(pop[[i]][,36]) == 1){
        locus_fits <- append(locus_fits, 1 - 0.5*s)
      }
      # If SAL36 genotype is 00, do nothing (appending 1 would not change product)
      
      # Repeat for SAL61
      if(sum(pop[[i]][,61]) == 2){
        locus_fits <- append(locus_fits, 1 - s)
      }else if(sum(pop[[i]][,61]) == 1){
        locus_fits <- append(locus_fits, 1 - 0.5*s)
      }
      
    }
    # Take product of fitness vector and add to vector of population fitness
    # values
    pop_fits[i] <- prod(locus_fits)
  }
  # Return vector of fitness values of each individual
  return(pop_fits)
}

# TODO Do we want to sample from a new gamma distribution each time a mutation
# occurs, or store a constant gamma distribution (as I have here) and keep
# drawing new numbers from it?
# 5 get a store distribution of fitness effects
GetDFE <- function(){
  dfe <- rgamma(5000, shape = .28, scale=113)
  # values less than -1 are not meaningful
  dfe[dfe>1] <- 1
  dfe <- 1-dfe
  return(dfe)  
}

# 6 lay down mutations
# mut_rate = mutation rate
# pop = list of N genomes
ActofGod <- function(pop, dfe){
  # Pick the individuals who get mutations, as well as the number of mutations
  # each mutated individual gets.
  # Probabilities of each number of mutations were derived from the following code:
  # as.data.frame(table(rbinom(n=100000, size=4*1587000,prob = 1.45*10^-8)))/100000
  hit <- sample(0:3, size= length(pop), 
                prob=c(0.91148, 0.08430, 0.00411, 0.00011),replace=T)
  
  # Iterate through each individual
  for(i in 1:length(pop)){
    # If this is a mutated individual, sample the indicated number of sites to
    # be mutated from the non-SDL, non-SAL
    if(hit[i]){
      mut_sites <- sample((1:100)[-c(13,36,61)], size = hit[i], replace = F)
    # At each site to be mutated, replace the value of a random homolog with
    # a randomly-selected value from dfe
      pop[[i]][sample(1:2, size = 1, replace = F), mut_sites] <-
        sample(dfe, size = hit[i])
    }
  }
  return(pop)
}

# 3 pick parents
GetParents <- function(pop, fits, N = length(pop)){
  sexes <- rep("fem",length(pop))
  for(i in 1:N){
    if(pop[[i]][2,13]){
      sexes[i] <- "mal"
    }
  }
  moms <- sample((1:N)[sexes=="fem"], prob=fits[sexes=="fem"], 
                 size=length(pop)/2, replace = T)
  dads <- sample((1:N)[sexes=="mal"], prob=fits[sexes=="mal"], 
                 size=length(pop)/2, replace = T)
  parents <- list(moms,dads)
  names(parents) <- c("Moms","Dads")
  return(parents)
}
# 4 make gametes
# pop = list of genomes
# parent = list of two sets of integers describing the individuals who get to 
# breed
# PARb = Starting locus of pseudoautosomal region of the male sex chromosomes
# Output (gametes) = list of two sets of vectors each describing the sequence
# of a haploid gamete
MakeGametes <- function(pop, parents, PARb){
  gametes <- vector(mode="list",length=2)
  names(gametes) <- c("eggs","sperm")
  
  # Iterate through each pair of parents
  for(i in 1:length(parents$Dads)){
  # For male...
    # Check whether individual is chiasmatic or not
    # TODO: Change this portion of code so individual is NOT always set as 
    # chiasmatic
    chiasm <- T
    # If individual is chiasmatic...
    if(chiasm){
      # Pick 3 random sites OTHER than 13, 36, and 61 in...
        # PAR ((PARb+1):24)
        # Autosome I (27:49)
        # Autosome II (52-99)
      # ... at which recombination occurs (recombination CANNOT occur at first or
      # final locus of a chromosome, so these loci are excluded)
      SexRec <- sample((PARb+1):24, 1)
      Chr1Rec <- sample(27:49, 1)
      Chr2Rec <- sample(52:99, 1)
      
      # For each chromosome, create two recombinant homologs and add to vectors
      # from which gametes will be selected
      SexChrGametes <- c(
        paste(c(pop[[parents$Dads[i]]][2,1:(SexRec-1)], 
                pop[[parents$Dads[i]]][1,SexRec:25]), collapse = ""),
        paste(c(pop[[parents$Dads[i]]][1,1:(SexRec-1)], 
                pop[[parents$Dads[i]]][2,SexRec:25]), collapse = "")
      )
      
      Chr1Gametes <- c(
        paste(c(pop[[parents$Dads[i]]][2,26:(Chr1Rec-1)], 
                pop[[parents$Dads[i]]][1,Chr1Rec:50]), collapse = ""),
        paste(c(pop[[parents$Dads[i]]][1,26:(Chr1Rec-1)], 
                pop[[parents$Dads[i]]][2,Chr1Rec:50]), collapse = "")
      )
      
      Chr2Gametes <- c(
        paste(c(pop[[parents$Dads[i]]][2,51:(Chr2Rec-1)], 
                pop[[parents$Dads[i]]][1,Chr2Rec:100]), collapse = ""),
        paste(c(pop[[parents$Dads[i]]][1,51:(Chr2Rec-1)], 
                pop[[parents$Dads[i]]][2,Chr2Rec:100]), collapse = "")
      )
      
      # For each chromosome, randomly select ONE homolog from the two 
      # recombinant homologs and stick the 3 chromosomes together to create a
      # haploid genome, adding this genome to the list of gametes under "sperm"
      gametes$sperm[i] <- paste(c(sample(SexChrGametes, 1),
                                  sample(Chr1Gametes, 1),
                                  sample(Chr2Gametes, 1)), collapse = "")
    }
  # For female..
    # Pick 3 random sites OTHER than 13, 36, and 61 in...
      # Sex chromosome (2:24)
      # Autosome I (27:49)
      # Autosome II (52-99)
    # ... at which recombination occurs (recombination CANNOT occur at first or
    # final locus of a chromosome, so these loci are excluded)
    SexRec <- sample(2:24, 1)
    Chr1Rec <- sample(27:49, 1)
    Chr2Rec <- sample(52:99, 1)
    
    # For each chromosome, create two recombinant homologs and add to vectors
    # from which gametes will be selected
    SexChrGametes <- c(
      paste(c(pop[[parents$Moms[i]]][2,1:(SexRec-1)], 
              pop[[parents$Moms[i]]][1,SexRec:25]), collapse = ""),
      paste(c(pop[[parents$Moms[i]]][1,1:(SexRec-1)], 
              pop[[parents$Moms[i]]][2,SexRec:25]), collapse = "")
    )
    
    Chr1Gametes <- c(
      paste(c(pop[[parents$Moms[i]]][2,26:(Chr1Rec-1)], 
              pop[[parents$Moms[i]]][1,Chr1Rec:50]), collapse = ""),
      paste(c(pop[[parents$Moms[i]]][1,26:(Chr1Rec-1)], 
              pop[[parents$Moms[i]]][2,Chr1Rec:50]), collapse = "")
    )
    
    Chr2Gametes <- c(
      paste(c(pop[[parents$Moms[i]]][2,51:(Chr2Rec-1)], 
              pop[[parents$Moms[i]]][1,Chr2Rec:100]), collapse = ""),
      paste(c(pop[[parents$Moms[i]]][1,51:(Chr2Rec-1)], 
              pop[[parents$Moms[i]]][2,Chr2Rec:100]), collapse = "")
    )
    # For each chromosome, randomly select one homolog from the two 
    # recombinant homologs and stick these chromosomes together to create a
    # haploid genome, adding this genome to the list of gametes under "eggs"
    gametes$eggs[i] <- paste(c(sample(SexChrGametes, 1),
                                sample(Chr1Gametes, 1),
                                sample(Chr2Gametes, 1)), collapse = "")
  }
  
  # Return list of eggs and sperm
  return(gametes)
}

# 7 make next gen
# gametes = List of two vectors describing the haplotypes of available eggs and
# sperm

# 8 return to step 2


