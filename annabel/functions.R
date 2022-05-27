# The purpose of this code is to compare the impact of
# achiasmatic meiosis on fusions of large or small autosomes
# with sex chromosomes we will build a canonically inspired
# model where the sex chromosome has an SDL and the autosomes have 
# SA loci

#  sex chromosome         autosome 1         autosome 2
#  1-12 SDL (PAR)14-25    26-27 SA 29-50     51-54 SA 56-100

# SDL 0=X 1=Y
# SA 0=female benefit; 1=male benefit
# All other loci are general fitness loci with values from 0-1
# reflecting their fitness. Fitness will be multiplicative

# Fusions are indicated at locus 25 by having a:
# 0 = no fusion
# 1 = small fusion
# 2 = large fusion

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
    pop[[i]][1,13] <- 0
    pop[[i]][2,13] <- sample(0:1, size = 1, replace = F, prob = rep(0.5,2))
    # Start each individual with NO sex chr-autosome fusions
    pop[[i]][,25] <- 0
    
    pop[[i]][,28] <- sample(0:1, size = 2, replace = T, prob = rep(0.5,2))
    pop[[i]][,55] <- sample(0:1, size = 2, replace = T, prob = rep(0.5,2))
  }
  return(pop)
}

# 2 assess fitness
# pop = list of genomes from GetPop()
# s = selection coefficient on sexually-antagonistic loci
GetFit <- function(pop, s){
  # Initialize a vector of fitness values for all individuals in the population
  pop_fits <- vector(mode = "numeric", length = length(pop))
  
  # Iterate through each individual in the population
  for(i in 1:length(pop)){
    # Create a vector of averaged fitness values of ALL homologous loci, 
    # excluding the SDL and SAL
    locus_fits <- colSums(pop[[i]][,-c(13,25,28,55)])/2
    
    # If the individual is a male...
    if(1 %in% pop[[i]][,13]){
      # If SAL28 genotype is 00, append 1-s to vector
      if(sum(pop[[i]][,28]) == 0){
        locus_fits <- append(locus_fits, 1 - s)
        # If SAL28 genotype is 01, append 1 - 0.5*s to vector
      }else if(sum(pop[[i]][,28]) == 1){
        locus_fits <- append(locus_fits, 1 - 0.5*s)
      }
      # If SAL28 genotype is 11, do nothing (appending 1 would not change product)
      
      # Repeat for SAL55
      if(sum(pop[[i]][,55]) == 0){
        locus_fits <- append(locus_fits, 1 - s)
      }else if(sum(pop[[i]][,55]) == 1){
        locus_fits <- append(locus_fits, 1 - 0.5*s)
      }
      # If the individual is a female...
    }else{
      # If SAL28 genotype is 11, append 1-s to vector
      if(sum(pop[[i]][,28]) == 2){
        locus_fits <- append(locus_fits, 1 - s)
        # If SAL28 genotype is 01, append 1 - 0.5*s to vector
      }else if(sum(pop[[i]][,28]) == 1){
        locus_fits <- append(locus_fits, 1 - 0.5*s)
      }
      # If SAL28 genotype is 00, do nothing (appending 1 would not change product)
      
      # Repeat for SAL55
      if(sum(pop[[i]][,55]) == 2){
        locus_fits <- append(locus_fits, 1 - s)
      }else if(sum(pop[[i]][,55]) == 1){
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
# fus.large: 
#         T = fusions can occur only between sex and LARGE autosomes
#         F = fusions can occur only between sex and SMALL autosomes
ActofGod <- function(pop, dfe, fus.type, mu, fus.large){
  mu.table <- as.data.frame(
    table(rbinom(n=100000, size=4*1587000, prob = mu))
  )[1:4,2]/sum(as.data.frame(table(rbinom(n=100000, size=4*1587000, prob = mu)))[1:4,2])
  
  # Pick the individuals who get mutations, as well as the number of mutations
  # each mutated individual gets.
  hit <- sample(0:3, size = length(pop), replace = T, prob = mu.table)
  # Iterate through each individual
  for(i in which(hit != 0)){
    # If this is a mutated individual, sample the indicated number of sites to
    # be mutated from the non-SDL, non-SAL
      mut_sites <- sample((1:100)[-c(13,25,28,55)], size = hit[i], replace = F)
      # At each site to be mutated, replace the value of a random homolog with
      # a randomly-selected value from dfe
      pop[[i]][sample(1:2, size = 1, replace = F), mut_sites] <-
        sample(dfe, size = hit[i])
  }
  # Get indices of males and females to later sample from
  F.indx <- c()
  M.indx <- c()
  
  
  
  checker <- function(x){
    x[2,13] == 0
  }
  
  females <- unlist(lapply(pop, FUN=checker))
  
  # for(i in 1:length(pop)){
  #   if(pop[[i]][2, 13]){
  #     M.indx <- append(M.indx, i)
  #   }else{
  #     F.indx <- append(F.indx, i)
  #   }
  # }
  # Each generation has 10% chance of having 1 individual with a fusion
  fuse.bool <- sample(0:1, size = 1, prob = c(0.9,0.1))
  fuse.bool <- F
  # If this is one of the generations with a fusion, randomly sample an 
  # applicable chromosome (X or Y) to which the fusion should be introduced
  if(fuse.bool){
    if(fus.type == "X"){
      # Choose whether a male (1) or female (0) gets the fusion. Females are
      # 2x as likely to get X fusions because they have 2 X chr
      fuse.sex <- sample(0:1, size = 1, prob = c(2, 1))
      # Sample a random individual of the select sex
      if(fuse.sex){
        fuse.indv <- sample(which(!females), size = 1)
        # If male was selected to be fused, sample only the X chr
        fuse.chr <- 1
      }else{
        fuse.indv <- sample(which(females), size = 1)
        # If female was selected to be fused, sample 1 of the X chr
        fuse.chr <- sample(1:2, size = 1)
      }
        
      # Check if the sampled X is already fused
      if(!pop[[fuse.indv]][fuse.chr, 25]){
        # If not, introduce large or small fusion based on value of "fus.large"
        if(fus.large){
          pop[[fuse.indv]][fuse.chr, 25] <- 2
        }else{
          pop[[fuse.indv]][fuse.chr, 25] <- 1
        }
      }
    # If only Y is permitted to have fusions, select a male, check if he already
    # has a fusion, then introduce fusion to Y
    }else if(fus.type == "Y"){
      fuse.indv <- sample(which(!females), size = 1)
      if(!pop[[fuse.indv]][2, 25]){
        # introduce large or small fusion based on value of "fus.large"
        if(fus.large){
          pop[[fuse.indv]][2, 25] <- 2
        }else{
          pop[[fuse.indv]][2, 25] <- 1
        }
      }
    }
  }
  
  return(pop)
}

# 3 find parents
# N = Total number of COUPLES/individuals to contribute to next generation
Maury <- function(pop, fits, N = length(pop)){
  sexes <- rep("fem",N)
  for(i in 1:length(sexes)){
    if(pop[[i]][2,13]){
      sexes[i] <- "mal"
    }
  }
  # Checks to ensure that there is more than 1 female before sampling.
  # If you don't have this check, you will get an "incorrect number of 
  # probabilities" error in the event of a single female (same for males)
  if(sum(sexes=="fem") > 1){
    moms <- sample((1:N)[sexes=="fem"], prob=fits[sexes=="fem"], 
                   size=N, replace = T)
  }else{
    moms <- rep((1:N)[sexes=="fem"], N)
  }
  
  if(sum(sexes=="mal") > 1){
    dads <- sample((1:N)[sexes=="mal"], prob=fits[sexes=="mal"], 
                   size=N, replace = T)
  }else{
    dads <- rep((1:N)[sexes=="mal"], N)
  }
  
  parents <- list(moms,dads)
  names(parents) <- c("Moms","Dads")
  return(parents)
}
# 4 make gametes
# pop = list of genomes
# parent = list of two sets of integers describing the individuals in pop who 
# get to breed
# Output (gametes) = list of two sets of vectors each describing the sequence
# of a haploid gamete
MakeGametes <- function(pop, parents, chiasm = T){
  gametes <- vector(mode="list",length=2)
  names(gametes) <- c("eggs","sperm")
  
  # Iterate through each pair of parents
  for(i in 1:length(parents$Dads)){
    # For male...
    # If individual is chiasmatic...
    if(chiasm){
      # Pick 3 random sites OTHER than 13, 28, and 55 in...
      # PAR (2:13), no recombination from 14:25
      # Autosome I (27:49)
      # Autosome II (52-99)
      # ... at which recombination occurs (recombination CANNOT occur at first or
      # final locus of a chromosome, so these loci are excluded)
      SexRec <- sample(2:13, 1)
      Chr1Rec <- sample(27:49, 1)
      Chr2Rec <- sample(52:99, 1)
      
      # For each chromosome, create two recombinant homologs and add to vectors
      # from which gametes will be selected
      SexChrGametes <- c(
        paste(c(pop[[parents$Dads[i]]][2,1:(SexRec-1)], 
                pop[[parents$Dads[i]]][1,SexRec:25]), collapse = ","),
        paste(c(pop[[parents$Dads[i]]][1,1:(SexRec-1)], 
                pop[[parents$Dads[i]]][2,SexRec:25]), collapse = ",")
      )
      
      Chr1Gametes <- c(
        paste(c(pop[[parents$Dads[i]]][2,26:(Chr1Rec-1)], 
                pop[[parents$Dads[i]]][1,Chr1Rec:50]), collapse = ","),
        paste(c(pop[[parents$Dads[i]]][1,26:(Chr1Rec-1)], 
                pop[[parents$Dads[i]]][2,Chr1Rec:50]), collapse = ",")
      )
      
      Chr2Gametes <- c(
        paste(c(pop[[parents$Dads[i]]][2,51:(Chr2Rec-1)], 
                pop[[parents$Dads[i]]][1,Chr2Rec:100]), collapse = ","),
        paste(c(pop[[parents$Dads[i]]][1,51:(Chr2Rec-1)], 
                pop[[parents$Dads[i]]][2,Chr2Rec:100]), collapse = ",")
      )
      
      # For each chromosome, randomly select ONE homolog from the two 
      # recombinant homologs and stick the 3 chromosomes together to create a
      # haploid genome, adding this genome to the list of gametes under "sperm"
      if(sum(pop[[parents$Dads[i]]][,25]) == 0){
        gametes$sperm[i] <- paste(c(sample(SexChrGametes, 1),
                                    sample(Chr1Gametes, 1),
                                    sample(Chr2Gametes, 1)), collapse = ",")
      }else if(1 %in% pop[[parents$Dads[i]]][,25]){
        pick <- sample(1:2, 1)
        gametes$sperm[i] <- paste(c(SexChrGametes[pick],
                                    Chr1Gametes[pick],
                                    sample(Chr2Gametes, 1)), collapse = ",")
      }else if(2 %in% pop[[parents$Dads[i]]][,25]){
        pick <- sample(1:2, 1)
        gametes$sperm[i] <- paste(c(SexChrGametes[pick],
                                    sample(Chr1Gametes, 1),
                                    Chr2Gametes[pick]), collapse = ",")
      }
      
    }else{
      pick <- sample(1:2, 3, replace = T)
      gametes$sperm[i] <- paste(c(pop[[parents$Dads[i]]][pick[1],1:25],
                                  pop[[parents$Dads[i]]][pick[2],26:50],
                                  pop[[parents$Dads[i]]][pick[3],51:100]), 
                                collapse = ",")
    }
    # For female..
    # Pick 3 random sites OTHER than 13, 28, and 55 in...
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
              pop[[parents$Moms[i]]][1,SexRec:25]), collapse = ","),
      paste(c(pop[[parents$Moms[i]]][1,1:(SexRec-1)], 
              pop[[parents$Moms[i]]][2,SexRec:25]), collapse = ",")
    )
    
    Chr1Gametes <- c(
      paste(c(pop[[parents$Moms[i]]][2,26:(Chr1Rec-1)], 
              pop[[parents$Moms[i]]][1,Chr1Rec:50]), collapse = ","),
      paste(c(pop[[parents$Moms[i]]][1,26:(Chr1Rec-1)], 
              pop[[parents$Moms[i]]][2,Chr1Rec:50]), collapse = ",")
    )
    
    Chr2Gametes <- c(
      paste(c(pop[[parents$Moms[i]]][2,51:(Chr2Rec-1)], 
              pop[[parents$Moms[i]]][1,Chr2Rec:100]), collapse = ","),
      paste(c(pop[[parents$Moms[i]]][1,51:(Chr2Rec-1)], 
              pop[[parents$Moms[i]]][2,Chr2Rec:100]), collapse = ",")
    )
    # For each chromosome, randomly select one homolog from the two 
    # recombinant homologs and stick these chromosomes together to create a
    # haploid genome, adding this genome to the list of gametes under "eggs"
    if(sum(pop[[parents$Moms[i]]][,25]) == 0){
      gametes$eggs[i] <- paste(c(sample(SexChrGametes, 1),
                                 sample(Chr1Gametes, 1),
                                 sample(Chr2Gametes, 1)), collapse = ",")
    }else if(1 %in% pop[[parents$Moms[i]]][,25]){
      pick <- sample(1:2, 1)
      gametes$eggs[i] <- paste(c(SexChrGametes[pick],
                                 Chr1Gametes[pick],
                                 sample(Chr2Gametes, 1)), collapse = ",")
    }else if(2 %in% pop[[parents$Moms[i]]][,25]){
      pick <- sample(1:2, 1)
      gametes$eggs[i] <- paste(c(SexChrGametes[pick],
                                 sample(Chr1Gametes, 1),
                                 Chr2Gametes[pick]), collapse = ",")
    }
  }
  
  # Return list of eggs and sperm
  return(gametes)
}
# 7 make next gen
# gametes = List of strings describing the haplotypes of available eggs and
# sperm where each locus' fitness value is separated by a comma
# newgen = List of matrices describing each of the new individuals
MiracleOfLife <- function(gametes){
  # New generation should be same size as number of couples
  newgen <- vector(mode = "list", length = length(gametes$eggs))
  
  # Until list of gametes is empty, sample a new mom and dad gamete and fuse
  # together into a 2X100 matrix, adding to population
  for(i in 1:length(newgen)){
    # pick a mom and dad at random from the mating pool
    e <- sample(1:length(gametes$eggs), 1)
    s <- sample(1:length(gametes$sperm), 1)
    
    # sex
    newgen[[i]] <- matrix(as.numeric(c(strsplit(gametes$eggs[e], ",")[[1]],
                                       strsplit(gametes$sperm[s], ",")[[1]])),
                          nrow = 2, byrow = T)
    gametes$eggs <- gametes$eggs[-e]
    gametes$sperm <- gametes$sperm[-s]
  }
  
  return(newgen)
}

# 8 return to step 2

# Runs the simulation on "pop_size" number of individuals for "gen_no" 
# generations with an SAL selection coefficient of "s" 
# fus.type "X" or "Y" for which chrom is fused to.
# fus.large: 
#         T = fusions can occur only between sex and LARGE autosomes
#         F = fusions can occur only between sex and SMALL autosomes
Evolve <- function(pop_size, gen_no, s, chiasm, fus.type, mu.table, fus.large){
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
  # Output final genome after this many generations
  return(results)
}



genotyper <- function(pop, locus){
  one.count <- 0
  for(i in 1:length(pop)){
    one.count <- one.count + sum(pop[[i]][,locus])
  }
  one.count/(2*length(pop))
}



