############################ Annabel Perry June 2022 ###########################
#The purpose of this code is to identify why fusions are not fixing in any sims#
################################################################################

#  sex chromosome         autosome 1         autosome 2
#  1-12 SDL (PAR)14-25    26-27 SA 29-50     51-54 SA 56-100

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

##################################### Functions ################################
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

GetDFE <- function(){
  dfe <- rgamma(5000, shape = .28, scale=113)
  # values less than -1 are not meaningful
  dfe[dfe>1] <- 1
  dfe <- 1-dfe
  return(dfe)  
}

ActofGod <- function(pop, dfe, fus.type, mu, fus.large, FitnessTracker, 
                     PaternalGenealogy, gen){
  # Initialize a vector to store sites which receive mutations
  mut_sites <- rep(0, ncol(pop[[1]]))
  # Iterate through each individual
  for(i in 1:length(pop)){
    # Decide which sites in this individual, if any, get mutations
    mut_sites[-c(13,25,28,55)] <- sample(0:1, size = (length(mut_sites) - 4), prob = c(1-mu, mu), replace = T)
    print(paste(c("Sites ", paste(which(mut_sites == 1), collapse = " & "), " in individual ", i, " got mutations."), collapse = ""))
    # If the individual gets any mutations...
    if(sum(mut_sites)){
      # At each site to be mutated, replace the value of a random homolog with
      # a randomly-selected value from dfe
      
      # In case Heath asks why you are LOOPING through all mutated sites:
      # If you were to instead do this...
      # pop[[i]][sample(1:2, size = 1), mut_sites[locus]]
      # ...ALL mutations would be added to the SAME homolog (1 or 2)
      for(locus in which(mut_sites == 1)){
        pop[[i]][sample(1:2, size = 1), mut_sites[locus]] <-
          sample(dfe, size = 1)
      }
    }
  }
  # Get indices of males and females to later sample from
  F.indx <- c()
  M.indx <- c()
  
  checker <- function(x){
    x[2,13] == 0
  }
  
  females <- unlist(lapply(pop, FUN=checker))
  
  # Each generation has 10% chance of having 1 individual with a fusion
  fuse.bool <- sample(0:1, size = 1, prob = c(0.9,0.1))
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
        # Mark lineage of current individual as "fused"
        FitnessTracker$`Fused?`[(FitnessTracker$Generation == gen) &
                                  (FitnessTracker$`Paternal Lineage` == PaternalGenealogy[fuse.indv])] <- 1
        # introduce large or small fusion based on value of "fus.large"
        if(fus.large){
          pop[[fuse.indv]][2, 25] <- 2
        }else{
          pop[[fuse.indv]][2, 25] <- 1
        }
      }
    }
  }
  
  return(list(pop, FitnessTracker))
}

Maury <- function(pop, fits, N = length(pop), PaternalGenealogy){
  sexes <- rep("fem",N)
  for(i in 1:length(sexes)){
    if(pop[[i]][2,13]){
      sexes[i] <- "mal"
    }
  }
  
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
  # Record the lineages of all individuals who were selected as dads
  DadLineages <- PaternalGenealogy[dads]
  
  # Find all lineages which were selected multiple times
  Multi.Lineages <- DadLineages[duplicated(DadLineages)]
  
  # For each lineage which was selected multiple times
  for(d in 1:length(Multi.Lineages)){
    # If the lineage is ALREADY an offshoot...
    if(as.numeric(Multi.Lineages[d])%%1 != 0){
      # Find the number of times the lineage was selected
      Offshoots <- grep(Multi.Lineages[d], DadLineages)
      # For each number of times the lineage was selected, create a new lineage
      # and store the new lineage in DadLineages
      for(o in 1:length(Offshoots)){
        DadLineages[Offshoots[o]] <- as.numeric(paste(c(Multi.Lineages[d],
                                                        o), collapse = ""))
      }
      # If the lineage is not already an offshoot...
    }else{
      # Find the number of times the lineage was selected
      Offshoots <- grep(Multi.Lineages[d], DadLineages)
      # For each number of times the lineage was selected, create a new lineage
      # WITH A DECIMAL POINT and store the new lineage in DadLineages
      for(o in 1:length(Offshoots)){
        DadLineages[Offshoots[o]] <- as.numeric(paste(c(Multi.Lineages[d], ".",
                                                        o), collapse = ""))
      }
    }
  }
  
  parents <- list(moms,dads)
  names(parents) <- c("Moms","Dads")
  return(list(parents, DadLineages))
}

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

MiracleOfLife <- function(gametes, PaternalGenealogy, DadLineages, gen, FitnessTracker){
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
    # Check if the sperm is X or Y-bearing
    cur.sperm <- as.numeric(strsplit(gametes$sperm[s], ",")[[1]])
    
    # If sperm is Y-bearing... 
    if(cur.sperm[13]){
      # add the dad's lineage to corresponding element of PaternalLineages
      PaternalGenealogy[i] <- DadLineages[s]
      # create a new row in FitnessTracker with the dad's lineage AND
      # record SAL genotype in row corresponding to the current male
      new.row <- c((gen + 1), DadLineages[s], cur.sperm[25], cur.sperm[28], NA)
      FitnessTracker <- rbind(FitnessTracker, new.row)
      # If sperm is X-bearing... 
    }else{
      # add an NA to corresponding element of PaternalLineages
      PaternalGenealogy[i] <- NA
    }
    
    gametes$eggs <- gametes$eggs[-e]
    gametes$sperm <- gametes$sperm[-s]
  }
  
  return(list(newgen, PaternalGenealogy, FitnessTracker))
}

################################### Parameters #################################
num_sims <- 1
pop_size <- 1000
gen_no <- 100
s <- 0.3 
chiasm <- T 
fus.type <- "Y" 
mu <- 10^-4
fus.large <- F
##################### Edited Version of Evolve() Function ######################
gen_sims <- 1:4
for(sim in 1:num_sims){
  print(paste("Simulation: ", sim, collapse = ""))
  # Get a new population
  pop <- GetPop(pop_size)
  # Vector tracking PATERNAL lineage of MALES in the CURRENT population
  PaternalGenealogy <- rep(NA, pop_size)
  # Find the individuals in initial population who are male and initialize the 
  # paternal lineages
  for(i in 1:length(pop)){
    if(pop[[i]][2,13]){
      PaternalGenealogy[i] <- i
    }
  }
  FitnessTracker <- data.frame(
    "Generation" = rep(1, sum(!is.na(PaternalGenealogy))),
    "Paternal Lineage" = PaternalGenealogy[!is.na(PaternalGenealogy)],
    "Fused?" = rep(0, sum(!is.na(PaternalGenealogy))),
    "Fused SAL" = rep(NA, sum(!is.na(PaternalGenealogy))),
    "Fitness" = rep(1, sum(!is.na(PaternalGenealogy)))
  )
  names(FitnessTracker) <- c("Generation", "Paternal Lineage", "Fused?",
                             "Fused SAL", "Fitness")
  # Get a distribution of fitness effects
  dfe <- GetDFE()
  # For each generation in "gen_no"...
  for(gen in 1:gen_no){
    print(paste(c("Generation: ", gen), collapse = ""))
    # Mutate the starting population
    AoGres <- ActofGod(pop, dfe, fus.type, mu, fus.large, FitnessTracker, 
                       PaternalGenealogy, gen)
    pop <- AoGres[[1]]
    FitnessTracker <- AoGres[[2]]
    # Assess the fitness of the mutated population
    fits <- GetFit(pop, s)
    # Add fitnesses of all individuals who ARE MALE to the rows corresponding to
    # paternal lineages in the current generation
    FitnessTracker$Fitness[(FitnessTracker$Generation == gen) & 
                             (FitnessTracker$`Paternal Lineage` %in% PaternalGenealogy[!is.na(PaternalGenealogy)])] <- fits[!is.na(PaternalGenealogy)]
    # Select parents based on the fitness
    Maury.res <- Maury(pop, fits, length(pop), PaternalGenealogy)
    parents <- Maury.res[[1]]
    DadLineages <- Maury.res[[2]]
    # Select gametes from these parents
    gametes <- MakeGametes(pop, parents, chiasm=T)
    # Breed the parents
    MoL.res <- MiracleOfLife(gametes, PaternalGenealogy, DadLineages, gen, FitnessTracker)
    pop <- MoL.res[[1]]
    PaternalGenealogy <- MoL.res[[2]]
    FitnessTracker <- MoL.res[[3]]
    # Move to next gen/sim combo
    gen_sims <- gen_sims + 4
    }
}

##################################### Plotting #################################
# Change data type of Fused? column to enable plotting
edited.FitnessTracker <- FitnessTracker[,-4]
edited.FitnessTracker$`Fused?` <- FitnessTracker$`Fused?` == 1
edited.FitnessTracker$`Paternal Lineage` <- as.factor(round(FitnessTracker$`Paternal Lineage`))
# Change names for plotting purposes
# X = Generation
# Y = Fitness
# Color = Paternal lineage
# Shape = Fusion status
library(ggplot2)
ggplot(data = edited.FitnessTracker, aes(x = Generation, y = Fitness, 
                                         shape = `Fused?`, color = `Paternal Lineage`)) +
  geom_point()

# Only the descendents of male #9 survive to the final generation
Nine_FitnessTracker <- edited.FitnessTracker[edited.FitnessTracker$`Paternal Lineage` == 9,]
ggplot(data = Nine_FitnessTracker, aes(x = Generation, y = Fitness, 
                                         shape = `Fused?`, color = `Paternal Lineage`)) +
  geom_point()

# Suss: Only 6 possible fitness values in entire dataframe
levels(as.factor(FitnessTracker$Fitness))
# "0.49"   "0.595"  "0.7"    "0.7225" "0.85"   "1"   
