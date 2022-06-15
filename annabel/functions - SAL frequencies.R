# The purpose of this code is to determine whether permitting selection to
# establish sex-specific SAL ratios prior to introducing fusions enhances the
# probability of fusions fixing.

library(ggplot2)
library(viridis)

############################## Modified Functions ############################## 
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

ActofGod <- function(pop, dfe, fus.type, mu, fus.large, fuse.bool){
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
  if(fuse.bool == T){
    fuse.bool <- sample(0:1, size = 1, prob = c(0.9,0.1))
  }
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

Maury <- function(pop, fits, N = length(pop)){
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
  
  parents <- list(moms,dads)
  names(parents) <- c("Moms","Dads")
  return(parents)
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

Evolve <- function(num_sims, pop_size, gen_no, s, chiasm, fus.type, mu, fus.large, fuse.bool){
  # Initialize list of lists
  results <- vector(mode = "list", length = gen_no)
  # Create iterator for the FOUR rows corresponding to the current gensim combo
  gen_sims <- 1:4
  for(sim in 1:num_sims){
    print(paste("Simulation: ", sim, collapse = ""))
    # Get a new population
    pop <- GetPop(pop_size)
    
    # Get a distribution of fitness effects
    dfe <- GetDFE()
    
    # For each generation in "gen_no"...
    for(gen in 1:gen_no){
      print(paste(c("Generation: ", gen), collapse = ""))
      
      # Mutate the starting population
      pop <- ActofGod(pop, dfe, fus.type, mu, fus.large, fuse.bool)
      
      # Assess the fitness of the mutated population
      fits <- GetFit(pop, s)
      # Add population to output
      results[[gen]] <- pop
      # Select parents based on the fitness
      parents <- Maury(pop, fits, length(pop))
      # Select gametes from these parents
      gametes <- MakeGametes(pop, parents, chiasm=T)
      # Breed the parents
      pop <- MiracleOfLife(gametes)
      # Move to next gen/sim combo
      gen_sims <- gen_sims + 4
    }
  }
  # Output final populations after all sims and generations
  return(results)
}

########## Evolve Population for Several Generations with No Fusions ###########
num_sims <- 1
pop_size <- 1000
gen_no <- 100
s <- 0.3
chiasm <- T
fus.type <- "Y"
mu <- 10^-14
fus.large <- F
fuse.bool <- F

results <- Evolve(num_sims, pop_size, gen_no, s, chiasm, fus.type, mu, fus.large, fuse.bool)

# For each generation, count number of MALES and number of FEMALES with SAL
SAL.counter <- data.frame(
  Generation = rep(1:gen_no, each = 2),
  Sex = rep(c("F", "M"), times = 2*gen_no),
  "Count of F-Benefit SAL" = rep(0, times = 2*gen_no),
  "Count of M-Benefit SAL" = rep(0, times = 2*gen_no),
  "Number of Individuals" = rep(0, times = 2*gen_no)
)

for(g in 1:gen_no){
  for(indv in 1:pop_size){
    # Record sex and SAL of each individual
    if(sum(results[[g]][[indv]][,13]) == 0){
      SAL.counter$Number.of.Individuals[SAL.counter$Generation == g & SAL.counter$Sex == "F"] <- SAL.counter$Number.of.Individuals[SAL.counter$Generation == g & SAL.counter$Sex == "F"] + 1
      SAL.counter$Count.of.F.Benefit.SAL[SAL.counter$Generation == g & SAL.counter$Sex == "F"] <- SAL.counter$Count.of.F.Benefit.SAL[SAL.counter$Generation == g & SAL.counter$Sex == "F"] + sum(c(results[[g]][[indv]][,28] == 0, results[[g]][[indv]][,55] == 0))
      SAL.counter$Count.of.M.Benefit.SAL[SAL.counter$Generation == g & SAL.counter$Sex == "F"] <- SAL.counter$Count.of.M.Benefit.SAL[SAL.counter$Generation == g & SAL.counter$Sex == "F"] + sum(c(results[[g]][[indv]][,28] == 1, results[[g]][[indv]][,55] == 1))
    }else{
      SAL.counter$Number.of.Individuals[SAL.counter$Generation == g & SAL.counter$Sex == "M"] <- SAL.counter$Number.of.Individuals[SAL.counter$Generation == g & SAL.counter$Sex == "M"] + 1
      SAL.counter$Count.of.F.Benefit.SAL[SAL.counter$Generation == g & SAL.counter$Sex == "M"] <- SAL.counter$Count.of.F.Benefit.SAL[SAL.counter$Generation == g & SAL.counter$Sex == "M"] + sum(c(results[[g]][[indv]][,28] == 0, results[[g]][[indv]][,55] == 0))
      SAL.counter$Count.of.M.Benefit.SAL[SAL.counter$Generation == g & SAL.counter$Sex == "M"] <- SAL.counter$Count.of.M.Benefit.SAL[SAL.counter$Generation == g & SAL.counter$Sex == "M"] + sum(c(results[[g]][[indv]][,28] == 1, results[[g]][[indv]][,55] == 1))
    }
  }
}
SAL.counter$Frequency_F_Benefit <- SAL.counter$Count.of.F.Benefit.SAL/(SAL.counter$Number.of.Individuals * 4)
SAL.counter$Frequency_M_Benefit <- SAL.counter$Count.of.M.Benefit.SAL/(SAL.counter$Number.of.Individuals * 4)

plot(SAL.counter$Frequency_M_Benefit[SAL.counter$Sex == "M"]~SAL.counter$Generation[SAL.counter$Sex == "M"], 
     ylab = "Frequency of Male-Benefit Alleles in Males", xlab = "Generation",
     pch = 16, main = "1 Simulation with No Fusions")

# Check genotypic frequencies of each sexually-antagonistic locus over time
# 0 = Homozygous for FEMALE-BENEFIT allele
# 1 = Heterozygous
# 2 = Homozygous for MALE-BENEFIT allele
genotypes <- data.frame(
  Generation = rep(1:gen_no, each = pop_size),
  Individual = rep(1:pop_size, times = gen_no),
  Sex = rep(0, times = gen_no*pop_size),
  Small_SAL = rep(0, times = gen_no*pop_size),
  Large_SAL = rep(0, times = gen_no*pop_size)
)

# Record sex and genotype of each individual
for(gen in 1:gen_no){
  print(gen)
  for(indv in 1:pop_size){
    genotypes$Small_SAL[genotypes$Generation == gen & genotypes$Individual == indv] <- sum(results[[gen]][[indv]][,28])
    genotypes$Large_SAL[genotypes$Generation == gen & genotypes$Individual == indv] <- sum(results[[gen]][[indv]][,55])
    if(sum(results[[g]][[indv]][,13]) == 0){
      genotypes$Sex[genotypes$Generation == gen & genotypes$Individual == indv] <- "F"
    }else{
      genotypes$Sex[genotypes$Generation == gen & genotypes$Individual == indv] <- "M"
    }
  }
}

# Genotypic distribution of small autosome SAL in males at generation 0
plot(density(genotypes$Small_SAL[genotypes$Sex == "M" &
                                               genotypes$Generation == 1]),
     main = "Distribution of Genotypes of Small Autosome in Males (First Generation)")
# Genotypic distribution of small autosome SAL in males at final generation
plot(density(genotypes$Small_SAL[genotypes$Sex == "M" &
                                               genotypes$Generation == gen_no]),
     main = "Distribution of Genotypes of Small Autosome in Males (Final Generation)")

# Find frequency of each genotype in each generation
genotypic_frequencies <- data.frame(
  Generation = rep(1:gen_no, each = 6),
  Autosome = rep(c("Small", "Large"), times = gen_no, each = 3),
  Genotype = rep(c("F-Benefit Homo", "Hetero", "M-Benefit Homo"), times = gen_no*2),
  Frequency = rep(0, times = gen_no*6),
  "Frequency in Males" = rep(0, times = gen_no*6),
  "Frequency in Females" = rep(0, times = gen_no*6)
  )

for(gen in 1:gen_no){
  genotypic_frequencies$Frequency[genotypic_frequencies$Generation == gen & genotypic_frequencies$Autosome == "Small" & genotypic_frequencies$Genotype == "F-Benefit Homo"] <- sum(genotypes$Small_SAL[genotypes$Generation == gen] == 0)/pop_size
  genotypic_frequencies$Frequency[genotypic_frequencies$Generation == gen & genotypic_frequencies$Autosome == "Small" & genotypic_frequencies$Genotype == "Hetero"] <- sum(genotypes$Small_SAL[genotypes$Generation == gen] == 1)/pop_size
  genotypic_frequencies$Frequency[genotypic_frequencies$Generation == gen & genotypic_frequencies$Autosome == "Small" & genotypic_frequencies$Genotype == "M-Benefit Homo"] <- sum(genotypes$Small_SAL[genotypes$Generation == gen] == 2)/pop_size
  genotypic_frequencies$Frequency[genotypic_frequencies$Generation == gen & genotypic_frequencies$Autosome == "Large" & genotypic_frequencies$Genotype == "F-Benefit Homo"] <- sum(genotypes$Large_SAL[genotypes$Generation == gen] == 0)/pop_size
  genotypic_frequencies$Frequency[genotypic_frequencies$Generation == gen & genotypic_frequencies$Autosome == "Large" & genotypic_frequencies$Genotype == "Hetero"] <- sum(genotypes$Large_SAL[genotypes$Generation == gen] == 1)/pop_size
  genotypic_frequencies$Frequency[genotypic_frequencies$Generation == gen & genotypic_frequencies$Autosome == "Large" & genotypic_frequencies$Genotype == "M-Benefit Homo"] <- sum(genotypes$Large_SAL[genotypes$Generation == gen] == 2)/pop_size
  
  num_males <- sum(genotypes$Sex[genotypes$Generation == gen] == "M") 
  genotypic_frequencies$Frequency.in.Males[genotypic_frequencies$Generation == gen & genotypic_frequencies$Autosome == "Small" & genotypic_frequencies$Genotype == "F-Benefit Homo"] <- sum(genotypes$Small_SAL[genotypes$Generation == gen] == 0)/num_males
  genotypic_frequencies$Frequency.in.Males[genotypic_frequencies$Generation == gen & genotypic_frequencies$Autosome == "Small" & genotypic_frequencies$Genotype == "Hetero"] <- sum(genotypes$Small_SAL[genotypes$Generation == gen] == 1)/num_males
  genotypic_frequencies$Frequency.in.Males[genotypic_frequencies$Generation == gen & genotypic_frequencies$Autosome == "Small" & genotypic_frequencies$Genotype == "M-Benefit Homo"] <- sum(genotypes$Small_SAL[genotypes$Generation == gen] == 2)/num_males
  genotypic_frequencies$Frequency.in.Males[genotypic_frequencies$Generation == gen & genotypic_frequencies$Autosome == "Large" & genotypic_frequencies$Genotype == "F-Benefit Homo"] <- sum(genotypes$Large_SAL[genotypes$Generation == gen] == 0)/num_males
  genotypic_frequencies$Frequency.in.Males[genotypic_frequencies$Generation == gen & genotypic_frequencies$Autosome == "Large" & genotypic_frequencies$Genotype == "Hetero"] <- sum(genotypes$Large_SAL[genotypes$Generation == gen] == 1)/num_males
  genotypic_frequencies$Frequency.in.Males[genotypic_frequencies$Generation == gen & genotypic_frequencies$Autosome == "Large" & genotypic_frequencies$Genotype == "M-Benefit Homo"] <- sum(genotypes$Large_SAL[genotypes$Generation == gen] == 2)/num_males
  
  num_females <- sum(genotypes$Sex[genotypes$Generation == gen] == "F")
  genotypic_frequencies$Frequency.in.Females[genotypic_frequencies$Generation == gen & genotypic_frequencies$Autosome == "Small" & genotypic_frequencies$Genotype == "F-Benefit Homo"] <- sum(genotypes$Small_SAL[genotypes$Generation == gen] == 0)/num_females
  genotypic_frequencies$Frequency.in.Females[genotypic_frequencies$Generation == gen & genotypic_frequencies$Autosome == "Small" & genotypic_frequencies$Genotype == "Hetero"] <- sum(genotypes$Small_SAL[genotypes$Generation == gen] == 1)/num_females
  genotypic_frequencies$Frequency.in.Females[genotypic_frequencies$Generation == gen & genotypic_frequencies$Autosome == "Small" & genotypic_frequencies$Genotype == "M-Benefit Homo"] <- sum(genotypes$Small_SAL[genotypes$Generation == gen] == 2)/num_females
  genotypic_frequencies$Frequency.in.Females[genotypic_frequencies$Generation == gen & genotypic_frequencies$Autosome == "Large" & genotypic_frequencies$Genotype == "F-Benefit Homo"] <- sum(genotypes$Large_SAL[genotypes$Generation == gen] == 0)/num_females
  genotypic_frequencies$Frequency.in.Females[genotypic_frequencies$Generation == gen & genotypic_frequencies$Autosome == "Large" & genotypic_frequencies$Genotype == "Hetero"] <- sum(genotypes$Large_SAL[genotypes$Generation == gen] == 1)/num_females
  genotypic_frequencies$Frequency.in.Females[genotypic_frequencies$Generation == gen & genotypic_frequencies$Autosome == "Large" & genotypic_frequencies$Genotype == "M-Benefit Homo"] <- sum(genotypes$Large_SAL[genotypes$Generation == gen] == 2)/num_females
  
  }
# Plot frequency of each genotype over time
small_frequencies <- genotypic_frequencies[genotypic_frequencies$Autosome == "Small",]
names(small_frequencies)[3] <- "Genotype of Small Autosome SAL"
ggplot() +
  geom_line(data = small_frequencies, aes(x = Generation, y = Frequency, 
                                          color = `Genotype of Small Autosome SAL`)) +
  ggtitle(paste(c("Chiasmatic, no fusions. Population Size = ", pop_size, " mu = ", mu, " s = ", s), collapse = ""))

large_frequencies <- genotypic_frequencies[genotypic_frequencies$Autosome == "Large",]
names(large_frequencies)[3] <- "Genotype of Large Autosome SAL"
ggplot() +
  geom_line(data = large_frequencies, aes(x = Generation, y = Frequency, 
                                          color = `Genotype of Large Autosome SAL`)) +
  ggtitle(paste(c("Chiasmatic, no fusions. Population Size = ", pop_size, " mu = ", mu, " s = ", s), collapse = ""))
