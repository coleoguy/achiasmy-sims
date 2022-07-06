# Goals: Determine whether...
  # 1. individuals with fusions have more rapid fitness decay than individuals without
  # 2. SAL genotype of fused chromosome influences rapidity of fitness decay

################################## Hypotheses ##################################
### Hypothesis A: Fusions are not fixing because the decay in fitness for even 
  # male-benefit fusions is extremely deleterious

  # H1: Lineages with purportedly beneficial fusions still decay more rapidly 
  # than lineages with no fusions, so DFE must be made less severe

  # H0: Lineages with beneficial fusions have approximately equal decay in 
  # fitness as compared to lineages with __no__ fusions, so increased extinction 
  # amongst fused lineages is merely a function of relative rarity of these lineages


### Hypothesis B: Fusions are not fixing because Muller's ratchet is so 
  # deleterious that males with male-benefit fusions have equivalent decay in 
  # fitness as males with female-benefit fusions

  # H1: Lineages with fusions of female-benefit SAL to Y have more rapid 
  # decrease in fitness than lineages with fusions of male-benefit SAL to Y

  # H0: New mutations are so deleterious that SAL allele does not override the 
  # deleterious effects of Muller's ratchet and both male-benefit to Y fusions 
  # and female-benefit to Y fusions have the same fitness consequences

######## Create edited versions of functions needed for data collection ########
# Read initial functions
source("functions.R")
# Define "not in" operator
`%!in%` <- Negate(`%in%`)

# Edit ActofGod to record fusions AND SAL genotype of the new fusions as they occur
ActofGod <- function(pop, dfe, fus.type, mu.table, fus.large, FitnessTracker, 
                     PaternalGenealogy, gen){
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
      # has a fusion, then introduce fusion to Y and record that this male's
      # lineage has a fusion
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
# Edit Maury to track lineage offshoots as DISTINCT lineages
Maury <- function(pop, fits, N = length(pop), PaternalGenealogy){
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
# Edit MiracleofLife so lineage corresponding to new males is recorded AND
# so genotype of the sexually-antagonistic allele of Y-fused small chromosomes
# is recorded
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

################## Initialize parameters for run of Evolve() ###################
pop_size <- 1000
gen_no <- 1000
mu <- 0.00000001
fus.type <- "Y"
fus.large <- F
s <- 1

########################## Initialize Data Collection ##########################
pop <- GetPop(pop_size)
# Vector tracking PATERNAL lineage of MALES in the CURRENT population
PaternalGenealogy <- rep(NA, 1000)
# Find the individuals in initial population who are male and initialize the 
# paternal lineages
for(i in 1:length(pop)){
  if(pop[[i]][2,13]){
    PaternalGenealogy[i] <- i
  }
}

# Initialize dataframe to collect... 
# 1. Generation: the generation in which the male existed
# 2. Paternal Lineage: the unique lineage identifier of each male
# 3. Fused?: the fusion status of the male
# 4. Fused SAL: the genotype of the sexually antagonistic allele of the small
  #  chromosome fused to this male's Y (if applicable)
# 5. Fitness: the relative fitness of the male
FitnessTracker <- data.frame(
  "Generation" = rep(1, sum(!is.na(PaternalGenealogy))),
  "Paternal Lineage" = PaternalGenealogy[!is.na(PaternalGenealogy)],
  "Fused?" = rep(0, sum(!is.na(PaternalGenealogy))),
  "Fused SAL" = rep(NA, sum(!is.na(PaternalGenealogy))),
  "Fitness" = rep(1, sum(!is.na(PaternalGenealogy)))
)
names(FitnessTracker) <- c("Generation", "Paternal Lineage", "Fused?",
                           "Fused SAL", "Fitness")

################################# Run Evolve() #################################
dfe <- GetDFE()

for(gen in 1:gen_no){
  print(paste(c("Generation: ", gen), collapse = ""))
  
  # Mutate the starting population
  AoGres <- ActofGod(pop, dfe, fus.type, mu.table, fus.large, FitnessTracker, 
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
}

##################### Organize Data for Slope Calculations #####################
# Read in FitnessTracker from Terra
FitnessTracker <- readRDS("FitnessTracker.rds")
# Remove ALL rows corresponding to final generation (fitnesses in these rows
# are NA)
FitnessTracker <- FitnessTracker[!is.na(FitnessTracker$Fitness),]
  
# Create subsets of FitnessTracker for each fusion class 
# (No Fusion, Male-Benefit, Female-Benefit)
NF.FitnessTracker <- FitnessTracker[FitnessTracker$`Fused?` == 0,]
Fus.FitnessTracker <- FitnessTracker[(FitnessTracker$`Fused?` == 1) & 
                                       !is.na(FitnessTracker$`Fused SAL`),]

MB.FitnessTracker <- Fus.FitnessTracker[Fus.FitnessTracker$`Fused SAL` == 1,]
FB.FitnessTracker <- Fus.FitnessTracker[Fus.FitnessTracker$`Fused SAL` == 0,]

### Calculate slopes for lineages with male-benefit fusions
# find all unique male-benefit lineages
MB.lineages <- levels(as.factor(MB.FitnessTracker$`Paternal Lineage`))

# initialize Tips and HigherOrder
Tips <- c()
HigherOrder <- c()

# initialize number to report the number of decimals in the most recent lineage
MostDecimals <- 0

# record a complete lineage for each offshoot
for(i in 1:gen_no){
  # check if this origin ever got a male-benefit fusion and, if so, group 
  # together all lineages with this origin
  SameOrigin <- c()
  for(m in 1:length(MB.lineages)){
    if(substr(MB.lineages[m], 1, nchar(i) + 1) == paste(c(i, "."), collapse = "")){
      SameOrigin <- append(SameOrigin, MB.lineages[m])
    }
  }
  if(length(SameOrigin)){
    # Find number of decimal places present in most recent lineages 
    for(x in 1:length(SameOrigin)){
      if(nchar(strsplit(as.character(SameOrigin[x]), "\\.")[[1]][2]) > MostDecimals){
        MostDecimals <- nchar(strsplit(as.character(SameOrigin[x]), "\\.")[[1]][2])
      }
    }
    # Find the lineage(s) of this origin which were created most recently and
    # store in Tips
    for(x in 1:length(SameOrigin)){
      if(nchar(strsplit(as.character(SameOrigin[x]), "\\.")[[1]][2]) == MostDecimals){
        Tips <- append(Tips, SameOrigin[x])
      }
    }
    MostDecimals <- 0
    # check if ORIGINAL origin is present and, if so, add this to SameOrigin
    # vector
    if(i %in% MB.lineages){
      SameOrigin <- append(SameOrigin, i)
    }
    # Store all the lineages which CREATED these lineages in HigherOrder
    HigherOrder <- SameOrigin[SameOrigin %!in% Tips]
    
    # If there are HigherOrder lineages...
    if(length(HigherOrder)){
      # For each lineage in HigherOrder...
      for(p in 1:length(HigherOrder)){
        # REPLICATE all MB.FitnessTracker rows where the HigherOrder lineage is present
        # as many times as there are tips
        DuplicatedRows <- c()
        for(t in 1:length(Tips)){
          if(length(grep(HigherOrder[p], substr(Tips[t],1,nchar(Tips[t])-1)))){
            DuplicatedRows <- rbind(DuplicatedRows, MB.FitnessTracker[MB.FitnessTracker$`Paternal Lineage` == HigherOrder[p],])
            # Rename the duplicated HigherOrder rows with the name in Tips
            DuplicatedRows$`Paternal Lineage`[DuplicatedRows$`Paternal Lineage` == HigherOrder[p]] <- rep(Tips[t], sum(DuplicatedRows$`Paternal Lineage` == HigherOrder[p]))
            
          }
        }
        # Add the replicated + renamed rows to MB.FitnessTracker and delete the 
        # original HigherOrder rows from MB.FitnessTracker
        MB.FitnessTracker <- rbind(MB.FitnessTracker, DuplicatedRows)
        MB.FitnessTracker <- subset(MB.FitnessTracker, `Paternal Lineage` != HigherOrder[p])
      }
    }
    Tips <- c()
    HigherOrder <- c()
  }
}

### Calculate slopes for lineages with female-benefit fusions
# find all unique female-benefit lineages
FB.lineages <- levels(as.factor(FB.FitnessTracker$`Paternal Lineage`))

# initialize Tips and HigherOrder
Tips <- c()
HigherOrder <- c()

# initialize number to report the number of decimals in the most recent lineage
MostDecimals <- 0

# record a complete lineage for each offshoot
for(i in 1:gen_no){
  # check if this origin ever got a female-benefit fusion and, if so, group 
  # together all lineages with this origin
  SameOrigin <- c()
  for(f in 1:length(FB.lineages)){
    if(substr(FB.lineages[f], 1, nchar(i) + 1) == paste(c(i, "."), collapse = "")){
      SameOrigin <- append(SameOrigin, FB.lineages[f])
    }
  }
  if(length(SameOrigin)){
    # Find number of decimal places present in most recent lineages 
    for(x in 1:length(SameOrigin)){
      if(nchar(strsplit(as.character(SameOrigin[x]), "\\.")[[1]][2]) > MostDecimals){
        MostDecimals <- nchar(strsplit(as.character(SameOrigin[x]), "\\.")[[1]][2])
      }
    }
    # Find the lineage(s) of this origin which were created most recently and
    # store in Tips
    for(x in 1:length(SameOrigin)){
      if(nchar(strsplit(as.character(SameOrigin[x]), "\\.")[[1]][2]) == MostDecimals){
        Tips <- append(Tips, SameOrigin[x])
      }
    }
    MostDecimals <- 0
    # check if ORIGINAL origin is present and, if so, add this to SameOrigin
    # vector
    if(i %in% FB.lineages){
      SameOrigin <- append(SameOrigin, i)
    }
    # Store all the lineages which CREATED these lineages in HigherOrder
    HigherOrder <- SameOrigin[SameOrigin %!in% Tips]
    
    # If there are HigherOrder lineages...
    if(length(HigherOrder)){
      # For each lineage in HigherOrder...
      for(p in 1:length(HigherOrder)){
        # REPLICATE all FB.FitnessTracker rows where the HigherOrder lineage is present
        # as many times as there are tips
        DuplicatedRows <- c()
        for(t in 1:length(Tips)){
          if(length(grep(HigherOrder[p], substr(Tips[t],1,nchar(Tips[t])-1)))){
            DuplicatedRows <- rbind(DuplicatedRows, FB.FitnessTracker[FB.FitnessTracker$`Paternal Lineage` == HigherOrder[p],])
            # Rename the duplicated HigherOrder rows with the name in Tips
            DuplicatedRows$`Paternal Lineage`[DuplicatedRows$`Paternal Lineage` == HigherOrder[p]] <- rep(Tips[t], sum(DuplicatedRows$`Paternal Lineage` == HigherOrder[p]))
            
          }
        }
        # Add the replicated + renamed rows to FB.FitnessTracker and delete the 
        # original HigherOrder rows from FB.FitnessTracker
        FB.FitnessTracker <- rbind(FB.FitnessTracker, DuplicatedRows)
        FB.FitnessTracker <- subset(FB.FitnessTracker, `Paternal Lineage` != HigherOrder[p])
      }
    }
    Tips <- c()
    HigherOrder <- c()
  }
}


### Calculate slopes for lineages with no fusions
# find all unique no fusion lineages
NF.lineages <- levels(as.factor(NF.FitnessTracker$`Paternal Lineage`))

# initialize Tips and HigherOrder
Tips <- c()
HigherOrder <- c()

# initialize number to report the number of decimals in the most recent lineage
MostDecimals <- 0

# record a complete lineage for each offshoot
for(i in 1:gen_no){
  # check if this origin never got a fusion and, if so, group together all 
  # lineages with this origin
  SameOrigin <- c()
  for(n in 1:length(NF.lineages)){
    if(substr(NF.lineages[n], 1, nchar(i) + 1) == paste(c(i, "."), collapse = "")){
      SameOrigin <- append(SameOrigin, NF.lineages[n])
    }
  }
  if(length(SameOrigin)){
    # Find number of decimal places present in most recent lineages 
    for(x in 1:length(SameOrigin)){
      if(nchar(strsplit(as.character(SameOrigin[x]), "\\.")[[1]][2]) > MostDecimals){
        MostDecimals <- nchar(strsplit(as.character(SameOrigin[x]), "\\.")[[1]][2])
      }
    }
    # Find the lineage(s) of this origin which were created most recently and
    # store in Tips
    for(x in 1:length(SameOrigin)){
      if(nchar(strsplit(as.character(SameOrigin[x]), "\\.")[[1]][2]) == MostDecimals){
        Tips <- append(Tips, SameOrigin[x])
      }
    }
    MostDecimals <- 0
    # check if ORIGINAL origin is present and, if so, add this to SameOrigin
    # vector
    if(i %in% NF.lineages){
      SameOrigin <- append(SameOrigin, i)
    }
    # Store all the lineages which CREATED these lineages in HigherOrder
    HigherOrder <- SameOrigin[SameOrigin %!in% Tips]
    
    # If there are HigherOrder lineages...
    if(length(HigherOrder)){
      # For each lineage in HigherOrder...
      for(p in 1:length(HigherOrder)){
        # REPLICATE all NF.FitnessTracker rows where the HigherOrder lineage is present
        # as many times as there are tips
        DuplicatedRows <- c()
        for(t in 1:length(Tips)){
          if(length(grep(HigherOrder[p], substr(Tips[t],1,nchar(Tips[t])-1)))){
            DuplicatedRows <- rbind(DuplicatedRows, NF.FitnessTracker[NF.FitnessTracker$`Paternal Lineage` == HigherOrder[p],])
            # Rename the duplicated HigherOrder rows with the name in Tips
            DuplicatedRows$`Paternal Lineage`[DuplicatedRows$`Paternal Lineage` == HigherOrder[p]] <- rep(Tips[t], sum(DuplicatedRows$`Paternal Lineage` == HigherOrder[p]))
            
          }
        }
        # Add the replicated + renamed rows to NF.FitnessTracker and delete the 
        # original HigherOrder rows from NF.FitnessTracker
        NF.FitnessTracker <- rbind(NF.FitnessTracker, DuplicatedRows)
        NF.FitnessTracker <- subset(NF.FitnessTracker, `Paternal Lineage` != HigherOrder[p])
      }
    }
    Tips <- c()
    HigherOrder <- c()
  }
}

############################### Calculate Slopes ############################### 
# Read in data from Terra run
MB.FitnessTracker <- readRDS("MB.FitnessTracker.rds")
FB.FitnessTracker <- readRDS("FB.FitnessTracker.rds")
NF.FitnessTracker <- readRDS("NF.FitnessTracker.rds")
### Male-Benefit
# Create factor of all the TERMINAL male-benefit lineages
MB.Tips <- levels(as.factor(MB.FitnessTracker$`Paternal Lineage`))
# Initialize a dataframe in which to store slopes
  # Lineage = Unique identifier of the terminal lineage
  # First Fitness = Fitness of the lineage in the first generation in which the
    # lineage receives a fusion
  # Last Fitness = Fitness of the lineage in the final generation in which the
    # lineage appears
  # First Generation = Generation in which first member of lienage receives 
    # male-benefit fusion
  # Last Generation = Final generation in which this lineage appears
  # Slope = Change in fitness over change in generations
MB.Slopes <- data.frame(
  Lineage = MB.Tips,
  a = rep(NA, length(MB.Tips)),
  b = rep(NA, length(MB.Tips)),
  c = rep(NA, length(MB.Tips)),
  d = rep(NA, length(MB.Tips)),
  e = rep(NA, length(MB.Tips))
)
names(MB.Slopes) <- c("Lineage", "First Fitness", "Last Fitness", 
                      "First Generation", "Last Generation", "Slope")

# For each terminal male-benefit lineage...
for(i in 1:length(MB.Tips)){
  # Get rows of all entries of this lineage
  cur.lineage <- MB.FitnessTracker[
    MB.FitnessTracker$`Paternal Lineage` == MB.Tips[i],
  ]
  # Store first and last for both generation and fitness
  MB.Slopes$`First Generation`[i] <- cur.lineage$Generation[which.min(cur.lineage$Generation)]
  MB.Slopes$`First Fitness`[i] <- cur.lineage$Fitness[which.min(cur.lineage$Generation)]
  
  MB.Slopes$`Last Generation`[i] <- cur.lineage$Generation[which.max(cur.lineage$Generation)]
  MB.Slopes$`Last Fitness`[i] <- cur.lineage$Fitness[which.max(cur.lineage$Generation)]
  # Calculate the difference between first and last fitness and first and last
    # generation and store slope
  MB.Slopes$Slope[i] <- (MB.Slopes$`Last Fitness`[i] - MB.Slopes$`First Fitness`[i])/(MB.Slopes$`Last Generation`[i] - MB.Slopes$`First Generation`[i])
}
# Retain only non-NA slopes
MB.Slopes.NotNA <- MB.Slopes[!is.na(MB.Slopes$Slope),]

### Female-Benefit
# Create factor of all the TERMINAL female-benefit lineages
FB.Tips <- levels(as.factor(FB.FitnessTracker$`Paternal Lineage`))
# Initialize a dataframe in which to store slopes
# Lineage = Unique identifier of the terminal lineage
# First Fitness = Fitness of the lineage in the first generation in which the
# lineage receives a fusion
# Last Fitness = Fitness of the lineage in the final generation in which the
# lineage appears
# First Generation = Generation in which first member of lienage receives 
# female-benefit fusion
# Last Generation = Final generation in which this lineage appears
# Slope = Change in fitness over change in generations
FB.Slopes <- data.frame(
  Lineage = FB.Tips,
  a = rep(NA, length(FB.Tips)),
  b = rep(NA, length(FB.Tips)),
  c = rep(NA, length(FB.Tips)),
  d = rep(NA, length(FB.Tips)),
  e = rep(NA, length(FB.Tips))
)
names(FB.Slopes) <- c("Lineage", "First Fitness", "Last Fitness", 
                      "First Generation", "Last Generation", "Slope")

# For each terminal female-benefit lineage...
for(i in 1:length(FB.Tips)){
  # Get rows of all entries of this lineage
  cur.lineage <- FB.FitnessTracker[
    FB.FitnessTracker$`Paternal Lineage` == FB.Tips[i],
  ]
  # Store first and last for both generation and fitness
  FB.Slopes$`First Generation`[i] <- cur.lineage$Generation[which.min(cur.lineage$Generation)]
  FB.Slopes$`First Fitness`[i] <- cur.lineage$Fitness[which.min(cur.lineage$Generation)]
  
  FB.Slopes$`Last Generation`[i] <- cur.lineage$Generation[which.max(cur.lineage$Generation)]
  FB.Slopes$`Last Fitness`[i] <- cur.lineage$Fitness[which.max(cur.lineage$Generation)]
  # Calculate the difference between first and last fitness and first and last
  # generation and store slope
  FB.Slopes$Slope[i] <- (FB.Slopes$`Last Fitness`[i] - FB.Slopes$`First Fitness`[i])/(FB.Slopes$`Last Generation`[i] - FB.Slopes$`First Generation`[i])
}
# Retain only non-NA slopes
FB.Slopes.NotNA <- FB.Slopes[!is.na(FB.Slopes$Slope),]

### No Fusion
# Create factor of all the TERMINAL no fusion lineages
NF.Tips <- levels(as.factor(NF.FitnessTracker$`Paternal Lineage`))
# Initialize a dataframe in which to store slopes
# Lineage = Unique identifier of the terminal lineage
# First Fitness = Fitness of the lineage in the first generation in which the
# lineage receives a fusion
# Last Fitness = Fitness of the lineage in the final generation in which the
# lineage appears
# First Generation = Generation in which first member of lienage receives 
# no fusion fusion
# Last Generation = Final generation in which this lineage appears
# Slope = Change in fitness over change in generations
NF.Slopes <- data.frame(
  Lineage = NF.Tips,
  a = rep(NA, length(NF.Tips)),
  b = rep(NA, length(NF.Tips)),
  c = rep(NA, length(NF.Tips)),
  d = rep(NA, length(NF.Tips)),
  e = rep(NA, length(NF.Tips))
)
names(NF.Slopes) <- c("Lineage", "First Fitness", "Last Fitness", 
                      "First Generation", "Last Generation", "Slope")

# For each terminal no fusion lineage...
for(i in 1:length(NF.Tips)){
  # Get rows of all entries of this lineage
  cur.lineage <- NF.FitnessTracker[
    NF.FitnessTracker$`Paternal Lineage` == NF.Tips[i],
  ]
  # Store first and last for both generation and fitness
  NF.Slopes$`First Generation`[i] <- cur.lineage$Generation[which.min(cur.lineage$Generation)]
  NF.Slopes$`First Fitness`[i] <- cur.lineage$Fitness[which.min(cur.lineage$Generation)]
  
  NF.Slopes$`Last Generation`[i] <- cur.lineage$Generation[which.max(cur.lineage$Generation)]
  NF.Slopes$`Last Fitness`[i] <- cur.lineage$Fitness[which.max(cur.lineage$Generation)]
  # Calculate the difference between first and last fitness and first and last
  # generation and store slope
  NF.Slopes$Slope[i] <- (NF.Slopes$`Last Fitness`[i] - NF.Slopes$`First Fitness`[i])/(NF.Slopes$`Last Generation`[i] - NF.Slopes$`First Generation`[i])
}
# Retain only non-NA slopes
NF.Slopes.NotNA <- NF.Slopes[!is.na(NF.Slopes$Slope),]

################################ Run Statistics ################################
# Statistical Test A: Check if decay in fitness is greater for lineages with
# male-benefit fusions as opposed to lineages with no fusions
t.test(MB.Slopes.NotNA$Slope, NF.Slopes.NotNA$Slope)
# Statistical Test B: Check if decay in fitness is greater for lineages with
# female-benefit fusions as opposed to lineages with male-benefit fusions
t.test(abs(MB.Slopes.NotNA$Slope), abs(FB.Slopes.NotNA$Slope))

######################## Testing: Testing of Male-Benefit Code #######################
# Testing whether male benefit code correctly splices out appropriate lineages
test.full <- data.frame(
  a = c(1, 1, 2, 2, 2, 3, 3, 4, 4),
  b = c(1, 2, 2.1, 2.2, 1, 1.1, 1.2, 1.11, 1.12),
  c = c(0, 0, 1, 1, 1, 1, 1, 1, 1),
  d = c(NA, NA, 1, 1, 1,1,1,1,1),
  e = c(1, 1, 0.9, 0.8, 1, 0.95, 0.85, 0.9, 0.8)
)
colnames(test.full) <- c("Generation", "Paternal Lineage", "Fused?",
                         "Fused SAL", "Fitness")
test.MB <- test.full[!is.na(test.full$`Fused SAL`),]
test.MB <- test.MB[test.MB$`Fused SAL` == 1,]

MB.lineages <- levels(as.factor(test.MB$`Paternal Lineage`))

# initialize Tips and HigherOrder
Tips <- c()
HigherOrder <- c()

# initialize number to report the number of decimals in the most recent lineage
MostDecimals <- 0

# record a complete lineage for each offshoot
gen_no <- 4
for(i in 1:gen_no){
  # check if this origin ever got a male-benefit fusion
  origin <- paste(c(i, "."), collapse = "")
  if(length(grep(origin, MB.lineages))){
    # if so, group together all lineages with this origin
    SameOrigin <- MB.lineages[grep(origin, MB.lineages)]
    
    # Find number of decimal places present in most recent lineages 
    for(x in 1:length(SameOrigin)){
      if(nchar(strsplit(as.character(SameOrigin[x]), "\\.")[[1]][2]) > MostDecimals){
        MostDecimals <- nchar(strsplit(as.character(SameOrigin[x]), "\\.")[[1]][2])
      }
    }
    # Find the lineage(s) of this origin which were created most recently and
    # store in Tips
    for(x in 1:length(SameOrigin)){
      if(nchar(strsplit(as.character(SameOrigin[x]), "\\.")[[1]][2]) == MostDecimals){
        Tips <- append(Tips, SameOrigin[x])
      }
    }
    MostDecimals <- 0
    # check if ORIGINAL origin is present and, if so, add this to SameOrigin
    # vector
    if(i %in% MB.lineages){
      SameOrigin <- append(SameOrigin, i)
    }
    # Store all the lineages which CREATED these lineages in HigherOrder
    HigherOrder <- SameOrigin[SameOrigin %!in% Tips]
    
    # If there are HigherOrder lineages...
    if(length(HigherOrder)){
      # For each lineage in HigherOrder...
      for(p in 1:length(HigherOrder)){
        # REPLICATE all test.MB rows where the HigherOrder lineage is present
        # as many times as there are tips
        DuplicatedRows <- c()
        for(t in 1:length(Tips)){
          if(length(grep(HigherOrder[p], substr(Tips[t],1,nchar(Tips[t])-1)))){
            DuplicatedRows <- rbind(DuplicatedRows, test.MB[test.MB$`Paternal Lineage` == HigherOrder[p],])
            # Rename the duplicated HigherOrder rows with the name in Tips
            DuplicatedRows$`Paternal Lineage`[DuplicatedRows$`Paternal Lineage` == HigherOrder[p]] <- rep(Tips[t], sum(DuplicatedRows$`Paternal Lineage` == HigherOrder[p]))
            
          }
        }
        # Add the replicated + renamed rows to test.MB and delete the 
        # original HigherOrder rows from test.MB
        test.MB <- rbind(test.MB, DuplicatedRows)
        test.MB <- subset(test.MB, `Paternal Lineage` != HigherOrder[p])
      }
    }
    Tips <- c()
    HigherOrder <- c()
  }
}

