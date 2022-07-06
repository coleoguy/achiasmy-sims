# Goal: Determine whether fusions are simply not common enough to fix

################################## Hypotheses ##################################
# Hypothesis: Lineages which receive fusions are more likely to go extinct than 
# lineages without fusions
    # If true, then back to the drawing board
# Null Hypothesis: Lineages which receive fusions are not more likely to go 
# extinct than lineages without fusions
    # If true, then simply increase frequency of fusions

######## Create edited versions of functions needed for data collection ########
# Read initial functions
source("functions.R")
# Define "not in" operator
`%!in%` <- Negate(`%in%`)

# Edit ActofGod so fusions to Y are recorded in LineageStatus matrix
# AND give each generation a 95% chance of having a single individual get a 
# fusion
ActofGod <- function(pop, dfe, fus.type, mu.table, fus.large, LineageStatus, PaternalGenealogy){
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
  
  # Each generation has 95% chance of having 1 individual with a fusion
  fuse.bool <- sample(0:1, size = 1, prob = c(0.05,0.95))
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
        LineageStatus$`Fused?`[LineageStatus$`Paternal Lineage` == PaternalGenealogy[fuse.indv]] <- 1
        # introduce large or small fusion based on value of "fus.large"
        if(fus.large){
          pop[[fuse.indv]][2, 25] <- 2
        }else{
          pop[[fuse.indv]][2, 25] <- 1
        }
      }
    }
  }
  
  return(list(pop, LineageStatus))
}
# Edit Maury so extinct paternal lineages are recorded
Maury <- function(pop, fits, N = length(pop), LineageStatus, PaternalGenealogy){
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
  # Record the lineages of all individuals who were NOT selected as dads as extinct
  incels <- PaternalGenealogy[PaternalGenealogy %!in% dads & !is.na(PaternalGenealogy)]
  LineageStatus$`Extinct?`[LineageStatus$`Paternal Lineage` %in% incels] <- 1
  
  # Record the lineages of all individuals who were selected as dads
  DadLineages <- PaternalGenealogy[dads]
    
  parents <- list(moms,dads)
  names(parents) <- c("Moms","Dads")
  return(list(parents, LineageStatus, DadLineages))
}
# Edit MiracleofLife so lineage corresponding to new males is recorded
MiracleOfLife <- function(gametes, PaternalGenealogy, DadLineages){
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
    
    # If sperm is Y-bearing, add the dad's lineage to corresponding element of PaternalLineages
    if(cur.sperm[13]){
      PaternalGenealogy[i] <- DadLineages[s]
    # If sperm is X-bearing, add an NA to corresponding element of PaternalLineages
    }else{
      PaternalGenealogy[i] <- NA
    }
    
    gametes$eggs <- gametes$eggs[-e]
    gametes$sperm <- gametes$sperm[-s]
  }
  
  return(list(newgen, PaternalGenealogy))
}

################## Initialize parameters for run of Evolve() ###################
pop_size <- 1000
gen_no <- 1000
mu <- 0.00000001
fus.type <- "Y"
fus.large <- F
s <- 1

########################## Initialize Data Collection ##########################
# Initialize matrix tracking whether a lineage has gone extinct or received a
# fusion to its Y
pop <- GetPop(pop_size)
# Vector tracking PATERNAL lineage of MALES in the CURRENT population
PaternalGenealogy <- rep(NA, 100)
# Find the individuals in initial population who are male and initialize the 
# paternal lineages
for(i in 1:length(pop)){
  if(pop[[i]][2,13]){
    PaternalGenealogy[i] <- i
  }
}
n.males <- sum(!is.na(PaternalGenealogy))
LineageStatus <- data.frame(matrix(c(PaternalGenealogy[!is.na(PaternalGenealogy)], 
                                     rep(0, n.males), rep(0, n.males)), 
                                   nrow = n.males, ncol = 3))
colnames(LineageStatus) <- c("Paternal Lineage", "Extinct?", "Fused?")

################################# Run Evolve() #################################
dfe <- GetDFE()
gen_rows <- 1:4

for(gen in 1:gen_no){
  print(paste(c("Generation: ", gen), collapse = ""))
  
  # Mutate the starting population
  AoGres <- ActofGod(pop, dfe, fus.type, mu.table, fus.large, LineageStatus, PaternalGenealogy)
  pop <- AoGres[[1]]
  LineageStatus <- AoGres[[2]]
  
  # Assess the fitness of the mutated population
  fits <- GetFit(pop, s)
  
  # Select parents based on the fitness
  Maury.res <- Maury(pop, fits, length(pop), LineageStatus, PaternalGenealogy)
  parents <- Maury.res[[1]]
  LineageStatus <- Maury.res[[2]]
  DadLineages <- Maury.res[[3]]
  
  # Select gametes from these parents
  gametes <- MakeGametes(pop, parents, chiasm=T)
  
  # Breed the parents
  MoL.res <- MiracleOfLife(gametes, PaternalGenealogy, DadLineages)
  pop <- MoL.res[[1]]
  PaternalGenealogy <- MoL.res[[2]]
  
  # Move to next two generation's rows
  gen_rows <- gen_rows + 4
}

# Edit data for chi-square test
LineageExtinction <- matrix(nrow = 2, ncol = 2)
LineageExtinction[1,1] <- sum(LineageStatus$`Extinct?` & LineageStatus$`Fused?`)
LineageExtinction[1,2] <- sum(LineageStatus$`Extinct?` & !LineageStatus$`Fused?`)
LineageExtinction[2,1] <- sum(!LineageStatus$`Extinct?` & LineageStatus$`Fused?`)
LineageExtinction[2,2] <- sum(!LineageStatus$`Extinct?` & !LineageStatus$`Fused?`)

# Run chi-square test
chisq.test(LineageExtinction)


