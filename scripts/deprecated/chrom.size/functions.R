GetPop <- function(N){
  pop <- matrix(0, N, 100)
  # this matrix represents the population
  # each row is an individual the values represent
  # the alleles present at a locus
  # 0 = 0/0
  # 1 = 0/1
  # 2 = 1/0
  # 3 = 1/1
  # columns 1:10 will represent autosome 1
  # columns 11:50 will represent autosome 2
  # columns 51:100 will represent the sex chromosomes
  # columns 10 and 50 will be SAL with 0 being female
  # benefit and 1 being male benefit.
  # column 100 will be the sex determining region with
  # 0 being the X and 1 being the Y


  # set half of individuals to be male
  pop[(N/2+1):N, 100] <- 1

  # set up SAL loci to have random alleles present
  pop[, c(10, 50)] <- sample(0:3, N*2, replace=T)
  return(pop)
}
Getdfe <- function(loci){
  dfe <- rgamma(loci, shape = .28, scale=113)
  # values less than -1 are not meaningful
  dfe[dfe>1] <- 1
  dfe <- 1-dfe
  dfe <- dfe-1
  # this vector holds our fitness values for
  # mutations ancestral alleles have fitness
  # of 1 while the mutations have fitness
  # distribution
}
GetFit <- function(pop, dfe){
  tpop <- pop
  pop[pop %in% c(1, 2)] <- .5
  pop[pop == 3] <- 1
  pop[tpop[, 100] == 1 & tpop[, 10] == 0, 10] <- 1
  pop[tpop[, 100] == 1 & tpop[, 10] == 3, 10] <- 0
  pop[tpop[, 100] == 1 & tpop[, 50] == 0, 50] <- 1
  pop[tpop[, 100] == 1 & tpop[, 50] == 3, 50] <- 0
  loci.fit <- 1 + pop %*% diag(dfe)
  g.fit <- apply(loci.fit[,-100], 1, prod)
  return(g.fit)
}
ConvertGenome <- function(genome){
  fgen <- matrix(NA, 2, 100)
  for(i in 1:100){
    switch(genome[i]+1,
           fgen[1:2, i] <- c(0, 0),
           fgen[1:2, i] <- c(0, 1),
           fgen[1:2, i] <- c(1, 0),
           fgen[1:2, i] <- c(1, 1))
  }
  return(fgen)
}

GetGametes <- function(pop, w, achiasmy = F){
  moms <- sample(1:500, 1000,
                 prob = w[1:500],
                 replace=T)
  dads <- sample(501:1000, 1000,
                 prob = w[501:1000],
                 replace = T)
  moms <- pop[moms, ]
  dads <- pop[dads, ]
  if(achiasmy == F){
    eggs <- matrix(0, 1000, 100)
    for(i in 1:1000){
      re.spots <- c(sample(2:9, 1),
                    sample(12:49, 1),
                    sample(51:99, 1))
      strand <- sample(1:2, 6, replace=T)
      cgen <- ConvertGenome(moms[i,])
      eggs[i, ] <- c(cgen[strand[1], 1:re.spots[1]],
                     cgen[strand[2], (re.spots[1] + 1):10],
                     cgen[strand[3], 11:re.spots[2]],
                     cgen[strand[4], (re.spots[2] + 1):50],
                     cgen[strand[5], 51:re.spots[3]],
                     cgen[strand[6], (re.spots[3] + 1):100])
    }
    sperm <- matrix(0, 1000, 100)
    for(i in 1:1000){
      re.spots <- c(sample(2:9, 1),
                    sample(12:49, 1),
                    sample(51:99, 1))
      strand <- sample(1:2, 6, replace=T)
      cgen <- ConvertGenome(dads[i,])
      sperm[i, ] <- c(cgen[strand[1], 1:re.spots[1]],
                      cgen[strand[2], (re.spots[1] + 1):10],
                      cgen[strand[3], 11:re.spots[2]],
                      cgen[strand[4], (re.spots[2] + 1):50],
                      cgen[strand[5], 51:re.spots[3]],
                      cgen[strand[6], (re.spots[3] + 1):100])
    }
  }
  gamete.list <- list(eggs, sperm)
  names(gamete.list) <- c("eggs","sperm")
  return(gamete.list)
}


GetNewGen <- function(gametes, pop){
  new.pop <- pop
  # get 500 X or Y bearing sperm
  ysperm <- sample(which(gametes$sperm[,100] == 1), 500, replace=T)
  xsperm <- sample(which(gametes$sperm[,100] == 0), 500, replace=T)
  # make females
  eggs.females <- sample(1:1000, 500)
  for(j in 1:500){
    # put two gametes together
    newgen <- rbind(gametes$eggs[eggs.females[j], ],
                    gametes$sperm[xsperm[j],])
    new.pop[j, colSums(newgen) == 0] <- 0
    new.pop[j, colSums(newgen) == 2] <- 3
    for(i in which(! colSums(newgen) %in% c(0,3))){
      if(all(newgen[,i] == c(0,1))){
        new.pop[j, i] <- 1
      }
      if(all(newgen[,i] == c(1,0))){
        new.pop[j,i] <- 2
      }
    }
  }
  # make males
  eggs.males <- sample(1:1000, 500)
  for(j in 1:500){
    # put two gametes together
    newgen <- rbind(gametes$eggs[eggs.males[j], ],
                    gametes$sperm[ysperm[j],])
    new.pop[(j+500), colSums(newgen) == 0] <- 0
    new.pop[(j+500), colSums(newgen) == 2] <- 3
    for(i in which(! colSums(newgen) %in% c(0,3))){
      if(all(newgen[,i] == c(0,1))){
        new.pop[(j+500), i] <- 1
      }
      if(all(newgen[,i] == c(1,0))){
        new.pop[(j+500),i] <- 2
      }
    }
  }
  return(new.pop)
}


