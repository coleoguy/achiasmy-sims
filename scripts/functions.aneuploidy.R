
# this function just gives us all our model pars
get.pars <- function(N=NULL, # the size of the populations
                     s.ant=NULL,    # selection on sexual antagonism
                     site.sdr=NULL, # site of the sex determining gene
                     site.sal=NULL, # site of the sexually antagonisitic gene
                     site.ach=NULL, # site of the mutation that causes achiasmatic meiosis
                     mut.rate=NULL, # mutation rate for all other genes
                     h.ant=NULL,    # genetic architecture of antagonistic locus
                     chromosome=NULL,
                     aneuploidy.rate=NULL,
                     parb = NULL){ #first value in par boundary

  # set chromosome to have mutation occur on
  if(is.null(chromosome)) chromosome <- "Y"
  # number of sites on sex chromosome
  sites <- 100

  # population size
  if(is.null(N)) N <- 500

  # architecture
  if(is.null(site.sdr)) site.sdr <- 5  # the SDR will have alleles of 0 and 1 representing X and Y resp.
  if(is.null(site.sal)) site.sal <- 10 # 0 good for females 1 good for males, site 50 also used
  if(is.null(site.ach)) site.ach <- 100 #site 60 used also

  if(is.null(parb)) parb <- 10 #region of PAR from parb-100

  # crossover the X chromosomes is about 197 centimorgans long
  # as such we should expect an average of two crossovers in our
  # model this would mean that the cM distance between loci should
  # be .02 this leads to an average of 1.94 crossovers per chromosome

  site.dist <- .02

  # we want the chance to unlink the achiasmatic mutation
  # from the sex chromosome so we will have a seperate
  # linkage factor that can be between 0 or 0.5 with 0.5
  # representing autosomal and <0.5 sex linked
  dist.ach <- .5

  #### Explanation of mutation rate choice ####
  # to maintain inference power to humans and other
  # real genomes we want a mutation rate that will
  # provide the correct number of mutations for an
  # x chromosome that is "real" to do this lets pattern
  # things after the human X chromosome. which has 1588410
  # bp of exonic DNA.
  # CITE: Sakharkar, M.K., Chow, V.T. and Kangueane, P.,
  # 2004. Distributions of exons and introns in the human genome.
  # In silico biology, 4(4), pp.387-393.
  # We expect the germline exome mutation
  # rate to be on the order of 1.45 X 10-8  per bp per year
  #CITE: Scally, A., The mutation rate in human evolution and demographic inference. 2016 Current opinion
  # in genetics & development, 41, pp.36-43.
  # CITE: Narasimhan VM, Rahbari R, Scally A, Wuster A, Mason D, Xue Y, Wright J, Trembath RC,
  # Maher ER, van Heel DA, et al.: A direct multi-generational estimate of
  # the human mutation rate from autozygous segments seen in thousands of
  # parentally related individuals. bioRxiv 2016.
  # if we take these numbers and then apply a generation time of 25 years
  # we can us a binomial distribution to find the number of mutations that
  # individuals will contain:
  # the prob argument is created by multiplying the annual mutation rate
  # by a generation time of 25y plus times two since we are dealing with diploids
  #### Explanation of mutation rate choice ####
  # normal and realistic 1.45e-08
  if(is.null(mut.rate)) mut.rate <- 1.45e-8
  mut.rate <- mut.rate * 2 #generation time and diploid accounted for

  # distribution of mutational effects to create something realistic
  # we will use a compound normal distribution to describe the value of
  # s where fitness is equal to 1+s

  #previously used dfe from functions.R
  # dfe <- rexp(5000,rate=3.5)-1
  # values less than -1 are not meaninful
  #dfe <- dfe[dfe<=1]


  #dfe <- rexp(5000,rate=3.5)-1
  dfe <- rgamma(5000, shape = .28, scale=113)
  # values less than -1 are not meaningful
  dfe[dfe>1] <- 1
  dfe <- 1-dfe
  dfe <- dfe-1
  #hist(dfe)
  #sum(dfe <= -.99)/length(dfe)
  #sum(dfe >= -.01)/length(dfe)


  # strength of selection on sexually antagonistic locus
  if(is.null(s.ant)) s.ant <- 0.5

  # genetic architecture of the sexually antagonistic locus
  if(is.null(h.ant)) h.ant <- 0.5

  # rate of aneuploidy in males
  if(is.null(aneuploidy.rate)) aneuploidy.rate <- 0.1

  pars <- list(N,
               sites,
               site.dist,
               site.sdr,
               site.sal,
               site.ach,
               dist.ach,
               h.ant,
               s.ant,
               dfe,
               mut.rate,
               chromosome,
               aneuploidy.rate,
               parb)

  names(pars) <- c("N",
                   "sites",
                   "site.dist",
                   "site.sdr",
                   "site.sal",
                   "site.ach",
                   "dist.ach",
                   "h.ant",
                   "s.ant",
                   "dfe",
                   "mut.rate",
                   "chromosome",
                   "aneuploidy.rate",
                   "parb")
  return(pars)
}

# this function initializes the population
get.initial.pop <- function(pars){
  pop.xy <- list()
  for(i in 1:pars$N){
    pop.xy[[i]] <- matrix(1, 2, pars$sites)

    # the SDR will have alleles of 0 and 1 representing X and Y resp.
    pop.xy[[i]][, pars$site.sdr] <- sample(c(0,0,0,1), size=2)

    # the SAL will have values of 0 and 1 representing female and male
    # benefit alleles respectively we will start off with all allele 0
    # and allow for the evolution of allele 1.
    pop.xy[[i]][, pars$site.sal] <- sample(0:1, 2, replace=T)
  }
  return(pop.xy)
}

# this function calculates the fitness of a genome
### TODO WORK OUT FITNESS FUNCTION THEN COME BACK AND CODE THIS
get.fit <- function(genome, pars){
  gen.fit <- (1:pars$sites)[-c(pars$site.sdr, pars$site.sal, pars$site.ach)]
  #general fitness includes all sites except sdr, sal, and ach
  gen.fit <- prod(colMeans(genome[,gen.fit]))
  m.f <- 1 #meiotic fitness
  # males
  if(sum(genome[,pars$site.sdr]) == 1){
    switch(as.character(sum(genome[,pars$site.sal])),
           "0" = s.f <- 1,
           "1" = s.f <- 1 + pars$h.ant * pars$s.ant,
           "2" = s.f <- 1 + pars$s.ant)

    switch(as.character(sum(genome[,pars$site.ach])),
           "0" = m.f <- 1 - pars$aneuploidy.rate,
           "1" = m.f <- 1,
           "2" = m.f <- 1)
  }
  # female
  if(sum(genome[,pars$site.sdr]) == 0){
    switch(as.character(sum(genome[,pars$site.sal])),
           "0" = s.f <- 1,
           "1" = s.f <- 1 / (1 + pars$h.ant * pars$s.ant),
           "2" = s.f <- 1 / (1 + pars$s.ant))
  }
  t.fit <- gen.fit * s.f * m.f #total fitness
  t.fit[t.fit < 0] <- 0
  return(t.fit)
}

# this function gives a vector of sexes
get.dams <- function(pop.xy){
  females <- c()
  for(i in 1:pars$N){
    switch(as.character(sum(pop.xy[[i]][,pars$site.sdr])),
           "0" = females[i] <- T,
           "1" = females[i] <- F)
  }
  return(females)
}

# make parent table
get.parents <- function(pop.xy, pars){
  parents <- matrix(, pars$N, 2)
  colnames(parents) <- c("sire","dam")
  fitness <- unlist(lapply(pop.xy, get.fit, pars))
  dams <- get.dams(pop.xy)
  fitness[!dams] <- fitness[!dams]/max(fitness[!dams])
  fitness[dams] <- fitness[dams]/max(fitness[dams])
  parents[, 1] <- sample(x = (1:pars$N)[!dams],
                         prob = fitness[!dams],
                         size = pars$N,
                         replace = T)
  parents[, 2] <- sample(x = (1:pars$N)[dams],
                         prob = fitness[dams],
                         size = pars$N,
                         replace = T)
  return(parents)
}

# this lays down mutations, same code as "young chromosome"
radiate <- function(pop.xy, pars){
  # 1587926 bp of exonic DNA in human X Chromosome
  hit.indiv <- rbinom(n=pars$N,size=1587926, prob=pars$mut.rate)
  for(i in 1:pars$N){
    if(hit.indiv[i] != 0){
      hit.sites <- (1:pars$sites)[-c(pars$site.sdr,
                                   pars$site.sal,
                                   pars$site.ach)]
      hit.sites <- sample(x=hit.sites, size = hit.indiv[i])
      for(j in 1:hit.indiv[i]){
        pop.xy[[i]][sample(x=0:1, 1), hit.sites[j]] <- 1+sample(pars$dfe, 1)
      }
    }
  }
  # TODO think about whether we want a more nuanced way of deciding if and
  # when SAL mutations occur this introduces mutations among the two alleles at the SAL
  if(sample(1:99, 1) >= pars$parb){
    pop.xy[[sample(1:pars$N, 1)]][sample(c(1,2), 1), pars$site.sal] <- sample(c(0,1), 1)
  }
  return(pop.xy)
}
# two alleles for site.ach, 0 for achiasmatic, 1 for chiasmatic meiosis
get.gamete <- function(genome, pars){
  # this will hold strand information
  strand <- c()
  # case female meiosis
  if(sum(genome[,pars$site.sdr])==0){
    # find crossover spots
    #TODO
    x <- which(sample(x=0:1,prob=c(.98,.02),replace=T,size=99)==1)+1 #Julia: Should this be 1:2
    #instead of 0:1
    # if crossover is happening on this gamete
    if(length(x) > 0){
      # loop through each crossover spot
      for(i in 1:length(x)){
        # case of a single crossover event
        if(i == 1) strand <- rep(sample(x=1:2,1), times = x[i])
        # case of more than one crossover event
        if(i > 1){
          z <- x[i] - x[(i-1)]
          # this switches from the current
          # chromosome to the alternate chromosome
          switch(strand[length(strand)],
                 "1" = strand <- c(strand, rep(2, times=z)),
                 "2" = strand <- c(strand, rep(1, times=z)))
        }
      }
      # this finishes the end of the chromosome
      if(length(strand) < pars$sites){
        z <- pars$sites - length(strand)
        switch(strand[length(strand)],
               "1" = strand <- c(strand, rep(2, times=z)),
               "2" = strand <- c(strand, rep(1, times=z)))
      }
    }
    # if no recombination just picking a chomosome to pass on
    if(length(x) == 0) strand <- rep(sample(1:2, 1), times=pars$sites)
    # shouldnt have been doing achiasmy in homogametic sex
    # this is for the case of achiasmatic meiosis
    # if(sum(genome[, pars$site.ach]) < 2){
    #   strand <- rep(sample(1:2, 1), times = pars$sites)
    #   # for the case of autosomal achiasmatic mutation
    #   if(pars$site.ach == 100){
    #     if(pars$dist.ach == .5){
    #       strand[length(strand)] <- sample(1:2, 1)
    #     }
    #   }
    # }
  }

  # case for males
  if(sum(genome[,pars$site.sdr])==1){
    # greater than zero we have chiasmatic meiosis
    if(sum(genome[, pars$site.ach]) == 2){
      # defines the pairing region
      x <- sample(pars$parb:99, size=1) #region able to have crossover event, smaller than
      #in functions.R so aneuploidy event is more likely to occur
      start.strand <- sample(1:2,1)
      strand <- rep(start.strand, times = x)
      z <- pars$sites - length(strand)
      # this is the same as the switch statement above
      if(start.strand==1) strand[(x+1):pars$sites] <- 2
      if(start.strand==2) strand[(x+1):pars$sites] <- 1
    }
    # this is for the case of achiasmatic meiosis
    if(sum(genome[, pars$site.ach]) < 2){
      strand <- rep(sample(1:2, 1), times = pars$sites)
      # for the case of autosomal achiasmatic mutation
      if(pars$site.ach == 100){
        if(pars$dist.ach == .5){
          strand[length(strand)] <- sample(1:2, 1)
        }
      }
    }
  }
  # we are now at the point that we have strand information
  # for our gamete
  l.genome <- c(genome[1, ], genome[2, ])
  gamete <- l.genome[(strand - 1) * pars$sites + 1:pars$sites]
  return(gamete)
}

# this makes offspring
get.offspring <- function(pop.xy, parents, pars){
  pop.xy2 <- pop.xy
  for(i in 1:length(pop.xy2)){
    for(j in 1:ncol(parents)){
      pop.xy2[[i]][1, ] <- get.gamete(pop.xy[[parents[i,1]]], pars)
      pop.xy2[[i]][2, ] <- get.gamete(pop.xy[[parents[i,2]]], pars)
    }
  }
  return(pop.xy2)
}

# this runs one generation
generation <- function(pop.xy, pars){
  pop.xy.big <- do.call(rbind, pop.xy)
  check.cols <- 1:pars$sites[-c(pars$site.sdr,
                                pars$site.sal,
                                pars$site.ach)]
  for(i in check.cols){
    if(length(unique(pop.xy.big[, i])) == 1){
      pop.xy.big[, i] <- 1
    }
  }
  pop.xy <- lapply(split(pop.xy.big, rep(1:(nrow(pop.xy.big)/2), each=2)), matrix, ncol=pars$sites)
  pop.xy <- radiate(pop.xy, pars)
  parents <- get.parents(pop.xy, pars)
  pop.xy <- get.offspring(pop.xy, parents, pars)
  return(pop.xy)
}

#frequency of males vs females
get.freqs <- function(pop.xy, pars){
  x.alleles <- y.alleles <- c()
  for(i in 1:pars$N){
    genome <- pop.xy[[i]]
    xs <- genome[, pars$site.sdr] == 0
    x.alleles <- c(x.alleles, genome[xs, pars$site.sal])
    ys <- genome[, pars$site.sdr] == 1
    y.alleles <- c(y.alleles, genome[ys, pars$site.sal])
  }
  y <- sum(y.alleles)/length(y.alleles)
  x <- sum(x.alleles)/length(x.alleles)
  report <- data.frame(x,y)
  colnames(report) <- c("X","Y")
  row.names(report) <- "frequency of male benefit"
  return(report)
}

# calculates heterozygosity probability that a
# site is heterozygous in population
get.het <- function(pop.xy, pars){
  pop.xy.big <- do.call(rbind, pop.xy)
  keep <- (1:pars$sites)[-c(pars$site.sdr,
             pars$site.sal,
             pars$site.ach)]
  x <- 0
  for(i in keep){
    if(length(unique(pop.xy.big[,i])) > 1) x <- x + 1
  }
  return(x/97)
}



  # get frequency of achiasmatic meiosis mutation
  get.freq.ach <- function(pop.xy, pars){
    males <- !get.dams(pop.xy)
    ach.alleles <- 0
    for(i in 1:pars$N){
      genome <- pop.xy[[i]]
      if(pars$chromosome=="Y"){
        if(males[i]){
          ach.alleles <- ach.alleles + sum(genome[1, pars$site.ach] == 0)
        }
      }
      if(pars$chromosome=="X"){
        if(!males[i]){
          ach.alleles <- ach.alleles + sum(genome[1:2, pars$site.ach] == 0)
        }
        if(males[i]){
          ach.alleles <- ach.alleles + sum(genome[2, pars$site.ach] == 0)
        }
      }
      if(pars$chromosome=="A"){
        if(!males[i]){
          ach.alleles <- ach.alleles + sum(genome[1:2, pars$site.ach] == 0)
        }
        if(males[i]){
          ach.alleles <- ach.alleles + sum(genome[1:2, pars$site.ach] == 0)
        }
      }
    }
    x.count <- sum(get.dams(pop.xy))*2 + sum(!get.dams(pop.xy))
    y.count <- sum(!get.dams(pop.xy))
    if(pars$chromosome=="Y") ach.freq <- ach.alleles/y.count
    if(pars$chromosome=="X") ach.freq <- ach.alleles/x.count
    if(pars$chromosome=="A") ach.freq <- ach.alleles/(pars$N*2)
    return(ach.freq)
  }


########
# # get frequency of achiasmatic meiosis mutation don't use
# get.freq.ach <- function(pop.xy, pars){
#   males <- !get.dams(pop.xy)
#   ach.alleles <- 0
#   for(i in 1:pars$N){
#     genome <- pop.xy[[i]]
#     if(pars$chromosome=="Y"){
#       if(males[i]){
#         ach.alleles <- ach.alleles + sum(genome[1, pars$site.ach] == 0)
#       }
#     }
#     if(pars$chromosome=="X"){
#       if(!males[i]){
#         ach.alleles <- ach.alleles + sum(genome[1:2, pars$site.ach] == 0)
#       }
#       if(males[i]){
#         ach.alleles <- ach.alleles + sum(genome[2, pars$site.ach] == 0)
#       }
#     }
#     if(pars$chromosome=="A"){
#       if(!males[i]){
#         ach.alleles <- ach.alleles + sum(genome[1:2, pars$site.ach] == 0)
#       }
#       if(males[i]){
#         ach.alleles <- ach.alleles + sum(genome[1:2, pars$site.ach] == 0)
#       }
#     }
#   }
#   x.count <- sum(get.dams(pop.xy))*2 + sum(!get.dams(pop.xy))
#   y.count <- sum(!get.dams(pop.xy))
#   if(pars$chromosome=="Y") ach.freq <- ach.alleles/y.count
#   if(pars$chromosome=="X") ach.freq <- ach.alleles/x.count
#   if(pars$chromosome=="A") ach.freq <- ach.alleles/(pars$N*2)
#   return(ach.freq)
# }
#
# get frequency of achiasmatic meiosis mutation
#######

achMut <- function(pop.xy, pars){
  males <- !get.dams(pop.xy)
  hit <- sample(which(males), 1)
  if(pars$chromosome == "Y") pop.xy[[hit]][1 ,pars$site.ach] <- 0
  if(pars$chromosome == "X") pop.xy[[hit]][2 ,pars$site.ach] <- 0
  if(pars$chromosome == "A") pop.xy[[hit]][sample(1:2, 1) ,pars$site.ach] <- 0
  return(pop.xy)
}




