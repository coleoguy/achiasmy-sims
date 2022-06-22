## code from final.loop.r use to do trial plots of y chromosome alone

source("functions.aneuploidy.R")
library(doParallel)
library(parallel)
no_cores <- detectCores(logical = TRUE)
cl <- makeCluster(no_cores-1)
registerDoParallel(cl)

mut.rates <- c(1.45e-8,
               1.45e-9,
               1.45e-10,
               1.45e-11,
               1.45e-12)

N <- 1000
gens <- 500
iter <- 1000
chroms <- c("Y")
# sites of achiasmatic mutation
positions <- c(99)


#young chromosome
all.results.young <- vector(length=1, mode="list")
names(all.results.young) <- chroms
# k cycles through chromosomes and positions of achiasmatic mutation
for(k in 1:1){
  print(paste("running", chroms[k], "chromosome model"))
  # set up container for results
  results.time <- vector(length=length(mut.rates), mode="list")
  names(results.time) <- mut.rates
  
  # loops through and tests each mutation rate
  for(j in 1:length(mut.rates)){
    print(paste("runnning mutation rate:", j, "of", length(mut.rates)))
    pars <- get.pars(N = N,
                     s.ant = .5,
                     h.ant= 0.5,
                     site.sdr = 5,
                     site.sal = 95,
                     aneuploidy.rate <- 0,
                     parb = 10, #81 in aneuploidy/ old chrom, 10 in young
                     site.ach = positions[k],
                     mut.rate = mut.rates[j],
                     chromosome = chroms[k]) #1.45e-8
    starting.pop <- get.initial.pop(pars)
    x <- foreach(k=1:iter, .combine = "rbind") %dopar%{
      time.table <- c()
      pop <- starting.pop
      segregating <- T
      counter <- 1
      while(segregating){
        pop <- generation(pop, pars)
        cur.freq <- get.freq.ach(pop, pars)
        time.table[counter] <- cur.freq
        if(counter == gens) segregating <- F
        if(cur.freq == 1){
          time.table[length(time.table):gens] <- 1
          segregating <- F
        }
        if(counter %% 10 == 0) pop <- achMut(pop, pars)
        counter <- counter + 1
      }
      time.table
    }
    results.time[[j]] <- x
  }
  all.results.young[[k]] <- results.time
}
save(all.results.young,file = "all.results.young")

#young chromosome, parb=0
all.results.parb0 <- vector(length=1, mode="list")
names(all.results.parb0) <- chroms
# k cycles through chromosomes and positions of achiasmatic mutation
for(k in 1:1){
  print(paste("running", chroms[k], "chromosome model"))
  # set up container for results
  results.time <- vector(length=length(mut.rates), mode="list")
  names(results.time) <- mut.rates
  
  # loops through and tests each mutation rate
  for(j in 1:length(mut.rates)){
    print(paste("runnning mutation rate:", j, "of", length(mut.rates)))
    pars <- get.pars(N = N,
                     s.ant = .5,
                     h.ant= 0.5,
                     site.sdr = 5,
                     site.sal = 95,
                     aneuploidy.rate <- 0,
                     parb = 0, #81 in aneuploidy/ old chrom, 10 in young
                     site.ach = positions[k],
                     mut.rate = mut.rates[j],
                     chromosome = chroms[k]) #1.45e-8
    starting.pop <- get.initial.pop(pars)
    x <- foreach(k=1:iter, .combine = "rbind") %dopar%{
      time.table <- c()
      pop <- starting.pop
      segregating <- T
      counter <- 1
      while(segregating){
        pop <- generation(pop, pars)
        cur.freq <- get.freq.ach(pop, pars)
        time.table[counter] <- cur.freq
        if(counter == gens) segregating <- F
        if(cur.freq == 1){
          time.table[length(time.table):gens] <- 1
          segregating <- F
        }
        if(counter %% 10 == 0) pop <- achMut(pop, pars)
        counter <- counter + 1
      }
      time.table
    }
    results.time[[j]] <- x
  }
  all.results.parb0[[k]] <- results.time
}
save(all.results.parb0,file = "all.results.parb0")


#young chromosome, sal=10
all.results.sal10 <- vector(length=1, mode="list")
names(all.results.sal10) <- chroms
# k cycles through chromosomes and positions of achiasmatic mutation
for(k in 1:1){
  print(paste("running", chroms[k], "chromosome model"))
  # set up container for results
  results.time <- vector(length=length(mut.rates), mode="list")
  names(results.time) <- mut.rates
  
  # loops through and tests each mutation rate
  for(j in 1:length(mut.rates)){
    print(paste("runnning mutation rate:", j, "of", length(mut.rates)))
    pars <- get.pars(N = N,
                     s.ant = .5,
                     h.ant= 0.5,
                     site.sdr = 5,
                     site.sal = 10,
                     aneuploidy.rate <- 0,
                     parb = 10, #81 in aneuploidy/ old chrom, 10 in young
                     site.ach = positions[k],
                     mut.rate = mut.rates[j],
                     chromosome = chroms[k]) #1.45e-8
    starting.pop <- get.initial.pop(pars)
    x <- foreach(k=1:iter, .combine = "rbind") %dopar%{
      time.table <- c()
      pop <- starting.pop
      segregating <- T
      counter <- 1
      while(segregating){
        pop <- generation(pop, pars)
        cur.freq <- get.freq.ach(pop, pars)
        time.table[counter] <- cur.freq
        if(counter == gens) segregating <- F
        if(cur.freq == 1){
          time.table[length(time.table):gens] <- 1
          segregating <- F
        }
        if(counter %% 10 == 0) pop <- achMut(pop, pars)
        counter <- counter + 1
      }
      time.table
    }
    results.time[[j]] <- x
  }
  all.results.sal10[[k]] <- results.time
}
save(all.results.sal10,file = "all.results.sal10")


#young chromosome, s.ant0.8
all.results.sant0.8 <- vector(length=1, mode="list")
names(all.results.sant0.8) <- chroms
# k cycles through chromosomes and positions of achiasmatic mutation
for(k in 1:1){
  print(paste("running", chroms[k], "chromosome model"))
  # set up container for results
  results.time <- vector(length=length(mut.rates), mode="list")
  names(results.time) <- mut.rates
  
  # loops through and tests each mutation rate
  for(j in 1:length(mut.rates)){
    print(paste("runnning mutation rate:", j, "of", length(mut.rates)))
    pars <- get.pars(N = N,
                     s.ant = .8,
                     h.ant= 0.5,
                     site.sdr = 5,
                     site.sal = 10,
                     aneuploidy.rate <- 0,
                     parb = 10, #81 in aneuploidy/ old chrom, 10 in young
                     site.ach = positions[k],
                     mut.rate = mut.rates[j],
                     chromosome = chroms[k]) #1.45e-8
    starting.pop <- get.initial.pop(pars)
    x <- foreach(k=1:iter, .combine = "rbind") %dopar%{
      time.table <- c()
      pop <- starting.pop
      segregating <- T
      counter <- 1
      while(segregating){
        pop <- generation(pop, pars)
        cur.freq <- get.freq.ach(pop, pars)
        time.table[counter] <- cur.freq
        if(counter == gens) segregating <- F
        if(cur.freq == 1){
          time.table[length(time.table):gens] <- 1
          segregating <- F
        }
        if(counter %% 10 == 0) pop <- achMut(pop, pars)
        counter <- counter + 1
      }
      time.table
    }
    results.time[[j]] <- x
  }
  all.results.sant0.8[[k]] <- results.time
}
save(all.results.sant0.8,file = "all.results.sant0.8")




