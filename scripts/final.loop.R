source("functions.aneuploidy.R")
library(doParallel)
library(parallel)
no_cores <- detectCores(logical = TRUE)
cl <- makeCluster(no_cores-1)
registerDoParallel(cl)

# mutation rates
mut.rates <- c(1.45e-8,
               1.45e-9,
               1.45e-10,
               1.45e-11,
               1.45e-12)
# selection coef of sex ant mutation
# run with this dropped at .5 basically always fixed with above
# mutation rates.
s.ant = .5

# dominance factor of the male benefit allele
h.ant= 0.5

# population size
N <- 1000

# generations
gens <- 500

#iterations
iter <- 1000

# chromosomes to try achiasmatic mutations on
chroms <- c("X","Y","A")

# sites of achiasmatic mutation
positions <- c(99,99,100)


#young chromosome
all.results.young <- vector(length=5, mode="list")
names(all.results.young) <- chroms
# k cycles through chromosomes and positions of achiasmatic mutation
for(k in 1:3){
  print(paste("running", chroms[k], "chromosome model"))
  # set up container for results
  results.time <- vector(length=length(mut.rates), mode="list")
  names(results.time) <- mut.rates

  # loops through and tests each mutation rate
  for(j in 1:length(mut.rates)){
    print(paste("runnning mutation rate:", j, "of", length(mut.rates)))
    pars <- get.pars(N = N,
                     s.ant = s.ant,
                     h.ant= h.ant,
                     site.sdr = 5,
                     site.sal = 95,
                     # line to change for aneuploidy = prob of 
                     # a aneuploidy event in males
                     aneuploidy.rate <- 0,
                     #81 in aneuploidy/ old chrom, 10 in young
                     parb = 10, 
                     site.ach = positions[k],
                     mut.rate = mut.rates[j],
                     chromosome = chroms[k]) #1.45e-8
    starting.pop <- get.initial.pop(pars)
    x <- foreach(m=1:iter, .combine = "rbind") %dopar%{
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

