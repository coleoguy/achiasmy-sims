# generation time change
# commented out older getfreq ach and deleted line that messed up
# newer version
source("../scripts/functions.aneuploidy.R")
mut.rates <- c(1.45e-8,
               1.45e-10,
               1.45e-28)
N <- 1000
gens <- 500
chroms <- c("Y")
# sites of achiasmatic mutation
positions <- c(99)
pars <- get.pars(N = N,
                 s.ant = .3,
                 h.ant= 0.5,
                 site.sdr = 5,
                 site.sal = 95,
                 aneuploidy.rate <- 0,
                 parb = 10, #81 in aneuploidy/ old chrom, 10 in young
                 site.ach = positions,
                 mut.rate = mut.rates[3],
                 chromosome = chroms)
starting.pop <- get.initial.pop(pars)
pop <- starting.pop
segregating <- T
counter <- 1
time.table <- c()
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
  print(counter)
  counter <- counter + 1
}
plot(time.table)
