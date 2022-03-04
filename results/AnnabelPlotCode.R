setwd("C:/Users/knigh/Documents/GitHub/achiasmy-sims/results")

source("../scripts/functions.aneuploidy.R")

load("C:/Users/knigh/Documents/GitHub/achiasmy-sims/results/results.young.RData")

load("C:/Users/knigh/Documents/GitHub/achiasmy-sims/results/results.old.RData")

library("doParallel")
library("parallel")

########################## Young Y Chromosome Results ########################## 
resY <- all.results.young$Y
# General Results: The HIGHER the mutation rate, the GREATER the deleterious
# impact of Mueller's ratchet and thus the LESS likely the achiasmatic allele
# is to fix

# Each row is a different simulation
# Each column is a generation
# Each cell is proportion of achiasmatic allele at that generation

# Results for different mutation rates on Y chromosome
resY08 <- resY$`1.45e-08`
resY10 <- resY$`1.45e-10`
resY12 <- resY$`1.45e-12`

# Purpose: Explore location in genome/parameters which influence the conditions
# under which/location in which achiasmatic allele fixes

# Plot EVERY simulation for mutation rate 10^-8 on Y chr and find proportion of
# sims where it fixed after gen 500
library(viridisLite)
col.vec <- viridis(500)
plot(resY08[1,], type = "line", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Achiasmatic Allele",
     col = col.vec[1])

for(i in 2:nrow(resY08)){
  lines(y = resY08[i,], x = 1:500, col = col.vec[i])
}
plot(colMeans(resY08), pch = 16, cex = 0.5, ylim = c(0,1))
propFixedY08 <- sum(resY08[500,]==1)/500

# Plot for mutation rate 10^-10 on Y chr
plot(resY10[1,], type = "line", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Achiasmatic Allele",
     col = col.vec[1])

for(i in 2:nrow(resY10)){
  lines(y = resY10[i,], x = 1:500, col = col.vec[i])
}

# Plot for mutation rate 10^-10 on Y chr
plot(resY10[1,], type = "line", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Achiasmatic Allele",
     col = col.vec[1])

for(i in 2:nrow(resY10)){
  lines(y = resY10[i,], x = 1:500, col = col.vec[i])
}
plot(colMeans(resY10), pch = 16, cex = 0.5, ylim = c(0,1))
propFixedY10 <- sum(resY10[500,]==1)/500

# Plot for mutation rate 10^-12 on Y chr
plot(resY12[1,], type = "line", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Achiasmatic Allele",
     col = col.vec[1])

for(i in 2:nrow(resY10)){
  lines(y = resY10[i,], x = 1:500, col = col.vec[i])
}
plot(colMeans(resY12), pch = 16, cex = 0.5, ylim = c(0,1))
propFixedY12 <- sum(resY12[500,]==1)/500

############################ Young X Chromosome ################################
resX <- all.results.young$X
resX12 <- resX$`1.45e-12`
# Plot for mutation rate 10^-12 on X chr
plot(resX12[1,], type = "line", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Achiasmatic Allele",
     col = col.vec[1])

for(i in 2:nrow(resX12)){
  lines(y = resX12[i,], x = 1:500, col = col.vec[i])
}
plot(colMeans(resX12), pch = 16, cex = 0.5, ylim = c(0,1))
propFixedX12 <- sum(resX12[500,]==1)/500

# Visualize proportion fixed over time for all three chromosome classes at the
# highest mutation rate
fixationX <- c()
fixationY <- c()
fixationA <- c()
for(i in 1:500){
  fixationX[i] <- sum(resX12[,i]==1)/500
  fixationY[i] <- sum(resY12[,i]==1)/500
  fixationA[i] <- sum(all.results.young$A$`1.45e-12`[,i]==1)/500
}

# Achiasmy is most likely to fix on Y chromosome becasue fixation of achiasmy 
# allele on Y chromosome ensures that achiasmy allele is only ever inherited by
# males.
# Fixation of achiasmy allele on Y ensures that achiasmy allele is always linked
# with the male-beneift allele of the sexually-antagonistic locus
# Achiasmy mutation does not ever fix on autosome because autosomes segregate
# independently of sex chromosomes
plot(fixationX, type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Achiasmatic Allele",
     col = col.vec[1])

lines(y = fixationY, x = 1:500, col = "blue")
lines(y = fixationA, x = 1:500, col = "red")

################################### Adjust DFE #################################
# We want the mean DFE to have fewer extremely deleterious alleles
dfe <- rgamma(5000, shape = .28, scale=113)
# values less than -1 are not meaningful
dfe[dfe>1] <- 1
dfe <- 1-dfe
dfe <- dfe-1
hist(dfe[dfe<1], breaks = 50)

# Overwrite get.pars() to have DFE with fewer extremely deleterious mutations
get.pars <- function(N=NULL, 
                     s.ant=NULL,    
                     site.sdr=NULL, 
                     site.sal=NULL, 
                     site.ach=NULL, 
                     mut.rate=NULL, 
                     h.ant=NULL,    
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
  dfe <- dfe-0.8
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


no_cores <- detectCores(logical = TRUE)
cl <- makeCluster(no_cores-1)
registerDoParallel(cl)
mut.rates <- c(1.45e-8,
               1.45e-9,
               1.45e-10,
               1.45e-11,
               1.45e-12)
s.ant = .5
h.ant= 0.5
N <- 1000
gens <- 500
iter <- 1000
chroms <- c("X","Y","A")
positions <- c(99,99,100)
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
save(all.results.young,file = "DFEedits_all.results.young")

