source("functions.R")
# get initial population
pop <- GetPop(1000)

# introduce sexual antagonistic variation
pop[ ,c(10,50)] <- sample(0:3, 2000,
                          prob=c(.6,.1,.1,.1),
                          replace=T)


# get mutation dfe
dfe <- Getdfe(100)
dfe[c(10,50)] <- -0.5

# sex specific fitness
# example of fitness function for s=.2, h=.5
#    Female         Male
# 0 00 1        1.0   1-s   0.8
# 1 01 1-hs     0.9   1-hs  0.9
# 2 10 1-hs     0.9   1-hs  0.9
# 3 11 1-s      0.8   1     1.0


# assess fitness
SALmat10 <- SALmat50 <- matrix(NA, 1000, 2000)
for(i in 1:2000){
  print(i)
  w <- GetFit(pop, dfe)
  gametes <- GetGametes(pop, w, achiasmy = F)
  pop <- GetNewGen(gametes, pop)
  SALmat10[, i] <- pop[,10]
  SALmat50[, i] <- pop[,50]
}
SALmat10[SALmat10 == 2] <- 1
SALmat10[SALmat10 == 3] <- 2
SALmat50[SALmat50 == 2] <- 1
SALmat50[SALmat50 == 3] <- 2

plot(colSums(SALmat10)/2000, type="l",
     ylab="Freq of male benefit allele")
lines(colSums(SALmat50)/2000, col="red")
