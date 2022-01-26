load("../results/Ychromosome.RData")
res.Y <- results.time
load("../results/Xchromosome.RData")
res.X <- results.time
load("../results/autosome.RData")
res.A <- results.time
rm(results.time)
fixation <- as.data.frame(matrix(,3,5))
row.names(fixation) <- c("y","x","a")
colnames(fixation) <- c(1.45e-8, 1.45e-9, 1.45e-10, 1.45e-11, 1.45e-12)
for(j in 1:5){
  fixation[1, j] <- sum(res.Y[[j]][,500]==1)
  fixation[2, j] <- sum(res.X[[j]][,500]==1)
  fixation[3, j] <- sum(res.A[[j]][,500]==1)
}


plot(x=1:5, y=fixation[1,],type="b",pch=16,cex=.75)
lines(x=1:5, y=fixation[2,],type="b",pch=16,cex=.75)
lines(x=1:5, y=fixation[3,],type="b",pch=16,cex=.75)

unlist(fixation)
