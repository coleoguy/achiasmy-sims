#use to make figure with 5 levels of mutation rates 
#(-8, -9, -10, -11, -12) on X, Y, and autosome as the
#frequency of meiosis mutation for given iteration 

library(viridis)
cols <- viridis(100) #number of iterations

par(pty="s")
par(mfcol=c(1,3)) # correspond to number of mutations (5 total)

#in final result plots the
#x axis= generation number (gens=500)
#y axis= probability of fixation 
#each line= one iteration (iter=1000)

# young model
# Y chromosome
# all.results corresponds to which model was simulated
for (j in 1:3){ #3 mut rates
  plot(x=1:100, y=all.results.young[[1]][[j]][1,], col=cols[1], type = "l", 
       xlab="Generation", ylab="Frequency of Fixation", xlim=c(0,150), ylim=c(0,1),lwd=.05)
  #i should go from 2 to the number of iterations
  for(i in 2:100){
    lines(x=1:100, y=all.results.young[[1]][[j]][i,], col=cols[i],lwd=.05)
  }
  fixed <- sum(all.results.young[[1]][[j]][,100]==1)
  lost <- sum(all.results.young[[1]][[j]][,100]==0)
  text(x=125, y=.95, fixed, pos=4)
  text(x=125, y=.05, lost, pos=4)
  #title(main = "Young Y chromosome")
}




#aneuploidy=0.5
# Y chromosome
for (j in 1:3){ #5 mut rates
  plot(x=1:100, y=all.results.aneu0.5[[1]][[j]][1,], col=cols[1], type = "l", 
       xlab="Generation", ylab="Frequency of Fixation", xlim=c(0,150), ylim=c(0,1),lwd=.05)
  #i should go from 2 to the number of iterations
  for(i in 2:100){
    lines(x=1:100, y=all.results.aneu0.5[[1]][[j]][i,], col=cols[i],lwd=.05)
  }
  fixed <- sum(all.results.aneu0.5[[1]][[j]][,100]==1)
  lost <- sum(all.results.aneu0.5[[1]][[j]][,100]==0)
  text(x=125, y=.95, fixed, pos=4)
  text(x=125, y=.05, lost, pos=4)
  #title(main = "Y chromosome, aneuploidy rate=0.5")
}


