#use to make figure with 5 levels of mutation rates 
#(-8, -9, -10, -11, -12) on X, Y, and autosome as the
#frequency of meiosis mutation for given iteration 

library(viridis)
cols <- viridis(1000) #number of iterations

par(pty="s")
par(mfcol=c(1,5)) # correspond to number of mutations (5 total)

#in final result plots the
#x axis= generation number (gens=500)
#y axis= probability of fixation 
#each line= one iteration (iter=1000)

# young model
# Y chromosome
# all.results corresponds to which model was simulated
for (j in 1:5){ #5 mut rates
  plot(x=1:500, y=all.results.young[[2]][[j]][1,], col=cols[1], type = "l", 
       xlab="Generation", ylab="Frequency of Fixation", xlim=c(0,600), ylim=c(0,1),lwd=.05)
  #i should go from 2 to the number of iterations
  for(i in 2:1000){
    lines(x=1:500, y=all.results.young[[2]][[j]][i,], col=cols[i],lwd=.05)
  }
  fixed <- sum(all.results.young[[2]][[j]][,500]==1)
  lost <- sum(all.results.young[[2]][[j]][,500]==0)
  text(x=550, y=.95, fixed, pos=4)
  text(x=550, y=.05, lost, pos=4)
  title(main = "Young Y chromosome")
}


# X chromosome
for (j in 1:5){ #5 mut rates
  plot(x=1:500, y=all.results.young[[1]][[j]][1,], col=cols[1], type = "l", 
       xlab="Generation", ylab="Frequency of Fixation", xlim=c(0,600), ylim=c(0,1),lwd=.05)
  #i should go from 2 to the number of iterations
  for(i in 2:1000){
    lines(x=1:500, y=all.results.young[[1]][[j]][i,], col=cols[i],lwd=.05)
  }
  fixed <- sum(all.results.young[[1]][[j]][,500]==1)
  lost <- sum(all.results.young[[1]][[j]][,500]==0)
  text(x=550, y=.95, fixed, pos=4)
  text(x=550, y=.05, lost, pos=4)
  title(main = "Young X chromosome")
}


# autosome
for (j in 1:5){ #5 mut rates
  plot(x=1:500, y=all.results.young[[3]][[j]][1,], col=cols[1], type = "l", 
       xlab="Generation", ylab="Frequency of Fixation", xlim=c(0,600), ylim=c(0,1),lwd=.05)
  #i should go from 2 to the number of iterations
  for(i in 2:1000){
    lines(x=1:500, y=all.results.young[[3]][[j]][i,], col=cols[i],lwd=.05)
  }
  fixed <- sum(all.results.young[[3]][[j]][,500]==1)
  lost <- sum(all.results.young[[3]][[j]][,500]==0)
  text(x=550, y=.95, fixed, pos=4)
  text(x=550, y=.05, lost, pos=4)
  title(main = "Young autosome")
}

#old model
# Y chromosome
for (j in 1:5){ #5 mut rates
  plot(x=1:500, y=all.results.old[[2]][[j]][1,], col=cols[1], type = "l", 
       xlab="Generation", ylab="Frequency of Fixation", xlim=c(0,600), ylim=c(0,1),lwd=.05)
  #i should go from 2 to the number of iterations
  for(i in 2:1000){
    lines(x=1:500, y=all.results.old[[2]][[j]][i,], col=cols[i],lwd=.05)
  }
  fixed <- sum(all.results.old[[2]][[j]][,500]==1)
  lost <- sum(all.results.old[[2]][[j]][,500]==0)
  text(x=550, y=.95, fixed, pos=4)
  text(x=550, y=.05, lost, pos=4)
  title(main = "Old Y chromosome")
}


# X chromosome
for (j in 1:5){ #5 mut rates
  plot(x=1:500, y=all.results.old[[1]][[j]][1,], col=cols[1], type = "l", 
       xlab="Generation", ylab="Frequency of Fixation", xlim=c(0,600), ylim=c(0,1),lwd=.05)
  #i should go from 2 to the number of iterations
  for(i in 2:1000){
    lines(x=1:500, y=all.results.old[[1]][[j]][i,], col=cols[i],lwd=.05)
  }
  fixed <- sum(all.results.old[[1]][[j]][,500]==1)
  lost <- sum(all.results.old[[1]][[j]][,500]==0)
  text(x=550, y=.95, fixed, pos=4)
  text(x=550, y=.05, lost, pos=4)
  title(main = "Old X chromosome")
}


# autosome
for (j in 1:5){ #5 mut rates
  plot(x=1:500, y=all.results.old[[3]][[j]][1,], col=cols[1], type = "l", 
       xlab="Generation", ylab="Frequency of Fixation", xlim=c(0,600), ylim=c(0,1),lwd=.05)
  #i should go from 2 to the number of iterations
  for(i in 2:1000){
    lines(x=1:500, y=all.results.old[[3]][[j]][i,], col=cols[i],lwd=.05)
  }
  fixed <- sum(all.results.old[[3]][[j]][,500]==1)
  lost <- sum(all.results.old[[3]][[j]][,500]==0)
  text(x=550, y=.95, fixed, pos=4)
  text(x=550, y=.05, lost, pos=4)
  title(main = "Old autosome")
}


#aneuploidy=0.1
# Y chromosome
for (j in 1:5){ #5 mut rates
  plot(x=1:500, y=all.results.aneu0.1[[2]][[j]][1,], col=cols[1], type = "l", 
       xlab="Generation", ylab="Frequency of Fixation", xlim=c(0,600), ylim=c(0,1),lwd=.05)
  #i should go from 2 to the number of iterations
  for(i in 2:1000){
    lines(x=1:500, y=all.results.aneu0.1[[2]][[j]][i,], col=cols[i],lwd=.05)
  }
  fixed <- sum(all.results.aneu0.1[[2]][[j]][,500]==1)
  lost <- sum(all.results.aneu0.1[[2]][[j]][,500]==0)
  text(x=550, y=.95, fixed, pos=4)
  text(x=550, y=.05, lost, pos=4)
  title(main = "Y chromosome, aneuploidy rate=0.1")
}



# X chromosome
for (j in 1:5){ #5 mut rates
  plot(x=1:500, y=all.results.aneu0.1[1]][[j]][1,], col=cols[1], type = "l", 
  xlab="Generation", ylab="Frequency of Fixation", xlim=c(0,600), ylim=c(0,1),lwd=.05)
  #i should go from 2 to the number of iterations
  for(i in 2:1000){
    lines(x=1:500, y=all.results.aneu0.1[[1]][[j]][i,], col=cols[i],lwd=.05)
  }
  fixed <- sum(all.results.aneu0.1[[1]][[j]][,500]==1)
  lost <- sum(all.results.aneu0.1[[1]][[j]][,500]==0)
  text(x=550, y=.95, fixed, pos=4)
  text(x=550, y=.05, lost, pos=4)
  title(main = "X chromosome, aneuploidy rate=0.1")
}


# autosome
for (j in 1:5){ #5 mut rates
  plot(x=1:500, y=all.results.aneu0.1[[3]][[j]][1,], col=cols[1], type = "l", 
       xlab="Generation", ylab="Frequency of Fixation", xlim=c(0,600), ylim=c(0,1),lwd=.05)
  #i should go from 2 to the number of iterations
  for(i in 2:1000){
    lines(x=1:500, y=all.results.aneu0.1[[3]][[j]][i,], col=cols[i],lwd=.05)
  }
  fixed <- sum(all.results.aneu0.1[[3]][[j]][,500]==1)
  lost <- sum(all.results.aneu0.1[[3]][[j]][,500]==0)
  text(x=550, y=.95, fixed, pos=4)
  text(x=550, y=.05, lost, pos=4)
  title(main = "Autosome, aneuploidy rate=0.1")
}

#aneuploidy=0.5
# Y chromosome
for (j in 1:5){ #5 mut rates
  plot(x=1:500, y=all.results.aneu0.5[[2]][[j]][1,], col=cols[1], type = "l", 
       xlab="Generation", ylab="Frequency of Fixation", xlim=c(0,600), ylim=c(0,1),lwd=.05)
  #i should go from 2 to the number of iterations
  for(i in 2:1000){
    lines(x=1:500, y=all.results.aneu0.5[[2]][[j]][i,], col=cols[i],lwd=.05)
  }
  fixed <- sum(all.results.aneu0.5[[2]][[j]][,500]==1)
  lost <- sum(all.results.aneu0.5[[2]][[j]][,500]==0)
  text(x=550, y=.95, fixed, pos=4)
  text(x=550, y=.05, lost, pos=4)
  title(main = "Y chromosome, aneuploidy rate=0.5")
}


# X chromosome
for (j in 1:5){ #5 mut rates
  plot(x=1:500, y=all.results.aneu0.5[[1]][[j]][1,], col=cols[1], type = "l", 
       xlab="Generation", ylab="Frequency of Fixation", xlim=c(0,600), ylim=c(0,1),lwd=.05)
  #i should go from 2 to the number of iterations
  for(i in 2:1000){
    lines(x=1:500, y=all.results.aneu0.5[[1]][[j]][i,], col=cols[i],lwd=.05)
  }
  fixed <- sum(all.results.aneu0.5[[1]][[j]][,500]==1)
  lost <- sum(all.results.aneu0.5[[1]][[j]][,500]==0)
  text(x=550, y=.95, fixed, pos=4)
  text(x=550, y=.05, lost, pos=4)
  title(main = "X chromosome, aneuploidy rate=0.5")
}


# autosome
for (j in 1:5){ #5 mut rates
  plot(x=1:500, y=all.results.aneu0.5[[3]][[j]][1,], col=cols[1], type = "l", 
       xlab="Generation", ylab="Frequency of Fixation", xlim=c(0,600), ylim=c(0,1),lwd=.05)
  #i should go from 2 to the number of iterations
  for(i in 2:1000){
    lines(x=1:500, y=all.results.aneu0.5[[3]][[j]][i,], col=cols[i],lwd=.05)
  }
  fixed <- sum(all.results.aneu0.5[[3]][[j]][,500]==1)
  lost <- sum(all.results.aneu0.5[[3]][[j]][,500]==0)
  text(x=550, y=.95, fixed, pos=4)
  text(x=550, y=.05, lost, pos=4)
  title(main = "Autosome, aneuploidy rate=0.5")
}