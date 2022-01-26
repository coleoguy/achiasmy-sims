#fraction of frequency that fixed was taken to make single plot for 
#comparison purposes

library(viridis)
cols <- viridis(5)

load(".../results/new. achiasmatic.Rdata.Rdata")



# Define the position of tick marks on x-axis
v1 <- c(1,2,3,4,5)
#define label for axis points
v2<- c("1.45e-12","1.45e-11", "1.45e-10","1.45e-9", "1.45e-8") 

#Y chromosome 
x<- c(1, 2, 3, 4, 5)
y<- c(0.735, 0.705, 0.266, 0.004, 0)

#X Chromosome
x2<- c(1, 2, 3, 4, 5)
y2<- c(0.079, 0.051, 0.007, 0, 0)

#autosome
x3<- c(1, 2, 3, 4, 5)
y3<- c(0, 0, 0, 0, 0)

plot(x,y, ylab="Probability", xlab = "Mutation Rate", 
     main = "Probability Achiasmatic Meiosis Mutation Fixes", col=cols[1],
     pch=15, type= "b", axes = FALSE) 
axis(side=1, at= v1, labels = v2)
axis(side=2, at=seq(0, 0.9, by=0.1))
box()

points(x2, y2, pch=15, col=cols[2], type="b")
points(x3, y3, pch=15, col=cols[3], type="b")


legend(3.5, 0.7, legend=c("Y Chromosome", "X Chromosome", "Autosome"),
       col=c(cols[1], cols[2], cols[3]), lty=1, cex=0.8)


