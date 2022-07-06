library("viridis")
# June 2022
# Annabel Perry
# This script plots the frequency of fusions over all generations in each 
# simulation run under a given set of parameters.

# Key to abbreviations are in filename_key.txt

# Each of these files has the frequency of fusions at each generation in each
# simulation run under the parameters specified in the filename (see 
# filename_key.txt for descriptions of parameters)
setwd("~/06-27-2022")
CYS_HmuLs <- readRDS("CYS_HmuLs.rds")
CYS_LmuHs <- readRDS("CYS_LmuHs.rds")
CYL_HmuLs <- readRDS("CYL_HmuLs.rds")
CYL_LmuHs <- readRDS("CYL_LmuHs.rds")
CXS_HmuLs <- readRDS("CXS_HmuLs.rds")
CXS_LmuHs <- readRDS("CXS_LmuHs.rds")
CXL_HmuLs <- readRDS("CXL_HmuLs.rds")
CXL_LmuHs <- readRDS("CXL_LmuHs.rds")
AYS_HmuLs <- readRDS("AYS_HmuLs.rds")
AYS_LmuHs <- readRDS("AYS_LmuHs.rds")
AYL_HmuLs <- readRDS("AYL_HmuLs.rds")
AYL_LmuHs <- readRDS("AYL_LmuHs.rds")
AXS_HmuLs <- readRDS("AXS_HmuLs.rds")
AXS_LmuHs <- readRDS("AXS_LmuHs.rds")
AXL_HmuLs <- readRDS("AXL_HmuLs.rds")
AXL_LmuHs <- readRDS("AXL_LmuHs.rds")

# Plot frequency of fusions as a function of generation, with simulations with 
# the same meiotic and fusion parameters plotted in same panel. On each plot,
# each line represents a distinct simulation.

# Sample a new color for each simulation
col.vec <- viridis(ncol(CYS_HmuLs))

###################################### CYS ##################################### 
par(mfrow = c(2,2), cex = 0.75)

plot(CYS_HmuLs[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "Chiasmatic Y-Small, High mu Low s",
     col = col.vec[1])
for(i in 2:ncol(CYS_HmuLs)){
  lines(y = CYS_HmuLs[,i], x = 1:nrow(CYS_HmuLs), col = col.vec[i])
}

plot(CYS_LmuHs[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "Chiasmatic Y-Small, Low mu High s",
     cex = 1,
     col = col.vec[1])
for(i in 2:ncol(CYS_LmuHs)){
  lines(y = CYS_LmuHs[,i], x = 1:nrow(CYS_LmuHs), col = col.vec[i])
}

###################################### CYL ##################################### 
plot(CYL_HmuLs[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "Chiasmatic Y-Large, High mu Low s",
     col = col.vec[1])
for(i in 2:ncol(CYL_HmuLs)){
  lines(y = CYL_HmuLs[,i], x = 1:nrow(CYL_HmuLs), col = col.vec[i])
}

plot(CYL_LmuHs[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "Chiasmatic Y-Large, Low mu High s",
     cex = 1,
     col = col.vec[1])
for(i in 2:ncol(CYL_LmuHs)){
  lines(y = CYL_LmuHs[,i], x = 1:nrow(CYL_LmuHs), col = col.vec[i])
}

###################################### CXS ##################################### 
par(mfrow = c(2,2), cex = 0.75)

plot(CXS_HmuLs[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "Chiasmatic X-Small, High mu Low s",
     col = col.vec[1])
for(i in 2:ncol(CXS_HmuLs)){
  lines(y = CXS_HmuLs[,i], x = 1:nrow(CXS_HmuLs), col = col.vec[i])
}

plot(CXS_LmuHs[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "Chiasmatic X-Small, Low mu High s",
     cex = 1,
     col = col.vec[1])
for(i in 2:ncol(CXS_LmuHs)){
  lines(y = CXS_LmuHs[,i], x = 1:nrow(CXS_LmuHs), col = col.vec[i])
}

###################################### CXL ##################################### 
#par(mfrow = c(1,2), cex = 0.75)

plot(CXL_HmuLs[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "Chiasmatic X-Large, High mu Low s",
     col = col.vec[1])
for(i in 2:ncol(CXL_HmuLs)){
  lines(y = CXL_HmuLs[,i], x = 1:nrow(CXL_HmuLs), col = col.vec[i])
}

plot(CXL_LmuHs[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "Chiasmatic X-Large, Low mu High s",
     cex = 1,
     col = col.vec[1])
for(i in 2:ncol(CXL_LmuHs)){
  lines(y = CXL_LmuHs[,i], x = 1:nrow(CXL_LmuHs), col = col.vec[i])
}

###################################### AYS ##################################### 
par(mfrow = c(2,2), cex = 0.75)

plot(AYS_HmuLs[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "Achiasmatic Y-Small, High mu Low s",
     col = col.vec[1])
for(i in 2:ncol(AYS_HmuLs)){
  lines(y = AYS_HmuLs[,i], x = 1:nrow(AYS_HmuLs), col = col.vec[i])
}

plot(AYS_LmuHs[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "Achiasmatic Y-Small, Low mu High s",
     cex = 1,
     col = col.vec[1])
for(i in 2:ncol(AYS_LmuHs)){
  lines(y = AYS_LmuHs[,i], x = 1:nrow(AYS_LmuHs), col = col.vec[i])
}

###################################### AYL ##################################### 
#par(mfrow = c(1,2), cex = 0.75)

plot(AYL_HmuLs[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "Achiasmatic Y-Large, High mu Low s",
     col = col.vec[1])
for(i in 2:ncol(AYL_HmuLs)){
  lines(y = AYL_HmuLs[,i], x = 1:nrow(AYL_HmuLs), col = col.vec[i])
}

plot(AYL_LmuHs[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "Achiasmatic Y-Large, Low mu High s",
     cex = 1,
     col = col.vec[1])
for(i in 2:ncol(AYL_LmuHs)){
  lines(y = AYL_LmuHs[,i], x = 1:nrow(AYL_LmuHs), col = col.vec[i])
}

###################################### AXS ##################################### 
par(mfrow = c(2,2), cex = 0.75)

plot(AXS_HmuLs[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "Achiasmatic X-Small, High mu Low s",
     col = col.vec[1])
for(i in 2:ncol(AXS_HmuLs)){
  lines(y = AXS_HmuLs[,i], x = 1:nrow(AXS_HmuLs), col = col.vec[i])
}

plot(AXS_LmuHs[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "Achiasmatic X-Small, Low mu High s",
     cex = 1,
     col = col.vec[1])
for(i in 2:ncol(AXS_LmuHs)){
  lines(y = AXS_LmuHs[,i], x = 1:nrow(AXS_LmuHs), col = col.vec[i])
}

###################################### AXL ##################################### 
#par(mfrow = c(1,2), cex = 0.75)

plot(AXL_HmuLs[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "Achiasmatic X-Large, High mu Low s",
     col = col.vec[1])
for(i in 2:ncol(AXL_HmuLs)){
  lines(y = AXL_HmuLs[,i], x = 1:nrow(AXL_HmuLs), col = col.vec[i])
}

plot(AXL_LmuHs[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "Achiasmatic X-Large, Low mu High s",
     cex = 1,
     col = col.vec[1])
for(i in 2:ncol(AXL_LmuHs)){
  lines(y = AXL_LmuHs[,i], x = 1:nrow(AXL_LmuHs), col = col.vec[i])
}
