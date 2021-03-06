library("viridis")
# July 2022
# Annabel Perry
# This script plots the frequency of fusions over all generations in each 
# simulation run with an empirically-derived mutation rate.

# Key to abbreviations are in filename_key.txt

# Each of these files has the frequency of fusions at each generation in each
# simulation run under the parameters specified in the filename (see 
# filename_key.txt for descriptions of parameters)

# Plot frequency of fusions as a function of generation.

################################## CYS - Pilot #################################
CYS_10 <- readRDS("CYS_5e-10.rds")
CYS_9 <- readRDS("CYS_5e-9.rds")
CYS_8 <- readRDS("CYS_5e-8.rds")

# Sample a new color for each simulation
col.vec <- viridis(ncol(CYS_10))

par(mfrow = c(1,3), cex = 0.6)

plot(CYS_10[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "Chiasmatic Y-Small: mu = 5e^(-10)",
     col = col.vec[1])
for(i in 2:ncol(CYS_10)){
  lines(y = CYS_10[,i], x = 1:nrow(CYS_10), col = col.vec[i])
}

plot(CYS_9[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "mu = 5e^(-9)",
     cex = 1,
     col = col.vec[1])
for(i in 2:ncol(CYS_9)){
  lines(y = CYS_9[,i], x = 1:nrow(CYS_9), col = col.vec[i])
}

plot(CYS_8[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "mu = 5e^(-8)",
     cex = 1,
     col = col.vec[1])
for(i in 2:ncol(CYS_8)){
  lines(y = CYS_8[,i], x = 1:nrow(CYS_8), col = col.vec[i])
}




###################################### CYS #####################################
par(mfrow = c(1,3), cex = 0.6)
CYS_10 <- readRDS("CYS_5e-10.rds")
CYS_9 <- readRDS("CYS_5e-09.rds")
CYS_8 <- readRDS("CYS_5e-08.rds")

plot(CYS_10[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "Chiasmatic Y-Small: mu = 5e^(-10)",
     col = col.vec[1])
for(i in 2:ncol(CYS_10)){
  lines(y = CYS_10[,i], x = 1:nrow(CYS_10), col = col.vec[i])
}

plot(CYS_9[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "mu = 5e^(-9)",
     cex = 1,
     col = col.vec[1])
for(i in 2:ncol(CYS_9)){
  lines(y = CYS_9[,i], x = 1:nrow(CYS_9), col = col.vec[i])
}

plot(CYS_8[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "mu = 5e^(-8)",
     cex = 1,
     col = col.vec[1])
for(i in 2:ncol(CYS_8)){
  lines(y = CYS_8[,i], x = 1:nrow(CYS_8), col = col.vec[i])
}


###################################### CXS #####################################
CXS_10 <- readRDS("CXS_5e-10.rds")
CXS_9 <- readRDS("CXS_5e-09.rds")
CXS_8 <- readRDS("CXS_5e-08.rds")

col.vec <- viridis(ncol(CXS_10))

par(mfrow = c(1,3), cex = 0.6)

plot(CXS_10[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "Chiasmatic X-Small: mu = 5e^(-10)",
     col = col.vec[1])
for(i in 2:ncol(CXS_10)){
  lines(y = CXS_10[,i], x = 1:nrow(CXS_10), col = col.vec[i])
}

plot(CXS_9[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "mu = 5e^(-9)",
     cex = 1,
     col = col.vec[1])
for(i in 2:ncol(CXS_9)){
  lines(y = CXS_9[,i], x = 1:nrow(CXS_9), col = col.vec[i])
}

plot(CXS_8[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "mu = 5e^(-8)",
     cex = 1,
     col = col.vec[1])
for(i in 2:ncol(CXS_8)){
  lines(y = CXS_8[,i], x = 1:nrow(CXS_8), col = col.vec[i])
}








###################################### CYL #####################################
CYL_10 <- readRDS("CYL_5e-10.rds")
CYL_9 <- readRDS("CYL_5e-09.rds")
CYL_8 <- readRDS("CYL_5e-08.rds")

col.vec <- viridis(ncol(CYL_10))

par(mfrow = c(1,3), cex = 0.6)

plot(CYL_10[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "Chiasmatic Y-Large: mu = 5e^(-10)",
     col = col.vec[1])
for(i in 2:ncol(CYL_10)){
  lines(y = CYL_10[,i], x = 1:nrow(CYL_10), col = col.vec[i])
}

plot(CYL_9[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "mu = 5e^(-9)",
     cex = 1,
     col = col.vec[1])
for(i in 2:ncol(CYL_9)){
  lines(y = CYL_9[,i], x = 1:nrow(CYL_9), col = col.vec[i])
}

plot(CYL_8[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "mu = 5e^(-8)",
     cex = 1,
     col = col.vec[1])
for(i in 2:ncol(CYL_8)){
  lines(y = CYL_8[,i], x = 1:nrow(CYL_8), col = col.vec[i])
}


###################################### CXL #####################################
CXL_10 <- readRDS("CXL_5e-10.rds")
CXL_9 <- readRDS("CXL_5e-09.rds")
CXL_8 <- readRDS("CXL_5e-08.rds")

col.vec <- viridis(ncol(CXL_10))

par(mfrow = c(1,3), cex = 0.6)

plot(CXL_10[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "Chiasmatic X-Large: mu = 5e^(-10)",
     col = col.vec[1])
for(i in 2:ncol(CXL_10)){
  lines(y = CXL_10[,i], x = 1:nrow(CXL_10), col = col.vec[i])
}

plot(CXL_9[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "mu = 5e^(-9)",
     cex = 1,
     col = col.vec[1])
for(i in 2:ncol(CXL_9)){
  lines(y = CXL_9[,i], x = 1:nrow(CXL_9), col = col.vec[i])
}

plot(CXL_8[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "mu = 5e^(-8)",
     cex = 1,
     col = col.vec[1])
for(i in 2:ncol(CXL_8)){
  lines(y = CXL_8[,i], x = 1:nrow(CXL_8), col = col.vec[i])
}


###################################### AYS #####################################
AYS_10 <- readRDS("AYS_5e-10.rds")
AYS_9 <- readRDS("AYS_5e-09.rds")
AYS_8 <- readRDS("AYS_5e-08.rds")

col.vec <- viridis(ncol(AYS_10))

par(mfrow = c(1,3), cex = 0.6)

plot(AYS_10[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "Achiasmatic Y-Sm: mu = 5e^(-10)",
     col = col.vec[1])
for(i in 2:ncol(AYS_10)){
  lines(y = AYS_10[,i], x = 1:nrow(AYS_10), col = col.vec[i])
}

plot(AYS_9[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "mu = 5e^(-9)",
     cex = 1,
     col = col.vec[1])
for(i in 2:ncol(AYS_9)){
  lines(y = AYS_9[,i], x = 1:nrow(AYS_9), col = col.vec[i])
}

plot(AYS_8[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "mu = 5e^(-8)",
     cex = 1,
     col = col.vec[1])
for(i in 2:ncol(AYS_8)){
  lines(y = AYS_8[,i], x = 1:nrow(AYS_8), col = col.vec[i])
}


###################################### AXS #####################################
AXS_10 <- readRDS("AXS_5e-10.rds")
AXS_9 <- readRDS("AXS_5e-09.rds")
AXS_8 <- readRDS("AXS_5e-08.rds")

col.vec <- viridis(ncol(AXS_10))

par(mfrow = c(1,3), cex = 0.6)

plot(AXS_10[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "Achiasmatic X-Sm: mu = 5e^(-10)",
     col = col.vec[1])
for(i in 2:ncol(AXS_10)){
  lines(y = AXS_10[,i], x = 1:nrow(AXS_10), col = col.vec[i])
}

plot(AXS_9[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "mu = 5e^(-9)",
     cex = 1,
     col = col.vec[1])
for(i in 2:ncol(AXS_9)){
  lines(y = AXS_9[,i], x = 1:nrow(AXS_9), col = col.vec[i])
}

plot(AXS_8[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "mu = 5e^(-8)",
     cex = 1,
     col = col.vec[1])
for(i in 2:ncol(AXS_8)){
  lines(y = AXS_8[,i], x = 1:nrow(AXS_8), col = col.vec[i])
}


###################################### AYL #####################################
AYL_10 <- readRDS("AYL_5e-10.rds")
AYL_9 <- readRDS("AYL_5e-09.rds")
AYL_8 <- readRDS("AYL_5e-08.rds")

col.vec <- viridis(ncol(AYL_10))

par(mfrow = c(1,3), cex = 0.6)

plot(AYL_10[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "Achiasmatic Y-La: mu = 5e^(-10)",
     col = col.vec[1])
for(i in 2:ncol(AYL_10)){
  lines(y = AYL_10[,i], x = 1:nrow(AYL_10), col = col.vec[i])
}

plot(AYL_9[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "mu = 5e^(-9)",
     cex = 1,
     col = col.vec[1])
for(i in 2:ncol(AYL_9)){
  lines(y = AYL_9[,i], x = 1:nrow(AYL_9), col = col.vec[i])
}

plot(AYL_8[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "mu = 5e^(-8)",
     cex = 1,
     col = col.vec[1])
for(i in 2:ncol(AYL_8)){
  lines(y = AYL_8[,i], x = 1:nrow(AYL_8), col = col.vec[i])
}


###################################### AXL #####################################
AXL_10 <- readRDS("AXL_5e-10.rds")
AXL_9 <- readRDS("AXL_5e-09.rds")
AXL_8 <- readRDS("AXL_5e-08.rds")

col.vec <- viridis(ncol(AXL_10))

par(mfrow = c(1,3), cex = 0.6)

plot(AXL_10[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "Achiasmatic X-La: mu = 5e^(-10)",
     col = col.vec[1])
for(i in 2:ncol(AXL_10)){
  lines(y = AXL_10[,i], x = 1:nrow(AXL_10), col = col.vec[i])
}

plot(AXL_9[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "mu = 5e^(-9)",
     cex = 1,
     col = col.vec[1])
for(i in 2:ncol(AXL_9)){
  lines(y = AXL_9[,i], x = 1:nrow(AXL_9), col = col.vec[i])
}

plot(AXL_8[,1], type = "l", ylim=c(0,1), 
     xlab = "Generation", ylab = "Frequency Fusions", 
     main = "mu = 5e^(-8)",
     cex = 1,
     col = col.vec[1])
for(i in 2:ncol(AXL_8)){
  lines(y = AXL_8[,i], x = 1:nrow(AXL_8), col = col.vec[i])
}


