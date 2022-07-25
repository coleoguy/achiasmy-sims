# July 2022
# Annabel Perry
# This script runs tests to answer the following two questions:
# 1. For each combination of sex chromosomes, selection coefficients, and auto-
#    somal sizes, do achiasmatic simulations have a LOWER proportion of fusions
#    than chiasmatic simulations?
# 2. For each combination sex chromosomes, selection coefficients, and chiasmy
#    categories, do simulations which permit fusions of large autosomes have a
#    LOWER proportion of fusions than simulations which permit fusions of small
#    autosomes?
# ... and visualizes the results

library("viridis")
library("dplyr")
library("stringr")

############################### Data Preparation ###############################
# Read in data files for each combination of meiotic type, sex chromosome, 
# autosome size, and selection coefficient
AXL0.3 <- readRDS("AXL_s=0.3_2022-07-23.rds")
AXL0.25 <- readRDS("AXL_s=0.25_2022-07-22.rds")
AXL0.35 <- readRDS("AXL_s=0.35_2022-07-24.rds")
AXS0.3 <- readRDS("AXS_s=0.3_2022-07-23.rds")
AXS0.25 <- readRDS("AXS_s=0.25_2022-07-22.rds")
AXS0.35 <- readRDS("AXS_s=0.35_2022-07-24.rds")
AYL0.3 <- readRDS("AYL_s=0.3_2022-07-23.rds")
AYL0.25 <- readRDS("AYL_s=0.25_2022-07-22.rds")
AYL0.35 <- readRDS("AYL_s=0.35_2022-07-24.rds")
AYS0.3 <- readRDS("AYS_s=0.3_2022-07-23.rds")
AYS0.25 <- readRDS("AYS_s=0.25_2022-07-22.rds")
AYS0.35 <- readRDS("AYS_s=0.35_2022-07-24.rds")
CXL0.3 <- readRDS("CXL_s=0.3_2022-07-23.rds")
CXL0.25 <- readRDS("CXL_s=0.25_2022-07-22.rds")
CXL0.35 <- readRDS("CXL_s=0.35_2022-07-24.rds")
CXS0.3 <- readRDS("CXS_s=0.3_2022-07-23.rds")
CXS0.25 <- readRDS("CXS_s=0.25_2022-07-21.rds")
CXS0.35 <- readRDS("CXS_s=0.35_2022-07-24.rds")
CYL0.3 <- readRDS("CYL_s=0.3_2022-07-22.rds")
CYL0.25 <- readRDS("CYL_s=0.25_2022-07-21.rds")
CYL0.35 <- readRDS("CYL_s=0.35_2022-07-23.rds")
CYS0.3 <- readRDS("CYS_s=0.3_2022-07-22.rds")
CYS0.25 <- readRDS("CYS_s=0.25_2022-07-21.rds")
CYS0.35 <- readRDS("CYS_s=0.35_2022-07-23.rds")

# Create a dataframe with the frequency of fusions in the final generation
# of each INDIVIDUAL simulation
num_sims <- 1000
gen_no <- 1000
s_coeffs <- c(0.25, 0.3, 0.35)

PropFusions_EachSim <- data.frame(
  Simulation = rep(1:num_sims, each = 2*2*2*length(s_coeffs)),
  Meiotic_Type = rep(rep(c("A", "C"), each = 2*2*length(s_coeffs)), 
                     times = num_sims),
  Sex_Chr = rep(rep(c("X", "Y"), each = 2*length(s_coeffs)), times = 2*num_sims),
  Autosome_Size = rep(c("L", "S"), times = 2*2*length(s_coeffs)*num_sims),
  s = rep(rep(s_coeffs, each = 2), times = 2*2*num_sims),
  PropFusions = rep(0, times = 2*2*2*length(s_coeffs)*num_sims)
)

# Create dataframe in which to store avg proportion of FUSIONS at final
# generation across ALL simulations. The order of the variables within 
# Meiotic_Type and Autosome_Size are arranged in a specific order to enable
# plotting later, so do not change the order
PropFusions_Avg <- data.frame(
  Meiotic_Type = rep(c("A", "C"), each = 2*2*length(s_coeffs)),
  Sex_Chr = rep(rep(c("X", "Y"), each = 2*length(s_coeffs)), times = 2),
  Autosome_Size = rep(c("L", "S"), times = 2*2*length(s_coeffs)),
  s = rep(rep(s_coeffs, each = 2), times = 4),
  PropFusions_Avg = rep(0, times = 2*2*2*length(s_coeffs)),
  LowerCI = rep(0, times = 2*2*2*length(s_coeffs)),
  UpperCI = rep(0, times = 2*2*2*length(s_coeffs))
)

# Fill each row of the data frame with the proportion of fusions of the 
# corresponding simulation
num_conditions <- 8
for(i in 1:length(s_coeffs)){
  s <- s_coeffs[i]
  for(cond in 1:num_conditions){ 
    # Retrieve a unique combination of meiotic, fusion, and autosomal size 
    # parameters
    if(cond <= 4){
      meiotic_type <- "C"
    }else{
      meiotic_type <- "A" 
    }
    
    if(cond %in% c(1, 2, 5, 6)){
      fus.type <- "Y"
    }else{
      fus.type <- "X"
    } 
    
    if(cond%%2){
      size <- "S"
    }else{
      size <- "L"
    } 
    
    # Prep dataframe name based on current parameters 
    dataname <- paste(c(meiotic_type, fus.type, size, s), collapse = "")
    print(dataname)
    
    
    
    # For each simulation run under the current combination of parameters, 
    # retrieve the proportion of fusions at the FINAL generation
    for(sim in 1:num_sims){
      # Retrieve and store the row of the fusion proportion dataframe 
      # corresponding to the current simulation and set of parameters
      current_row <- which(PropFusions_EachSim$Simulation == sim &
                             PropFusions_EachSim$Meiotic_Type == meiotic_type &
                             PropFusions_EachSim$Sex_Chr == fus.type &
                             PropFusions_EachSim$Autosome_Size == size &
                             PropFusions_EachSim$s == s)
      # Retrieve and store the proportion of fusions at the final generation of
      # the current simulation and dataframe
      prop_fusions <- eval(parse(text = paste(dataname, "[[", sim, "]]", 
                                              sep = "")))[gen_no]
      PropFusions_EachSim$PropFusions[current_row] <- prop_fusions
    }
    
    # Calculate and store the average proportion of fusions across all 
    # simulations as well as the lower and upper confidence intervals for this
    # current set of conditions
    
    # First, average the proportion of fusions at final generation for all 
    # simulations run under the current set of parameters
    all_props <- PropFusions_EachSim$PropFusions[
      PropFusions_EachSim$Meiotic_Type == meiotic_type &
        PropFusions_EachSim$Sex_Chr == fus.type &
        PropFusions_EachSim$Autosome_Size == size &
        PropFusions_EachSim$s == s]
    avg_props <- sum(all_props)/num_sims
    # Store average proportion of fusions across all simulations which were
    # run under the current set of parameters
    cur_pars <- which(PropFusions_Avg$Meiotic_Type == meiotic_type &
                        PropFusions_Avg$Sex_Chr == fus.type &
                        PropFusions_Avg$Autosome_Size == size &
                        PropFusions_Avg$s == s)
    PropFusions_Avg$PropFusions_Avg[cur_pars] <- avg_props
    
    # Finally, use the avg proportion of fusions at final generation  and the 
    # standard error of the mean to calculate the upper and lower 95% Bayesian 
    # confidence intervals for the count of fixations
    PropFusions_Avg$LowerCI[cur_pars] <- avg_props - 1.96*(sd(all_props)/sqrt(4))
    
    PropFusions_Avg$UpperCI[cur_pars] <- avg_props + 1.96*(sd(all_props)/sqrt(4))
  }
}

# Create dataframe where parameters are merged together
stat_df <- data.frame(
  Parameters = apply(PropFusions_EachSim[,2:5], MARGIN = 1, FUN = paste, 
                     collapse = ""),
  PropFusions = PropFusions_EachSim$PropFusions
)
# Run ANOVA to compare all sets of parameters
all_tests <- TukeyHSD(aov(lm(stat_df$PropFusions ~ stat_df$Parameters)))[[1]]

pars <- dimnames(all_tests)[[1]]

################################## Question 1 ################################## 
# Retain only the p-values for achiasmy-chiasmy comparisons where all other 
# parameters are held constant
Q1 <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(Q1) <- c("Parameters", "P_Value")
for(row in 1:length(pars)){
  # If this is an achiasmy-chiasmy comparison BUT all other parameters are
  # constant, add to the Question 1 dataframe
  if(!is.na(str_locate(pars[row], "A")[1]) & 
     !is.na(str_locate(pars[row], "C")[1]) &
     (substr(pars[row], start = 2, stop = str_locate(pars[row], "-") - 1) ==
      substr(pars[row], start = str_locate(pars[row], "-") + 2, 
             stop = nchar(pars[row])))){
       Q1 <- rbind(Q1, data.frame(
         Parameters = pars[row], 
         P_Value = all_tests[row, 4])
         )
  }
}

################################## Question 2 ################################## 
# Retain only the p-values for large-small comparsions where all other 
# parameters are held constant
Q2 <- data.frame(matrix(nrow = 0, ncol = 2))
colnames(Q2) <- c("Parameters", "P_Value")
for(row in 1:length(pars)){
  # If this is an large-small comparison BUT all other parameters are
  # constant, add to the Question 2 dataframe
  if(!is.na(str_locate(pars[row], "L")[1]) & 
     !is.na(str_locate(pars[row], "S")[1]) &
     (substr(pars[row], start = 1, stop = 2) ==
      substr(pars[row], start = str_locate(pars[row], "-") + 1, 
             stop = str_locate(pars[row], "-") + 2)) &
     (substr(pars[row], start = 4, stop = str_locate(pars[row], "-") - 1) ==
      substr(pars[row], start = str_locate(pars[row], "-") + 4, 
             stop = nchar(pars[row])))){
    Q2 <- rbind(Q2, data.frame(
      Parameters = pars[row], 
      P_Value = all_tests[row, 4])
    )
  }
}

################################### Plotting ################################### 
# Initialize the X plot with the average proportion of fusions where the
# selection coefficient is 0 and the sex chromosome is X. The script later will 
# add points for every selection coefficient.
cur_pars <- which(PropFusions_Avg$Sex_Chr == "X" & 
                    PropFusions_Avg$s == s_coeffs[1])
props <- PropFusions_Avg$PropFusions_Avg[cur_pars]
plot_listX <- list("Large" = props[1], 
                   "Small" = props[2], 
                   "Large" = props[3],
                   "Small" = props[4]
)

# Get a unique color for each selection coefficient
point_colors <- viridis(length(s_coeffs))
  
# Plot the average proportion of fusions at final generation of all simulations 
# where s = 0 and sex_chr = "X"
par(mar = c(5.1, 4.1, 4.1, 8.1), xpd = TRUE, cex = 0.8, cex.axis = 1)
stripchart(plot_listX,
           pch = 16,
           yaxt = "n",
           vertical = TRUE,
           col = point_colors[1],
           # Set the limits of the plot based on the most extreme possible 
           # confidence interval values
           ylim = c(PropFusions_Avg$LowerCI[which.min(PropFusions_Avg$LowerCI)],
                    PropFusions_Avg$UpperCI[which.max(PropFusions_Avg$UpperCI)])
)
axis(side = 2, at = (0:10)/10)
title(main = "Proportion of Fusions to X Chromosome at Generation 1000, 
      Averaged Across 1000 Simulations",
      sub = paste(c("Achiasmatic","Chiasmatic"), 
                  collapse = "                                   "),
      adj = 0.5)

legend(x = 4.25, y = 0.8, inset=c(-0.2,0), legend = s_coeffs, col = point_colors, 
       pch = 16, cex = 0.8, title = "Selection Coefficient")

# Add confidence intervals for the average proportion of fusions
CIs <- data.frame(
  lower = PropFusions_Avg$LowerCI[cur_pars],
  upper = PropFusions_Avg$UpperCI[cur_pars]
)
for(i in 1:nrow(CIs)){
  # If the upper and lower confidence intervals have the same value, simply
  # draw a line through the mean
  if(CIs$lower[i] == CIs$upper[i]){
    lines(x = c(i - 0.1, i + 0.1), y = c(CIs$lower[i],CIs$upper[i]),
          lwd = 1, point_colors[1])
  }else{
    arrows(y0 = CIs$lower[i], y1 = CIs$upper[i], x0 = i, x1 = i, 
           angle = 90, code = 3, lwd = 1, length = 0.1,
           col = point_colors[1])
    
  }
}

# Add points and confidence intervals corresponding to each selection 
# coefficient
# Skip the first selection coefficient (since you have already added points
# corresponding to this selection coefficient)
for(sc in 2:length(s_coeffs)){
  # Set current parameters based on current selection coefficient
  cur_pars <- which(PropFusions_Avg$Sex_Chr == "X" &
                      PropFusions_Avg$s == s_coeffs[sc])
  # Obtain the proportion of fusions for each combination of chiasmy type and 
  # autosome size for the X chromosome and current selection coefficient. Use 
  # the color appropriate for the current s
  points(y = PropFusions_Avg$PropFusions_Avg[cur_pars], x = 1:4, pch = 16, 
         col = point_colors[sc])
  CIs <- data.frame(
    lower = PropFusions_Avg$LowerCI[cur_pars],
    upper = PropFusions_Avg$UpperCI[cur_pars]
    )
  for(i in 1:nrow(CIs)){
    # If the upper and lower confidence intervals have the same value, simply
    # draw a line through the mean
    if(CIs$lower[i] == CIs$upper[i]){
      lines(x = c(i - 0.1, i + 0.1), y = c(CIs$lower[i],CIs$upper[i]),
            lwd = 1, col = point_colors[sc])
    }else{
      arrows(y0 = CIs$lower[i], y1 = CIs$upper[i], x0 = i, x1 = i, 
             angle = 90, code = 3, lwd = 1, length = 0.1,
             col = point_colors[sc])
      
    }
  }
}


# Repeat the initialization process for the Y chromosome plot
cur_pars <- which(PropFusions_Avg$Sex_Chr == "Y" & 
                    PropFusions_Avg$s == s_coeffs[1])
props <- PropFusions_Avg$PropFusions_Avg[cur_pars]
plot_listY <- list("Large" = props[1], 
                   "Small" = props[2], 
                   "Large" = props[3],
                   "Small" = props[4]
)

# Plot the average proportion of fusions at final generation of all simulations 
# where s = 0 and sex_chr = "Y"
par(mar = c(5.1, 4.1, 4.1, 8.1), xpd = TRUE, cex = 0.8, cex.axis = 1)
stripchart(plot_listY,
           pch = 16,
           yaxt = "n",
           vertical = TRUE,
           col = point_colors[1],
           # Set the limits of the plot based on the most extreme possible 
           # confidence interval values
           ylim = c(PropFusions_Avg$LowerCI[which.min(PropFusions_Avg$LowerCI)],
                    PropFusions_Avg$UpperCI[which.max(PropFusions_Avg$UpperCI)])
)
axis(side = 2, at = (0:10)/10)
title(main = "Proportion of Fusions to Y Chromosome at Generation 1000, 
      Averaged Across 1000 Simulations",
      sub = paste(c("Achiasmatic","Chiasmatic"), 
                  collapse = "                                   "),
      adj = 0.5)

legend(x = 4.25, y = 0.8, inset=c(-0.2,0), legend = s_coeffs, col = point_colors, 
       pch = 16, cex = 0.8, title = "Selection Coefficient")

# Add confidence intervals for the average proportion of fusions
CIs <- data.frame(
  lower = PropFusions_Avg$LowerCI[cur_pars],
  upper = PropFusions_Avg$UpperCI[cur_pars]
)
for(i in 1:nrow(CIs)){
  # If the upper and lower confidence intervals have the same value, simply
  # draw a line through the mean
  if(CIs$lower[i] == CIs$upper[i]){
    lines(x = c(i - 0.1, i + 0.1), y = c(CIs$lower[i],CIs$upper[i]),
          lwd = 1, point_colors[1])
  }else{
    arrows(y0 = CIs$lower[i], y1 = CIs$upper[i], x0 = i, x1 = i, 
           angle = 90, code = 3, lwd = 1, length = 0.1,
           col = point_colors[1])
    
  }
}

# Add points and confidence intervals corresponding to each selection 
# coefficient
# Skip the first selection coefficient (since you have already added points
# corresponding to this selection coefficient)
for(sc in 2:length(s_coeffs)){
  # Set current parameters based on current selection coefficient
  cur_pars <- which(PropFusions_Avg$Sex_Chr == "Y" &
                      PropFusions_Avg$s == s_coeffs[sc])
  # Obtain the proportion of fusions for each combination of chiasmy type and 
  # autosome size for the Y chromosome and current selection coefficient. Use 
  # the color appropriate for the current s
  points(y = PropFusions_Avg$PropFusions_Avg[cur_pars], x = 1:4, pch = 16, 
         col = point_colors[sc])
  CIs <- data.frame(
    lower = PropFusions_Avg$LowerCI[cur_pars],
    upper = PropFusions_Avg$UpperCI[cur_pars]
  )
  for(i in 1:nrow(CIs)){
    # If the upper and lower confidence intervals have the same value, simply
    # draw a line through the mean
    if(CIs$lower[i] == CIs$upper[i]){
      lines(x = c(i - 0.1, i + 0.1), y = c(CIs$lower[i],CIs$upper[i]),
            lwd = 1, col = point_colors[sc])
    }else{
      arrows(y0 = CIs$lower[i], y1 = CIs$upper[i], x0 = i, x1 = i, 
             angle = 90, code = 3, lwd = 1, length = 0.1,
             col = point_colors[sc])
      
    }
  }
}

