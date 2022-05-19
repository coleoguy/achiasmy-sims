# Find why ALL males are marked as going extinct

# From GetPop()
N = 10
s = 1
pop <- GetPop(N)
fits <- GetFit(pop, s)

PaternalGenealogy <- rep(NA, N)
for(i in 1:length(pop)){
  if(pop[[i]][2,13]){
    PaternalGenealogy[i] <- i
  }
}
n.males <- sum(!is.na(PaternalGenealogy))
LineageStatus <- data.frame(matrix(c(PaternalGenealogy[!is.na(PaternalGenealogy)], 
                                     rep(0, n.males), rep(0, n.males)), 
                                   nrow = n.males, ncol = 3))
colnames(LineageStatus) <- c("Paternal Lineage", "Extinct?", "Fused?")

# From Maury()
sexes <- rep("fem",N)
for(i in 1:length(sexes)){
  if(pop[[i]][2,13]){
    sexes[i] <- "mal"
  }
}
if(sum(sexes=="fem") > 1){
  moms <- sample((1:N)[sexes=="fem"], prob=fits[sexes=="fem"], 
                 size=N, replace = T)
}else{
  moms <- rep((1:N)[sexes=="fem"], N)
}

if(sum(sexes=="mal") > 1){
  dads <- sample((1:N)[sexes=="mal"], prob=fits[sexes=="mal"], 
                 size=N, replace = T)
}else{
  dads <- rep((1:N)[sexes=="mal"], N)
}
# Record the lineages of all individuals who were NOT selected as dads as extinct
incels <- PaternalGenealogy[PaternalGenealogy %!in% dads & !is.na(PaternalGenealogy)]
LineageStatus$`Extinct?`[LineageStatus$`Paternal Lineage` %in% incels] <- 1
