# The purpose of this code is to compare the impact of
# achiasmatic meiosis on fusions of large or small autosomes
# with sex chromosomes we will build a canonically inspired
# model where the sex chromosome has an SDL and the autosomes have 
# SA loci

#  sex chromosome         autosome 1         autosome 2
#  1-12(PAR) SDL 14-25    26-35 SA 37-50     51-60 SA 62-100
#

# SDL 0=X 1=Y
# SA 0=female benefit; 1=male benefit
# All other loci are general fitness loci with values from 0-1
# reflecting their fitness. Fitness will be multiplicative

# genotype at SA    male      female
#  00                1-s      1
#  01                1-.5s    1-.5s
#  11                1        1-s
#

# for all models females will recombine on all chromosomes
# for achiasmatic model males will have no recombination on
# any chromosomes but females will recombine as normal
# males that are chiasmatic will recombine only in the PAR 
# region on the sex chromosome.



# variables
# s: the selection coefficient for SA .5
# dfe: vector of selection coefficients for general fitness loci
# gen: number of generations to run
# iter: number of trials to use
# 

# 1 make random starting genomes

# 2 assess fitness

# 3 pick parents

# 4 make gametes

# 5 lay down mutations
dfe <- rgamma(5000, shape = .28, scale=113)
# values less than -1 are not meaningful
dfe[dfe>1] <- 1
dfe <- 1-dfe
hist(dfe)
# 6 make next gen

# 7 return to step 2


