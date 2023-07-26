# wget https://github.com/kwesterman/ESPRESSO.GxE.RV/raw/master/data/env.params.rda
# wget https://github.com/kwesterman/ESPRESSO.GxE.RV/raw/master/data/geno.params.rda
# wget https://github.com/kwesterman/ESPRESSO.GxE.RV/raw/master/data/pheno.params.rda
# wget https://github.com/kwesterman/ESPRESSO.GxE.RV/raw/master/data/simulation.params.rda
# 
# load(file='env.params.rda')
# load(file='geno.params.rda')
# load(file='pheno.params.rda')
# load(file='simulation.params.rda')

# devtools::install_github("https://github.com/kwesterman/ESPRESSO.GxE.RV")

# install.packages("remotes")
# remotes::install_github("agaye/ESPRESSO.GxE")

library(ESPRESSO.GxE.RV)
# library(ESPRESSO.GxE)

data(simulation.params)
data(pheno.params)
data(geno.params)
data(env.params)

str(simulation.params)
str(pheno.params)
str(geno.params)
str(env.params)

####### @@@@@
geno.params$geno.nvariants=1
####### @@@@@
geno.params$geno.frac_causal=0.1

simulation.params_test=simulation.params
simulation.params_test$seed.val[3]=1
simulation.params_test$numsubjects[3]=260417
simulation.params_test$p.val[3]=5e-8

####### @@@@@  TEST ####### @@@@@ 
simulation.params_test$interaction.efkt[3]=0.2

pheno.params_test=pheno.params
####### @@@@@
pheno.params_test$pheno.mean[3]=5.5
####### @@@@@
pheno.params_test$pheno.sd[3]=0.5
####### @@@@@
pheno.params_test$pheno.reliability[3]=1

geno.params_test=geno.params
####### @@@@@
geno.params_test$MAF[3]=0.0025
geno.params_test$geno.efkt[3]=0

env.params_test=env.params
####### @@@@@
env.params_test$env.model[3]=1
####### @@@@@
env.params_test$env.prevalence[3]=0.31
####### @@@@@
env.params_test$env.OR[3]=1.105171
env.params_test$env.efkt[3]=0.1
####### @@@@@
env.params_test$env.mean[3]=2.0
####### @@@@@
env.params_test$env.sd[3]=1.2
####### @@@@@
env.params_test$env.reliability[3]=1

# geno.params_test$MAF[3]=0.005
# env.params_test$env.efkt[3]=0.15

# Remove  Columns in List
simulation.params_test <- simulation.params_test[,!names(simulation.params_test) %in% c("numcases","numcontrols",
                                                                                        "interaction.OR","RR.5.95")]

pheno.params_test <- pheno.params_test[,!names(pheno.params_test) %in% c("disease.prev","pheno.sensitivity",
                                                                         "pheno.specificity")]

geno.params_test <- geno.params_test[,!names(geno.params_test) %in% c("geno.OR","geno.sensitivity",
                                                                      "geno.specificity")]

env.params_test <- env.params_test[,!names(env.params_test) %in% c("env.low.lim","env.up.lim",
                                                                   "env.sensitivity","env.specificity")]

run.espresso.GxE.RV(simulation.params_test, pheno.params_test, geno.params_test, env.params_test, scenarios2run=c(3))

str(simulation.params_test)
str(pheno.params_test)
str(geno.params_test)
str(env.params_test)
