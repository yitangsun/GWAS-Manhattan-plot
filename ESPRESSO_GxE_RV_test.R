wget https://github.com/kwesterman/ESPRESSO.GxE.RV/raw/master/data/env.params.rda
wget https://github.com/kwesterman/ESPRESSO.GxE.RV/raw/master/data/geno.params.rda
wget https://github.com/kwesterman/ESPRESSO.GxE.RV/raw/master/data/pheno.params.rda
wget https://github.com/kwesterman/ESPRESSO.GxE.RV/raw/master/data/simulation.params.rda

load(file='env.params.rda')
load(file='geno.params.rda')
load(file='pheno.params.rda')
load(file='simulation.params.rda')


devtools::install_github("https://github.com/kwesterman/ESPRESSO.GxE.RV")

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

geno.params$geno.nvariants=1
geno.params$geno.frac_causal=0.1

geno.params_test=geno.params
geno.params_test$MAF[4]=0.005
env.params_test=env.params
env.params_test$env.efkt[4]=0.15

geno.params_test$MAF[3]=0.005
env.params_test$env.efkt[3]=0.15


run.espresso.GxE.RV(simulation.params, pheno.params, geno.params_test, env.params_test, scenarios2run=c(3,4))

str(simulation.params)
str(pheno.params)
str(geno.params_test)
str(env.params_test)





# run the function for the first two scenarios, two binomial models
run.espresso.GxE(simulation.params, pheno.params, geno.params, env.params, scenarios2run=c(1,2))

# run the function for the last two scenarios, two gaussian models
run.espresso.GxE(simulation.params, pheno.params, geno.params, env.params, scenarios2run=c(3,4))

run.espresso.GxE(simulation.params, pheno.params, geno.params, env.params, scenarios2run=c(4))


# run the function for the first two scenarios, two binomial models
run.espresso.GxE.RV(simulation.params, pheno.params, geno.params, env.params, scenarios2run=c(1,2))

# run the function for the last two scenarios, two gaussian models
run.espresso.GxE.RV(simulation.params, pheno.params, geno.params, env.params, scenarios2run=c(3,4))

run.espresso.GxE.RV(simulation.params, pheno.params, geno.params, env.params, scenarios2run=c(4))


run.espresso.GxE.RV(simulation.params=NULL, pheno.params=NULL, geno.params=NULL, 
                                env.params=NULL, scenarios2run=1)

run.espresso.GxE(simulation.params=NULL, pheno.params=NULL, geno.params=NULL, 
                    env.params=NULL, scenarios2run=1)

?ESPRESSO.GxE.RV::run.espresso.GxE()
?ESPRESSO.GxE.RV::get.critical.results.GxE

run.espresso.GxE(simulation.params, pheno.params, geno.params, env.params, scenarios2run=c(4))

str(simulation.params)
str(pheno.params)
str(geno.params)
str(env.params)

geno.params_test=geno.params
geno.params_test$MAF[4]=0.005
env.params_test=env.params
env.params_test$env.efkt[4]=0.125

