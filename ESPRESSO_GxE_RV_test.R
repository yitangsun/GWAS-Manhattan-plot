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

# run the function for the first two scenarios, two binomial models
run.espresso.GxE(simulation.params, pheno.params, geno.params, env.params, scenarios2run=c(1,2))

# run the function for the last two scenarios, two gaussian models
run.espresso.GxE(simulation.params, pheno.params, geno.params, env.params, scenarios2run=c(3,4))

?ESPRESSO.GxE.RV::run.espresso.GxE()


run.espresso.GxE(simulation.params, pheno.params, geno.params, env.params, scenarios2run=c(4))

str(simulation.params)
str(pheno.params)
str(geno.params)
str(env.params)

geno.params_test=geno.params
geno.params_test$MAF[4]=0.005
env.params_test=env.params
env.params_test$env.efkt[4]=0.125

run.espresso.GxE(simulation.params, pheno.params, geno.params_test, env.params_test, scenarios2run=c(4))
