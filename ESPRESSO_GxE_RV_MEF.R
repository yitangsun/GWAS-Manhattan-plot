# sim.QTL.data.GxE
#'
#' @title Simulates subjects for continuous outcome
#' @description Generates the specified number of subjects
#' @param n Number of subjects to simulate
#' @param ph.mean statistical mean
#' @param ph.sd standard deviation
#' @param g.nvar Integer number of genetic variants to simulate
#' @param g.frac_causal Fraction of causal genetic variants
#' @param freq Minor allele frequencies of the genetic variants
#' @param g.model Genetic model; 0 for binary and 1 for continuous
#' @param g.efkt Effects of the genetic variant
#' @param e.model Model of the environmental exposure
#' @param e.efkt Effects of the environment determinats
#' @param e.prev Prevalences of the environmental determinants
#' @param e.mean Mean under quantitative-normal model
#' @param e.sd Standard deviation under quantitative-normal model
#' @param e.low.lim lower limit under quantitative-uniform model
#' @param e.up.lim upper limit under quantitative-uniform model
#' @param i.efkt Interaction effect
#' @param pheno.rel reliability of the assessment for a quantitative outcome.
#' @return A matrix
#' @keywords internal
#' @author Gaye A.; Westerman K.

sim.QTL.data.GxE <- function(
    n=NULL,ph.mean=NULL,ph.sd=NULL,
    g.nvar=NULL,g.frac_causal=NULL,freq=NULL,g.model=NULL,g.efkt=NULL,
    e.model=NULL,e.efkt=NULL,e.prev=NULL,e.mean=NULL,e.sd=NULL,e.low.lim=NULL,
    e.up.lim=NULL,i.efkt=NULL,pheno.rel=NULL
) {
  
  # GENERATE THE TRUE GENOTYPE DATA
  geno <- matrix(nrow=n, ncol=g.nvar)
  for (j in 1:g.nvar) {
    geno.data <- sim.geno.data(num.obs=n, geno.model=g.model, MAF=freq)
    allele.A <- geno.data$allele.A  # Is this info ultimately needed? Can it be the same across all variants?
    allele.B <- geno.data$allele.B
    geno[, j] <- geno.data$genotype
  }
  
  # GENERATE THE TRUE ENVIRONMENTAL EXPOSURE DATA			
  env <- sim.env.data(num.obs=n,env.model=e.model,env.prev=e.prev,env.mean=e.mean,
                      env.sd=e.sd,env.low.lim=e.low.lim,env.up.lim=e.up.lim)
  
  # GENERATE THE TRUE INTERACTION DATA           
  interaction <- geno * env  # N x g.nvar matrix
  
  # GENERATE THE TRUE OUTCOME DATA
  pheno.data <- sim.pheno.qtl.GxE(numsubjects=n,pheno.mean=ph.mean,pheno.sd=ph.sd,
                                  genotype=geno,geno.efkt=g.efkt,environment=env,env.efkt=e.efkt,
                                  interaction=interaction,int.efkt=i.efkt,
                                  geno.frac_causal=g.frac_causal)
  true.phenotype <- pheno.data
  
  # GENERATE THE OBSERVED OUTCOME DATA 
  obs.phenotype <- get.obs.pheno(phenotype=true.phenotype, pheno.model=1, 
                                 pheno.sd=ph.sd, pheno.reliability=pheno.rel)
  pheno <- obs.phenotype
  
  # STORE THE GENERATED TRUE DATA INTO AN OUTPUT MATRIX 
  sim.matrix <- cbind(pheno,env,geno,interaction,allele.A,allele.B)
  
  # ADD IDs (JUST A ROW COUNT)
  sim.matrix <- cbind(1:nrow(sim.matrix), sim.matrix)
  
  # ADD COLUMN NAMES AND RETURN A DATAFRAME
  colnames(sim.matrix) <- c("id","phenotype","environment",
                            paste0("genotype", 1:ncol(geno)),
                            paste0("interaction", 1:ncol(interaction)),
                            "allele.A","allele.B")
  mm <- data.frame(sim.matrix)
}

# sim.geno.data
#' 
#' @title Generates genotypes for a genetic variant
#' @description Generates two alleles and combines them to form the genotype of 
#' a SNP under a binary or additive genetic model.
#' @param num.obs number of observations to generate.
#' @param geno.model genetic model; binary=0 and additive=1.
#' @param MAF minor allele frequency of the SNP (in ESPRESSO this is the frequency of the 'at risk' allele)
#' @return a dataframe that contains the following variables:
#' \code{allele.A} major allele
#' \code{allele.B} minor allele
#' \code{genotype} genotype
#' @keywords internal
#' @author Gaye A.
#' 
sim.geno.data <- function(num.obs=NULL, geno.model=NULL, MAF=NULL){
  
  
  # CORRECTION TERM FOR MEAN CENTERING FOR ADDITIVE 
  mean.geno.additive <- (2*MAF*(1-MAF)+2*(MAF^2))
  
  # CORRECTION TERM FOR MEAN CENTERING FOR BINARY GENE
  mean.geno.binary <- (2*MAF*(1-MAF)+(MAF^2))
  
  # CREATE, CENTRE AND ROUND AN ADDITIVE GENETIC GENOTYPE COVARIATE WITH
  # APPROPRIATE MAF 
  allele.A <- rbinom(num.obs,1,MAF)
  allele.B <- rbinom(num.obs,1,MAF)
  
  # ACTUAL GENOTYPE IS SUM OF ALLELES IN ADDITIVE GENETIC MODEL
  # AND IS 1 IF SUM OF ALLELES IS 1 OR GREATER IN THE BINARY MODEL 
  genotype <- allele.A+allele.B
  
  if(geno.model==0){
    geno.U <- genotype > 0
    genotype <- geno.U - mean.geno.binary
  }
  if(geno.model==1){
    geno.U <- genotype
    genotype <- geno.U - mean.geno.additive
  }
  
  # return the data as a dataframe
  dtframe <- data.frame(allele.A, allele.B, genotype)
}

# sim.env.data
#' 
#' @title Simulates cases and controls
#' @description Generates data for a binary, quantitative-normal, or quantitative-uniform environmental determinant.
#' @param num.obs number of observations to simulate.
#' @param env.model model of the exposure: binary=0, quantitative-normal=1, quantitative-uniform=2.
#' @param env.prev prevalence of the environmental exposure.
#' @param env.mean statisitical man under quantitative-normal model.
#' @param env.sd standard deviation under quantitative-normal model.
#' @param env.low.lim lower limit under quantitative-uniform model.
#' @param env.up.lim upper limit under quantitative-uniform model.
#' @return a vector of continuous or binary values.
#' @keywords internal
#' @author Gaye A.
#' 
sim.env.data <- function(num.obs=NULL, env.model=NULL, env.prev=NULL, env.mean=NULL,
                         env.sd=NULL,env.low.lim=NULL,env.up.lim=NULL){
  
  # CREATE THE FIRST ENVIRONMENTAL COVARIATE 
  if(env.model==0){   # BINARY DISTRIBUTION
    env.U <- rbinom(num.obs, 1, env.prev)
    e.mean <- env.prev
    env.U <- env.U-e.mean
  }
  if(env.model==1){   # NORMAL DISTRIBUTION
    env.U <- rnorm(num.obs, env.mean, env.sd)
    env.U <- env.U-mean(env.U)     # mean centering
  }
  if(env.model==2){  # UNIFORM DISTRIBUTION
    if(env.low.lim >= env.up.lim){
      stop("\n\nALERT!\n Uniform Distribution: The upper limit must be greater than the lower limit\n\n")
    }else{
      env.U <- runif(num.obs, env.low.lim, env.up.lim)
      env.U <- env.U-mean(env.U)     # mean centering
    }
  }
  
  # return a vector 
  return(env.U)
}

# sim.pheno.qtl.GxE
#' 
#' @title Simulates continuous outcome data
#' @description Uses the effects data of the determinants to construct a linear predictor(LP). The outcome is normally distributed variable generated with a mean equal to LP and a standard deviation of 1. Some error is then added to the simulated outcome to obtained the observed outcome.
#' @param num.subjects Number of subjects to simulate
#' @param pheno.mean statistical mean
#' @param pheno.sd standard deviation
#' @param genotype Genotype matrix (N x g.nvar)
#' @param geno.efkt Effect of the genetic variant
#' @param geno.frac_causal Fraction of causal genetic variants
#' @param environment Exposure data for environment
#' @param env.efkt Effect of the environmental determinants
#' @param interaction Interaction matrix (N x g.nvar)
#' @param int.efkt Interaction effect
#' @return A dataframe of phenotype
#' @keywords internal
#' @author Gaye A.; Westerman K.
#'
sim.pheno.qtl.GxE <- function(
    numsubjects=NULL,pheno.mean=NULL,pheno.sd=NULL,
    genotype=NULL,geno.efkt=NULL,geno.frac_causal=NULL,
    environment=NULL,env.efkt=NULL,interaction=NULL,int.efkt=NULL
) {  
  
  # ALPHA IS EQUAL TO THE MEAN OF THE TRAIT, WHICH IS 0
  num.obs <- numsubjects
  alpha <- pheno.mean
  
  # SET FINAL VARIANT EFFECTS BASED ON THE CAUSAL FRACTION
  num_causal <- ceiling(geno.frac_causal * ncol(genotype))  # Round up
  geno.efkts <- c(rep(geno.efkt, num_causal), 
                  rep(0, ncol(genotype) - num_causal))
  int.efkts <- c(rep(int.efkt, num_causal), 
                 rep(0, ncol(genotype) - num_causal))
  
  # GENERATE THE LINEAR PREDICTOR
  lp <- alpha + 
    (genotype %*% geno.efkts) + 
    (environment * env.efkt) + 
    (interaction %*% int.efkts)
  
  # GENERATE THE TRUE PHENOTYPE DATA
  phenotype <- rnorm(num.obs, mean=lp, sd=sqrt(pheno.sd^2 - var(lp)))
  
  return(phenotype)
}

# get.obs.pheno
#' 
#' @title Generates observed outcome data
#' @description Adds a set level of error to error free binary or quantitative data (the true phenotype data) 
#' to obtain data with a larger variance (the observed phenotype data).
#' @param seed 
#' @param phenotype outcome status.
#' @param pheno.model distribution of the outcome variable: binary=0, normal=1
#' @param pheno.sd standard deviation of the outcome in the study population
#' @param pheno.model distribution of the outcome variable: binary=0, normal=1 or uniform=2.
#' @param pheno.sd standard deviation of the outcome in the study the population
#' @param pheno.error misclassification rates: 1-sensitivity and 1-specificity
#' @param pheno.reliability reliability of the assessment for a quantitative outcome.
#' @return A dataframe containing:
#' \code{true.phenotype} the error free outcome data (true data).
#' \code{observed.phenotype} the true outcome data with some added error (observed data).
#' @keywords internal
#' @author Gaye A.
#'
get.obs.pheno <- function (phenotype=NULL, pheno.model=NULL, pheno.sd=NULL, pheno.error=NULL, pheno.reliability=NULL){ 
  
  # GET THE OBSERVED OUTCOME DATA
  if(pheno.model==0){ # IF THE OUTCOME IS BINARY
    
    pheno.original<-phenotype
    disease.missed<-rbinom(length(phenotype),1,pheno.error[1])
    non.disease.missed<-rbinom(length(phenotype),1,pheno.error[2])
    pheno.U<-
      ((pheno.original==0) * (disease.missed==0) * (non.disease.missed==0))*pheno.original+
      ((pheno.original==1) * (disease.missed==0) * (non.disease.missed==0))*pheno.original+
      ((pheno.original==0) * (disease.missed==0) * (non.disease.missed==1))*(1-pheno.original)+
      ((pheno.original==1) * (disease.missed==0) * (non.disease.missed==1))*pheno.original+
      ((pheno.original==0) * (disease.missed==1) * (non.disease.missed==0))*pheno.original+
      ((pheno.original==1) * (disease.missed==1) * (non.disease.missed==0))*(1-pheno.original)+
      ((pheno.original==0) * (disease.missed==1) * (non.disease.missed==1))*(1-pheno.original)+
      ((pheno.original==1) * (disease.missed==1) * (non.disease.missed==1))*(1-pheno.original)
    observed.phenotype <- pheno.U
    
  }else{ # IF THE OUTCOME IS CONTINUOUS NORMAL
    # USE THE RELIABITLITY OF PHENOTYPE ASSESSMENT TO COMPUTE THE VARIANCE OF MEASURED PHENOTYPES.
    # RELIABITY = (VAR.true.data/VAR.obs.data) AND VAR.obs.data = (VAR.true.data + VAR.measurement)
    # IT FOLLOWS THAT VAR.true.data + VAR.of.estimate = VAR.true.data/RELIABILITY AND THEN:
    # VAR.measurement = (VAR.true.data/RELIABILITY) - VAR.true.data
    var.m <- (pheno.sd^2/pheno.reliability)-(pheno.sd^2)
    
    # ADD THE ERROR TO ORIGINAL PHENOTYPES TO GENERATE THE OBSERVED PHENOTYPE DATA
    n <- length(phenotype)
    observed.phenotype <- rnorm(n, phenotype, sqrt(var.m))
  } 
  
  # RETURN THE TRUE AND OBSERVED PHENOTYPE DATA AS A DATAFRAME
  df <- observed.phenotype
  
}

# get.observed.data.GxE
#'
#' @title Generates exposure data with some error
#' @description Uses functions make.obs.geno and make.obs.env to generate effect data with a set level of error
#' @param data Input table of simulated data considered as true data
#' @param g.error Misclassification rates in genetic assessment: 1-sensitivity and 1-specificity
#' @param g.model Genetic model; 0 for binary and 1 for additive
#' @param freq Minor allele frequency
#' @param e.error Misclassification rates in environmental exposures assessment: 1-sensitivity and 1-specificity
#' @param e.model Model of the exposure: binary=0, quantitative-normal=1 or quantitative-uniform=2
#' @param e.prev Prevalence of environmental exposure
#' @param e.sd Standard Deviation
#' @param e.reliability Reliability of the assessment of quantitative exposure
#' @return A matrix
#' @keywords internal
#' @author Gaye A.; Westerman K.
#'
get.observed.data.GxE <- function(
    data=NULL,g.error=NULL,g.model=NULL,freq=NULL,
    e.error=NULL,e.model=NULL,e.prev=NULL,e.sd=NULL,
    e.reliability=NULL
) {
  
  sim.df <- data
  
  # # GET THE OBSERVED GENOTYPES
  # # KW comment: ESPRESSO.GxE.RV won't add any genotyping error for now
  # for (genotype in grep("genotype", names(sim.df), value=TRUE)) {
  #     true.genotype <- sim.df[[genotype]]
  #     obs.genotype <- get.obs.geno(allele.A=sim.df$allele.A,allele.B=sim.df$allele.B,
  #                                  geno.model=g.model,MAF=freq,geno.error=g.error)
  # }
  obs.genotype <- sim.df[, grepl("genotype", names(sim.df)), drop=FALSE]
  
  # GET THE OBSERVED ENVIRONMENTAL EXPOSURE DATA
  true.environment <- sim.df$environment
  
  obs.environment <- get.obs.env(env.data=true.environment,env.model=e.model,env.prev=e.prev,
                                 env.sd=e.sd,env.error=e.error,env.reliability=e.reliability)
  
  # GET THE OBSERVED INTERACTION DATA
  obs.interaction <- obs.genotype * obs.environment
  colnames(obs.interaction) <- gsub("genotype", "interaction",
                                    colnames(obs.interaction))
  
  # REPLACE THE TRUE DATA BY THE NOW GENERATED OBSERVED GENOTYPES
  # IN THE INITIAL MATRIX THAT HELD THE TRUE DATA
  # sim.df$genotype <- obs.genotype$observed.genotype
  # sim.df$allele.A <- obs.genotype$observed.allele.A
  # sim.df$allele.B <- obs.genotype$observed.allele.B
  sim.df$environment <- obs.environment
  for (field in names(obs.interaction)) {
    sim.df[[field]] <- obs.interaction[[field]]
  }
  # sim.df$interaction <- obs.interaction
  
  # RETURN THE MATRIX WHICH NOW CONTAINS ONLY THE OBSERVED DATA TO ANALYSE BY GLM
  # colnames(sim.df) <- c("id", "phenotype", "genotype", "allele.A", "allele.B", "environment", "interaction")
  return(sim.df)
}

# get.obs.env
#'
#' @title Generates the observed environmental exposure data
#' @description Adds a set level of error to binary or quantitative data (the true data)
#' to obtain data with a larger variance (the observed data). The level of error is determined by
#' the misclassification rates in binary exposure and by the set level of variance in the quantitative
#' exposure.
#' @param env.data a vector of environmental measures that represents the true data.
#' @param env.model distribution of the exposure: binary=0 , normal=1 or uniform=2.
#' @param env.prev prevalence of the environmental exposure.
#' @param env.error misclassification rates: 1-sensitivity and 1-specificity.
#' @param env.reliability reliability of the assessment of quantitative exposure.
#' @return a dataframe with two coloumns/variables:
#' \code{true.environment} the error free exposure data (true data).
#' \code{observed.environment} the true esposure data with some added error (observed data).
#' @keywords internal
#' @author Amadou Gaye
#' 
get.obs.env <- function(env.data=NULL, env.model=NULL, env.sd=NULL, env.prev=NULL, env.error=NULL, env.reliability=NULL){
  
  if(env.model==0){
    numsubs <- length(env.data)
    environ.1.missed<-rbinom(numsubs,1,env.error[1])
    environ.0.missed<-rbinom(numsubs,1,env.error[2])
    
    environ <- env.data
    environ.test <- (environ > 0) * 1
    environ.new <- 
      ((environ.test==0) * (environ.1.missed==0) * (environ.0.missed==0) * environ.test)+
      ((environ.test==1) * (environ.1.missed==0) * (environ.0.missed==0) * environ.test)+
      ((environ.test==0) * (environ.1.missed==0) * (environ.0.missed==1) * (1-environ.test))+
      ((environ.test==1) * (environ.1.missed==0) * (environ.0.missed==1) * environ.test)+
      ((environ.test==0) * (environ.1.missed==1) * (environ.0.missed==0) * environ.test)+
      ((environ.test==1) * (environ.1.missed==1) * (environ.0.missed==0) * (1-environ.test))+
      ((environ.test==0) * (environ.1.missed==1) * (environ.0.missed==1) * (1-environ.test))+
      ((environ.test==1) * (environ.1.missed==1) * (environ.0.missed==1) * (1-environ.test))
    obs.env <- environ.new-env.prev
    
  }else{
    var.error <- (env.sd^2/env.reliability)-env.sd^2
    numsubs <- length(env.data)
    if(env.model==1){
      obs.env <- rnorm(numsubs, env.data, sqrt(var.error))
    }else{
      env.model.error <- rnorm(numsubs, 0, sqrt(var.error))
      obs.env <- env.data + env.model.error  
    }
  }
  
  return(obs.env)
}

# glm.analysis.GxE
#' 
#' @title Carries out regression analysis
#' @description Fits a conventional unconditional logistic regression model with a binary or continuous phentype as outcome and the genetic, environmental, interaction determinants as covariates.
#' @param pheno.model Type of outcome; 0=binary and 1=continuous
#' @param observed.data A dataframe that contains covariates and outcome data
#' @return A vector containing the beta, standard-error, z-statistic, and p-value for each of the covariates
#' @keywords internal
#' @author Gaye A.; Westerman K.
#'
glm.analysis.GxE <- function(pheno.model=NULL, observed.data=NULL, g.idx=NULL) {
  
  observed.data$genotype <- observed.data[[paste0("genotype", g.idx)]]
  observed.data$interaction <- observed.data[[paste0("interaction", g.idx)]]
  
  # BINARY OUTCOME
  if(pheno.model == 0){
    # FIT CONVENTIONAL UNCONDITIONAL LOGISTIC REGRESSION MODEL
    mod.glm <- glm(phenotype ~ 1+genotype+environment+interaction,
                   family=binomial,data=observed.data)
    mod.sum <- summary(mod.glm)
  }
  
  # QUANTITATIVE OUTCOME
  if(pheno.model == 1){
    # FIT A GLM FOR A GAUSSIAN OUTCOME
    mod.glm <- glm(phenotype ~ 1+genotype+environment+interaction,
                   family=gaussian,data=observed.data)
    mod.sum <- summary(mod.glm)     
  }
  
  beta.value <- mod.sum$coefficients[4,1]
  se.value <- mod.sum$coefficients[4,2]
  z.value <- mod.sum$coefficients[4,3]
  p.value <- mod.sum$coefficients[4,4]
  
  # RETURN A VECTOR
  return(list(beta=beta.value, se=se.value, z=z.value, p=p.value))
}

combine.variants <- function(variant_results_list) {
  # This implements Fisher's method for p-value combination
  pvec <- sapply(variant_results_list, function(x) x[["p"]])
  p <- pchisq(-2 * sum(log(pvec)), df=2 * length(pvec), lower.tail=FALSE)
  p
}

# power.calc
#' 
#' @title Calculates the empirical and theoretical power
#' @description The function determines the empirical and theoretical power. 
#' The empirical power is the proportion of simulations in which 
#' the z-statistic for the parameter of interest exceeds the z-statistic 
#' for the desured level if statistical significance. 
#' The theoretical power is the power of the study.
#' @param pval.thresh cut-off p-value defining statistical significance.
#' @param p.values Vector of p-values (one per genetic variant).
#' @return a list that contains the computed empirical power and theoretical power.
#' @keywords internal
#' @author Gaye A.; Westerman K.
#'
power.calc <- function(pval.thresh=NULL, p.values=NULL){
  
  #if(is.null(z.values)){
  #  cat("\n\n ALERT!\n")
  #  cat(" No z-statistics found\n")
  #  cat(" Check the argument 'z.values'.\n")
  #  stop(" End of process!\n\n", call.=FALSE)
  #}
  #
  #if(is.null(mean.model.z)){
  #  cat("\n\n ALERT!\n")
  #  cat(" The argument 'mean.model.z' is set to NULL.\n")
  #  cat(" This argument should be the ratio 'mean.beta/mean.se'.\n")
  #  stop(" End of process!\n\n", call.=FALSE)
  #}
  
  ## CALCULATE Z STATISTIC THRESHOLD FOR DESIRED P-VALUE 
  #z.pval <- qnorm(1-pval.thresh/2)
  
  # GET EMPIRICAL POWER: THE PROPORTION OF SIMULATIONS IN WHICH THE 
  # Z STATISTIC FOR THE PARAMETER OF INTEREST EXCEEDS THE Z STATISTIC 
  # FOR THE DESIRED LEVEL OF STATISTICAL SIGNIFICANCE
  #empirical.power <- round(mean((z.values > z.pval), na.rm=TRUE),3)
  empirical.power <- round(mean((p.values < pval.thresh), na.rm=TRUE), 3)
  
  # # GET THE MODELLED POWER
  # modelled.power <- 1 #pnorm(mean.model.z-z.pval)
  
  return(empirical.power)  
  # return(list(empirical=empirical.power, modelled=modelled.power))
}

# get.critical.results.GxE
#'
#' @title Summarizes the main results
#' @description Gets the number of cases and controls or subjects and the empirical and theoretical power under each model and prints a summary on the screen
#' @param scenario Scenario number
#' @param pheno.model Type of the outcome; 0 for binary and 2 for continuous
#' @param geno.model Genetic model; 0 for binary and 1 for additive
#' @param env.model Model of the enviromental explosure
#' @param empirical.power Estimated empirical power
##' @param sample.sizes.required Number of cases and controls or number of subjects required to achieve the desired power
#' @return A table containing the following items:
#' \code{genetic.model} Model of the genetic determinant
#' \code{environment.model} Model of the environmental determinant
##' \code{number.of.subjects.required} Number of subjects required to achieve the desired power under a quantitative outcome model
#' \code{empirical.power} Estimated empirical power under each model
#' @keywords internal
#' @author Gaye A.; Westerman K.
#'
get.critical.results.GxE <- function(
    scenario=NULL, pheno.model=NULL,geno.model=NULL,env.model=NULL,
    empirical.power=NULL#,sample.sizes.required=NULL
    # modelled.power=NULL,mean.beta=NULL)
) {
  
  if(geno.model==0){
    g.model <- "binary"
  }else{
    g.model <- "additive"
  }
  
  if(env.model==0){
    e.model <- "binary"
  }else{
    if(env.model==1){
      e.model <- "normally distributed"
    }else{
      e.model <- "uniformly distributed"
    }
  }
  #if(pheno.model == 0){
  #    				numcases <- sample.sizes.required[[1]]
  #    			 numcontrols <- sample.sizes.required[[2]]
  #}else{
  #    				numsubjects <- sample.sizes.required[[1]]
  #}
  
  if(pheno.model==0) {
    # estimated ORs
    # estimated.OR <- exp(mean.beta)
    # estimated.effect <- 'NA'
    
    cat("\n---- SUMMARY OF SCENARIO",scenario,"----\n")
    cat("\nModels\n")
    cat("------\n")
    cat("  Outcome: binary \n")
    cat("  Genetic determinant:",g.model,"\n")
    cat("  Environmental determinant:",e.model)
    
    #cat("\n\nNumber of cases required\n")
    #cat("------------------------\n")
    #cat(" ", numcases)
    #
    #cat("\n\nNumber of controls required\n")
    #cat("---------------------------\n")
    #cat(" ", numcontrols)
    
    cat("\n\nEmpirical power\n")
    cat("---------------\n")
    cat(" ",empirical.power)
    
    # cat("\n\nModel power\n")
    # cat("-----------\n")
    # cat(" ",round(modelled.power,2))
    
    # cat("\n\nEstimated interaction OR\n")
    # cat("-----------\n")
    # cat(" ",round(estimated.OR,2))
    
    cat("\n\n---- END OF SUMMARY ----\n")
    
    # crit.res <- c(g.model,e.model,numcases,numcontrols,round(empirical.power,2),round(modelled.power,2),
    #               round(estimated.OR,2), estimated.effect)
    # return(list(genetic.model=crit.res[1], environment.model=crit.res[2],number.of.cases.required=crit.res[3],
    #             number.of.controls.required=crit.res[4],empirical.power=crit.res[5], 
    #             modelled.power=crit.res[6],estimated.OR=crit.res[7], estimated.effect=crit.res[8]))
    return(list(genetic.model=g.model, environment.model=e.model, #numer.of.subjects.required=numcases+numcontrols,
                empirical.power=round(empirical.power,2)))
    
  } else {
    # estimated ORs
    # estimated.effect <- mean.beta
    # estimated.OR <- 'NA'
    
    cat("\n---- SUMMARY OF SCENARIO",scenario,"----\n")
    cat("\nModels\n")
    cat("------\n")
    cat(" Outcome: quantitative; ")
    cat(" Genetic determinant: ",g.model,"\n")
    cat(" Environmental determinant: ",e.model)
    
    #cat("\n\nNumber of subjects required\n")
    #cat("------------------------\n")
    #cat(" ",numsubjects)
    
    cat("\n\nEmpirical power\n")
    cat("---------------\n")
    cat(" ",empirical.power)
    
    # cat("\n\nModel power\n")
    # cat("-----------\n")
    # cat(" ",round(modelled.power,2))
    # 
    # cat("\n\nEstimated interaction effect\n")
    # cat("-----------\n")
    # cat(" ",estimated.effect)
    
    cat("\n\n---- END OF SUMMARY ----\n")
    
    # crit.res <- c(g.model,e.model,numsubjects,round(empirical.power,2),round(modelled.power,2),estimated.OR,round(estimated.effect,2))
    # return(list(genetic.model=crit.res[1], environment.model=crit.res[2],number.of.subjects.required=crit.res[3],
    #            empirical.power=crit.res[4], modelled.power=crit.res[5], estimated.OR=crit.res[6],estimated.effect=crit.res[7]))
    return(list(genetic.model=g.model, environment.model=e.model, #number.of.subjects.required=numsubjects,
                empirical.power=round(empirical.power,2)))
  }
}

################## yitang_run.espresso.GxE.RV
################## yitang_run.espresso.GxE.RV
################## yitang_run.espresso.GxE.RV

yitang_run.espresso.GxE.RV <- function (simulation.params = NULL, pheno.params = NULL, geno.params = NULL, 
          env.params = NULL, scenarios2run = 1) 
{
  if (is.null(simulation.params)) {
    cat("\n WARNING!\n")
    cat(" No simulation parameters supplied\n")
    cat(" The default simulation parameters will be used\n")
    simulation.params <- data(simulation.params)
  }
  if (is.null(pheno.params)) {
    cat("\n WARNING!\n")
    cat(" No outcome parameters supplied\n")
    cat(" The default outcome parameters will be used\n")
    pheno.params <- data(pheno.params)
  }
  if (is.null(geno.params)) {
    cat("\n WARNING!\n")
    cat(" No genotype parameters supplied\n")
    cat(" The default genotype parameters will be used\n")
    geno.params <- data(geno.params)
  }
  if (is.null(env.params)) {
    cat("\n WARNING!\n")
    cat(" No environmental parameters supplied\n")
    cat(" The default environmental parameters will be used\n")
    env.params <- data(env.params)
  }
  s.temp1 <- merge(simulation.params, pheno.params)
  s.temp2 <- merge(s.temp1, geno.params)
  s.parameters <- merge(s.temp2, env.params)
  trace.interval <- 10
  allowed.sample.size <- 2e+07
  block.size <- 20000
  output.file <- "output.RV.csv"
  output.matrix <- matrix(numeric(0), ncol = ncol(s.parameters) + 
                            2)
  column.names <- c(colnames(s.parameters), "exceeded.sample.size?", 
                    "empirical.power")
  write(t(column.names), output.file, dim(output.matrix)[2], 
        sep = ";")
  for (j in scenarios2run) {
    set.seed(s.parameters$seed.val[j])
    scenario.id <- s.parameters$scenario.id[j]
    seed.val <- s.parameters$seed.val[j]
    numsims <- s.parameters$numsims[j]
    numcases <- s.parameters$numcases[j]
    numcontrols <- s.parameters$numcontrols[j]
    numsubjects <- s.parameters$numsubjects[j]
    int.OR <- s.parameters$interaction.OR[j]
    int.efkt <- s.parameters$interaction.efkt[j]
    baseline.OR <- s.parameters$RR.5.95[j]
    pval <- s.parameters$p.val[j]
    power <- s.parameters$power[j]
    pheno.model <- s.parameters$pheno.model[j]
    pheno.mean <- s.parameters$pheno.mean[j]
    pheno.sd <- s.parameters$pheno.sd[j]
    disease.prev <- s.parameters$disease.prev[j]
    pheno.error <- c(1 - s.parameters$pheno.sensitivity[j], 
                     1 - s.parameters$pheno.specificity[j])
    pheno.reliability <- s.parameters$pheno.reliability[j]
    geno.model <- s.parameters$geno.model[j]
    MAF <- s.parameters$MAF[j]
    geno.OR <- s.parameters$geno.OR[j]
    geno.efkt <- s.parameters$geno.efkt[j]
    geno.error <- c(1 - s.parameters$geno.sensitivity[j], 
                    1 - s.parameters$geno.specificity[j])
    geno.nvariants <- s.parameters$geno.nvariants[j]
    geno.frac_causal <- s.parameters$geno.frac_causal[j]
    env.model <- s.parameters$env.model[j]
    env.prev <- s.parameters$env.prevalence[j]
    env.OR <- s.parameters$env.OR[j]
    env.efkt <- s.parameters$env.efkt[j]
    env.mean <- s.parameters$env.mean[j]
    env.sd <- s.parameters$env.sd[j]
    env.low.lim <- s.parameters$env.low.lim[j]
    env.up.lim <- s.parameters$env.up.lim[j]
    env.error <- c(1 - s.parameters$env.sensitivity[j], 1 - 
                     s.parameters$env.specificity[j])
    env.reliability <- s.parameters$env.reliability[j]
    p.values <- rep(NA, numsims)
    sample.size.excess <- 0
    for (s in 1:numsims) {
      if (pheno.model == 0) {
        sim.data <- sim.CC.data.GxE(n = block.size, ncases = numcases, 
                                    ncontrols = numcontrols, max.sample.size = allowed.sample.size, 
                                    pheno.prev = disease.prev, g.nvar = geno.nvariants, 
                                    g.frac_causal = geno.frac_causal, freq = MAF, 
                                    g.model = geno.model, g.OR = geno.OR, e.model = env.model, 
                                    e.prev = env.prev, e.mean = env.mean, e.sd = env.sd, 
                                    e.low.lim = env.low.lim, e.up.lim = env.up.lim, 
                                    e.OR = env.OR, i.OR = int.OR, b.OR = baseline.OR, 
                                    ph.error = pheno.error)
        true.data <- sim.data$data
      }
      else {
        true.data <- sim.QTL.data.GxE(n = numsubjects, 
                                      ph.mean = pheno.mean, ph.sd = pheno.sd, g.nvar = geno.nvariants, 
                                      g.frac_causal = geno.frac_causal, freq = MAF, 
                                      g.model = geno.model, g.efkt = geno.efkt, e.model = env.model, 
                                      e.efkt = env.efkt, e.prev = env.prev, e.mean = env.mean, 
                                      e.sd = env.sd, e.low.lim = env.low.lim, e.up.lim = env.up.lim, 
                                      i.efkt = int.efkt, pheno.rel = pheno.reliability)
      }
      observed.data <- get.observed.data.GxE(data = true.data, 
                                             g.error = geno.error, g.model = geno.model, freq = MAF, 
                                             e.error = env.error, e.model = env.model, e.prev = env.prev, 
                                             e.sd = env.sd, e.reliability = env.reliability)
      variant_results_list <- list()
      for (g.idx in 1:geno.nvariants) {
        glm.estimates.variant <- glm.analysis.GxE(pheno.model, 
                                                  observed.data, g.idx)
        variant_results_list[[g.idx]] <- glm.estimates.variant
      }
      p.combined <- combine.variants(variant_results_list)
      p.values[s] <- p.combined
      if (s%%trace.interval == 0) 
        cat("\n", s, "of", numsims, "runs completed in scenario", 
            scenario.id)
    }
    cat("\n\n")
    empirical.power <- power.calc(pval, p.values)
    critical.res <- get.critical.results.GxE(j, pheno.model, 
                                             geno.model, env.model, empirical.power)
    if (pheno.model == 0) {
      sample.size.excess <- sim.data$allowed.sample.size.exceeded
      if (sample.size.excess == 1) {
        excess <- "yes"
        cat("\nTO GENERATE THE NUMBER OF CASES SPECIFIED AT OUTSET\n")
        cat("THE SIMULATION EXCEEDED THE MAXIMUM POPULATION SIZE OF ", 
            allowed.sample.size, "\n")
      }
      else {
        excess <- "no"
      }
    }
    inparams <- s.parameters[j, ]
    cat(colnames(inparams))
    if (pheno.model == 0) {
      if (env.model == 0) {
        irrelevant_params <- c("geno.efkt", "env.efkt", 
                               "interaction.efkt", "env.mean", "env.sd", "env.low.lim", 
                               "env.up.lim", "env.reliability")
        inparams[irrelevant_params] <- "NA"
        inputs <- inparams
      }
      else if (env.model == 1) {
        irrelevant_params <- c("geno.efkt", "env.efkt", 
                               "interaction.efkt", "env.prevalence", "env.OR", 
                               "env.sensitivity", "env.specificity")
        inputs <- inparams
      }
      outputs <- c(sample.size.excess, critical.res$empirical.power)
    }
    else {
      mod <- "quantitative"
      if (env.model == 0) {
        irrelevant_params <- c("geno.OR", "env.OR", "interaction.OR", 
                               "env.mean", "env.sd", "env.low.lim", "env.up.lim", 
                               "env.reliability")
        inputs <- inparams
      }
      else if (env.model == 1) {
        irrelevant_params <- c("geno.OR", "env.OR", "interaction.OR", 
                               "env.prevalence", "env.OR", "env.sensitivity", 
                               "env.specificity")
        inputs <- inparams
      }
      outputs <- c(sample.size.excess, critical.res$empirical.power)
    }
    jth.row <- as.numeric(as.character(c(inputs, outputs)))
    write(t(jth.row), output.file, dim(output.matrix)[2], 
          append = TRUE, sep = ";")
    
    output_yitang = data.frame(t(jth.row))
    colnames(output_yitang) <- column.names
    return(output_yitang)
    
  }
}

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
library(dplyr)

Pathway_out=c("/scratch/ys98038/UKB/plink2_format/COVID_19/Analyses/SNP/UKBiobank/Lipid_PRS/MEF/")

data(simulation.params)
data(pheno.params)
data(geno.params)
data(env.params)

str(simulation.params)
str(pheno.params)
str(geno.params)
str(env.params)

simulation.params_test=simulation.params
simulation.params_test$seed.val[4]=1
simulation.params_test$numsubjects[4]=260417
simulation.params_test$p.val[4]=5e-8

pheno.params_test=pheno.params
pheno.params_test$pheno.mean[4]=0
pheno.params_test$pheno.sd[4]=1
pheno.params_test$pheno.reliability[4]=1

geno.params_test=geno.params
geno.params_test$geno.efkt[4]=0

env.params_test=env.params

env.params_test$env.OR[4]=1.105171
env.params_test$env.efkt[4]=0.1
env.params_test$env.reliability[4]=1

####### @@@@@
env.params_test$env.model[4]=0

PUFA_name="Fish_Oil_Supplementation"
####### @@@@@
env.params_test$env.prevalence[4]=0.31

PUFA_name="Vegetarianism"
####### @@@@@
env.params_test$env.prevalence[4]=0.039

####### @@@@@
env.params_test$env.model[4]=1

PUFA_name="dietary_omega_3"
####### @@@@@
env.params_test$env.mean[4]=2.0
####### @@@@@
env.params_test$env.sd[4]=1.2

PUFA_name="dietary_omega_6"
####### @@@@@
env.params_test$env.mean[4]=11.0
####### @@@@@
env.params_test$env.sd[4]=6.1

PUFA_name="total_dietary_PUFAs"
####### @@@@@
env.params_test$env.mean[4]=13.0
####### @@@@@
env.params_test$env.sd[4]=6.2

######## $ env.sensitivity: num  0.99 0.99 0.99 0.99
######## $ env.specificity: num  0.99 0.99 0.99 0.99

res=data.frame("Interaction effect","MAF",	"N-genetic variants",	"Causal fraction",	"Empirical power")

write.table(res, file= paste(Pathway_out,"minimum_effect_size_",PUFA_name,".txt", sep=""), col.names = FALSE, append = TRUE,
            row.names = F, quote = FALSE, na = "",sep='\t')

if (env.params_test$env.model[4]==0) {
  for (var_num3 in c(0.0025,0.005,0.0075,0.01)) {
    geno.params_test$MAF[4] <- var_num3
    for (var_num4 in c(1,5,10,20)) {
      geno.params_test$geno.nvariants <- var_num4
      for (var_num5 in c(0.1,0.25,0.5,1)) {
        geno.params_test$geno.frac_causal <- var_num5
        
        # Remove  Columns in List
        simulation.params_test <- simulation.params_test[,!names(simulation.params_test) %in% c("numcases","numcontrols",
                                                                                                "interaction.OR","RR.5.95")]
        
        pheno.params_test <- pheno.params_test[,!names(pheno.params_test) %in% c("disease.prev","pheno.sensitivity",
                                                                                 "pheno.specificity")]
        
        geno.params_test <- geno.params_test[,!names(geno.params_test) %in% c("geno.OR","geno.sensitivity",
                                                                              "geno.specificity")]
        
        env.params_test <- env.params_test[,!names(env.params_test) %in% c("env.low.lim","env.up.lim",
                                                                           "env.sd","env.mean")]
        simulation.params_test$interaction.efkt[4]=0
        
        res=yitang_run.espresso.GxE.RV(simulation.params_test, pheno.params_test, geno.params_test,
                                       env.params_test, scenarios2run=c(4))
        print(str(res))
        
        while (res$empirical.power < 0.8) {
          simulation.params_test$interaction.efkt[4] = simulation.params_test$interaction.efkt[4] + 0.1
          
          res=yitang_run.espresso.GxE.RV(simulation.params_test, pheno.params_test, geno.params_test,
                                         env.params_test, scenarios2run=c(4))
          print(str(res))
        }
        
        if (res$interaction.efkt <= 0.1) {
          simulation.params_test$interaction.efkt[4]=0
          
          res=yitang_run.espresso.GxE.RV(simulation.params_test, pheno.params_test, geno.params_test,
                                         env.params_test, scenarios2run=c(4))
          print(str(res))
          
          while (res$empirical.power < 0.8) {
            simulation.params_test$interaction.efkt[4] = simulation.params_test$interaction.efkt[4] + 0.01
            
            res=yitang_run.espresso.GxE.RV(simulation.params_test, pheno.params_test, geno.params_test,
                                           env.params_test, scenarios2run=c(4))
            print(str(res))
          } 
          
          if (res$interaction.efkt <= 0.01) {
            simulation.params_test$interaction.efkt[4]=0
            
            res=yitang_run.espresso.GxE.RV(simulation.params_test, pheno.params_test, geno.params_test,
                                           env.params_test, scenarios2run=c(4))
            print(str(res))
            
            while (res$empirical.power < 0.8) {
              simulation.params_test$interaction.efkt[4] = simulation.params_test$interaction.efkt[4] +  0.001
              res=yitang_run.espresso.GxE.RV(simulation.params_test, pheno.params_test, geno.params_test,
                                             env.params_test, scenarios2run=c(4))
              print(str(res))
            } 
            if (res$interaction.efkt <= 0.001) {
              simulation.params_test$interaction.efkt[4]=0
              
              res=yitang_run.espresso.GxE.RV(simulation.params_test, pheno.params_test, geno.params_test,
                                             env.params_test, scenarios2run=c(4))
              print(str(res))
              
              while (res$empirical.power < 0.8) {
                simulation.params_test$interaction.efkt[4] = simulation.params_test$interaction.efkt[4] +  0.0001
                res=yitang_run.espresso.GxE.RV(simulation.params_test, pheno.params_test, geno.params_test,
                                               env.params_test, scenarios2run=c(4))
                print(str(res))
              } 
              if (res$interaction.efkt <= 0.0001) {
                simulation.params_test$interaction.efkt[4]=0
                
                res=yitang_run.espresso.GxE.RV(simulation.params_test, pheno.params_test, geno.params_test,
                                               env.params_test, scenarios2run=c(4))
                print(str(res))
                
                while (res$empirical.power < 0.8) {
                  simulation.params_test$interaction.efkt[4] = simulation.params_test$interaction.efkt[4] +  0.00001
                  res=yitang_run.espresso.GxE.RV(simulation.params_test, pheno.params_test, geno.params_test,
                                                 env.params_test, scenarios2run=c(4))
                  print(str(res))
                } 
              }
            }
          }
        } else {
          simulation.params_test$interaction.efkt[4] = res$interaction.efkt -0.1
          
          res=yitang_run.espresso.GxE.RV(simulation.params_test, pheno.params_test, geno.params_test,
                                         env.params_test, scenarios2run=c(4))
          print(str(res))
          
          while (res$empirical.power < 0.8) {
            simulation.params_test$interaction.efkt[4] = simulation.params_test$interaction.efkt[4] + 0.01
            res=yitang_run.espresso.GxE.RV(simulation.params_test, pheno.params_test, geno.params_test,
                                           env.params_test, scenarios2run=c(4))
            print(str(res))
          }
        }
        
        res=res %>% select(interaction.efkt,MAF,geno.nvariants,geno.frac_causal,empirical.power)
        
        print(str(res))

        write.table(res, file= paste(Pathway_out,"minimum_effect_size_",PUFA_name,".txt", sep=""), col.names = FALSE, append = TRUE,
                    row.names = F, quote = FALSE, na = "",sep='\t')
        
      }
    }
  }
}

if (env.params_test$env.model[4]==1) {
  for (var_num3 in c(0.0025,0.005,0.0075,0.01)) {
    geno.params_test$MAF[4] <- var_num3
    for (var_num4 in c(1,5,10,20)) {
      geno.params_test$geno.nvariants <- var_num4
      for (var_num5 in c(0.1,0.25,0.5,1)) {
        geno.params_test$geno.frac_causal <- var_num5
        
        simulation.params_test$interaction.efkt[4]=0
        
        # Remove  Columns in List
        simulation.params_test <- simulation.params_test[,!names(simulation.params_test) %in% c("numcases","numcontrols",
                                                                                                "interaction.OR","RR.5.95")]
        
        pheno.params_test <- pheno.params_test[,!names(pheno.params_test) %in% c("disease.prev","pheno.sensitivity",
                                                                                 "pheno.specificity")]
        
        geno.params_test <- geno.params_test[,!names(geno.params_test) %in% c("geno.OR","geno.sensitivity",
                                                                              "geno.specificity")]
        
        env.params_test <- env.params_test[,!names(env.params_test) %in% c("env.low.lim","env.up.lim","env.prevalence",
                                                                           "env.OR",
                                                                           "env.sensitivity","env.specificity")]
        
        res=yitang_run.espresso.GxE.RV(simulation.params_test, pheno.params_test, geno.params_test,
                                       env.params_test, scenarios2run=c(4))
        print(str(res))

        while (res$empirical.power < 0.8) {
          simulation.params_test$interaction.efkt[4] = simulation.params_test$interaction.efkt[4] + 0.1
          
          res=yitang_run.espresso.GxE.RV(simulation.params_test, pheno.params_test, geno.params_test,
                                         env.params_test, scenarios2run=c(4))
          print(str(res))
        }
        
        if (res$interaction.efkt <= 0.1) {
          simulation.params_test$interaction.efkt[4]=0
          
          res=yitang_run.espresso.GxE.RV(simulation.params_test, pheno.params_test, geno.params_test,
                                         env.params_test, scenarios2run=c(4))
          print(str(res))
          
          while (res$empirical.power < 0.8) {
            simulation.params_test$interaction.efkt[4] = simulation.params_test$interaction.efkt[4] + 0.01
            
            res=yitang_run.espresso.GxE.RV(simulation.params_test, pheno.params_test, geno.params_test,
                                           env.params_test, scenarios2run=c(4))
            print(str(res))
          } 
          
          if (res$interaction.efkt <= 0.01) {
            simulation.params_test$interaction.efkt[4]=0
            
            res=yitang_run.espresso.GxE.RV(simulation.params_test, pheno.params_test, geno.params_test,
                                           env.params_test, scenarios2run=c(4))
            print(str(res))
            
            while (res$empirical.power < 0.8) {
              simulation.params_test$interaction.efkt[4] = simulation.params_test$interaction.efkt[4] +  0.001
              res=yitang_run.espresso.GxE.RV(simulation.params_test, pheno.params_test, geno.params_test,
                                             env.params_test, scenarios2run=c(4))
              print(str(res))
            } 
            
            if (res$interaction.efkt <= 0.001) {
              simulation.params_test$interaction.efkt[4]=0
              
              res=yitang_run.espresso.GxE.RV(simulation.params_test, pheno.params_test, geno.params_test,
                                             env.params_test, scenarios2run=c(4))
              print(str(res))
              
              while (res$empirical.power < 0.8) {
                simulation.params_test$interaction.efkt[4] = simulation.params_test$interaction.efkt[4] +  0.0001
                res=yitang_run.espresso.GxE.RV(simulation.params_test, pheno.params_test, geno.params_test,
                                               env.params_test, scenarios2run=c(4))
                print(str(res))
              } 

              if (res$interaction.efkt <= 0.001) {
                simulation.params_test$interaction.efkt[4]=0
                
                res=yitang_run.espresso.GxE.RV(simulation.params_test, pheno.params_test, geno.params_test,
                                               env.params_test, scenarios2run=c(4))
                print(str(res))
                
                while (res$empirical.power < 0.8) {
                  simulation.params_test$interaction.efkt[4] = simulation.params_test$interaction.efkt[4] +  0.00001
                  res=yitang_run.espresso.GxE.RV(simulation.params_test, pheno.params_test, geno.params_test,
                                                 env.params_test, scenarios2run=c(4))
                  print(str(res))
                } 
              }
            }
          }
        } else {
          simulation.params_test$interaction.efkt[4] = res$interaction.efkt -0.1
          
          res=yitang_run.espresso.GxE.RV(simulation.params_test, pheno.params_test, geno.params_test,
                                         env.params_test, scenarios2run=c(4))
          print(str(res))
          
          while (res$empirical.power < 0.8) {
            simulation.params_test$interaction.efkt[4] = simulation.params_test$interaction.efkt[4] + 0.01
            res=yitang_run.espresso.GxE.RV(simulation.params_test, pheno.params_test, geno.params_test,
                                           env.params_test, scenarios2run=c(4))
            print(str(res))
          }
        }
        
        res=res %>% select(interaction.efkt,MAF,geno.nvariants,geno.frac_causal,empirical.power)
        
        print(str(res))

        write.table(res, file= paste(Pathway_out,"minimum_effect_size_",PUFA_name,".txt", sep=""), col.names = FALSE, append = TRUE,
                    row.names = F, quote = FALSE, na = "",sep='\t')
      }
    }
  }
}
