
setwd( "C:/Users/James.Thorson/Desktop/Project_git/mixed-effects/linear_mixed_model" )

######################
# Simulate data
######################

Factor = rep( 1:1000, each=1000)
Z = rnorm( length(unique(Factor)), mean=0, sd=1)

X0 = 0

Y = Z[Factor] + X0 + rnorm( length(Factor), mean=0, sd=1)

######################
# Run in R
######################

library(lme4)

Lme = lmer( Y ~ 1|factor(Factor))

######################
# Run in TMB
######################

library(TMB)

Version = c( "linear_mixed_model", "linear_mixed_model_parallel")[2]
compile( paste0(Version,".cpp") )

Data = list( "n_data"=length(Y), "n_factors"=length(unique(Factor)), "Factor"=Factor-1, "Y"=Y)
Parameters = list( "X0"=-10, "log_SD0"=2, "log_SDZ"=2, "Z"=rep(0,Data$n_factor) )
Random = c("Z")

dyn.load( dynlib("linear_mixed_model") )
Obj = MakeADFun(data=Data, parameters=Parameters, random=Random)

# Prove that function and gradient calls work
Obj$fn( Obj$par )
Obj$gr( Obj$par )

# Optimize
start_time = Sys.time()
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, control=list("trace"=1) )
  Opt[["final_gradient"]] = Obj$gr( Opt$par )
  Opt[["total_time"]] = Sys.time() - start_time
  Report = Obj$report()
  SD = sdreport( Obj, bias.correct=TRUE )

# Bias correct experiment
SD$unbiased$value['SampleVarZ']
Report$SampleVarZ  
SD$unbiased$value['SampleSDZ']
Report$SampleSDZ
SD$unbiased$value['SDZ']
Report$SDZ  

######################
# Compare estimates
######################

# Global mean
c( fixef(Lme), Report$X0 )

# Global mean
cbind( "True"=Z, ranef(Lme)[['factor(Factor)']], Report$Z )

# Variances
summary(Lme)
unlist( Report[c("SDZ","SD0")] )

#########################
# Experiment with parallelization
#########################

# Parallelized experimentations
ben <- benchmark(Obj, n=100, cores=1:4)
plot(ben)

# Parallel
ben <- benchmark(Obj, n=1, cores=1:2, expr=expression(do.call("nlminb",obj)))
plot(ben)

