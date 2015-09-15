
######################
# Simulate data
######################

Factor = rep( 1:5, each=5)
Z = rnorm( 5, mean=0, sd=1)

X0 = 0

Y = Z[Factor] + X0 + rnorm( 25, mean=0, sd=1)

######################
# Run in R
######################

library(lme4)

Lme = lmer( Y ~ 1|factor(Factor))

######################
# Run in TMB
######################

library(TMB)

compile( "linear_mixed_model.cpp" )

Data = list( "n_data"=25, "n_factors"=5, "Factor"=Factor-1, "Y"=Y)
Parameters = list( "X0"=-10, "log_SD0"=2, "log_SDZ"=2, "Z"=rep(0,Data$n_factor) )
Random = c("Z")

dyn.load( dynlib("linear_mixed_model") )
Obj = MakeADFun(data=Data, parameters=Parameters, random=Random)

# Prove that function and gradient calls work
Obj$fn( Obj$par )
Obj$gr( Obj$par )

# Optimize
Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, control=list("trace"=1) )
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
