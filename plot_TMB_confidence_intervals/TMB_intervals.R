
setwd( "C:/Users/James.Thorson/Desktop/Git/mixed-effects/plot_TMB_confidence_intervals" )
library(TMB)

###########
# Delta-model for canary rockfish
###########

# Generate data
n_data = 30
X = runif( n_data )
Y = 1 + 0.2*X + rnorm( n_data )

# Step 1 -- make and compile template file
compile( "TMB_intervals.cpp" )

# Step 2 -- build inputs and object
dyn.load( dynlib("TMB_intervals") )
Data = list( "y_i"=Y, "X_ij"=cbind(1,X) )
Params = list( "b_j"=rep(0,ncol(Data$X_ij)), "log_sd"=0 )
Obj = MakeADFun( data=Data, parameters=Params, DLL="TMB_intervals")

# Step 3 -- test and optimize
Opt = nlminb( start=Obj$par, objective=Obj$fn, gr=Obj$gr )
SD = sdreport( Obj )

# Step 4 -- Simulate from predictive distribution
match_index = grep( "b_j", names(Opt$par) )
bhat_rj = mvtnorm::rmvnorm( n=1e4, mean=Opt$par[match_index], sigma=SD$cov.fixed[match_index,match_index] )

# predict response for new values
Xpred_z = seq( from=-10, to=10, length=1000 )
Ybounds_zj = matrix( NA, ncol=2, nrow=length(Xpred_z) )
for( z in 1:nrow(Ybounds_zj) ){
  ysim_r = bhat_rj[,1] + bhat_rj[,2]*Xpred_z[z]
  Ybounds_zj[z,] = quantile( ysim_r, prob=c(0.1,0.9) )
}

# plot results
plot( x=X, y=Y )
abline( a=Opt$par[match_index][1], b=Opt$par[match_index][2] )
polygon( x=c(Xpred_z,rev(Xpred_z)), y=c(Ybounds_zj[,1],rev(Ybounds_zj[,2])), col=rgb(1,0,0,0.2) )

