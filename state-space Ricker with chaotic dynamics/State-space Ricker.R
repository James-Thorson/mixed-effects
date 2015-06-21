
library( TMB )

# Parameters
r = 4 # Determinstic chaos r>2.692
k = 1
sigmap = 0.1
sigmam = 0.1
nobs = 1000

# Simulate
n = rep(NA,nobs)
n[1] = exp(rnorm(1,log(k),0.1))
for(t in 2:length(n)) n[t] = n[t-1]*exp(r*(1-n[t-1]/k))*exp(rnorm(1,0,sigmap))
nobs = n * exp(rnorm(1,0,sigmam))

# Visualize
par(mfrow=c(1,2))
plot(n, type="l")
plot( x=n[-length(n)], y=n[-1]-n[-length(n)] )

#################
# Innovations parameterization
#################
Version_innovations = "innovations_v1"
compile( paste0(Version_innovations,".cpp") )

# Fit the model
Data = list("y_t"=n)
Params = list("log_r"=log(1), "log_k"=log(1), "log_sigmap"=log(1), "log_sigmam"=log(1), "e_t"=rnorm(length(n)) )
dyn.load( dynlib(Version_innovations) )
Obj_0 = MakeADFun( data=Data, parameters=Params, random="e_t")
Opt_0 = nlminb(start=Obj_0$par, objective=Obj_0$fn, gradient=Obj_0$gr, control=list("trace"=1))
Opt_0[["final_gradient"]] = Obj_0$gr( Opt_0$par )

#################
# Innovations parameterization
#################
Version_statespace = "statespace_v1"
compile( paste0(Version_statespace,".cpp") )

# Fit the model
Data = list("y_t"=n)
Params = list("log_r"=log(1), "log_k"=log(1), "log_sigmap"=log(1), "log_sigmam"=log(1), "log_x_t"=rnorm(length(n)) )
dyn.load( dynlib(Version_statespace) )
Obj_1 = MakeADFun( data=Data, parameters=Params, random="log_x_t")
Opt_1 = nlminb(start=Obj_1$par, objective=Obj_1$fn, gradient=Obj_1$gr, control=list("trace"=1))
Opt_1[["final_gradient"]] = Obj_1$gr( Opt_1$par )
