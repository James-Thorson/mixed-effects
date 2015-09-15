
#########################
# Spatial Gompertz model
# SEE: James T. Thorson, Hans Skaug, Kasper Kristensen, Andrew O. Shelton, Eric J. Ward, John Harms, Jim Benante. In press. The importance of spatial models for estimating the strength of density dependence. Ecology.
#########################

setwd( "C:/Users/James.Thorson/Desktop/Project_git/mixed-effects/spatial_Gompertz_model" )

# load libraries
library(INLA)
library(TMB)
library(RandomFields)

source( "Sim_Gompertz_Fn.R" )

# Settings
Species = c('Saved_example', "Simulated_example")[2]
  MeshType = c("Samples", "Refined")[2] # Samples: faster; Refined: slower and more accurate approximation to GRMF

# Read data
if( Species=='Saved_example' ){
  Data = read.csv( file="spatial_gompertz_simulated_data.csv", header=TRUE)
}
if( Species=='Simulated_example' ){
  # n_years=10; n_stations=100; SpatialScale=0.1; SD_O=0.5; SD_E=0.2; SD_extra=0; rho=0.8; logMeanDens=1; phi=NULL; Loc=NULL
  set.seed( 1 )
  Sim_List = Sim_Gompertz_Fn( n_years=10, n_stations=100, SpatialScale=0.1, SD_O=0.4, SD_E=0.2, SD_extra=0, rho=0.5, logMeanDens=1, phi=0.0, Loc=NULL )
  Data = Sim_List[["DF"]]
}

# Length of data
n_years = length(unique(Data$Year))
n_stations = length(unique(Data$Site))
n_data = n_stations*n_years
x_stations = Data[match(unique(Data$Site),Data$Site),'Lon..DDD.DDDDD.']
y_stations = Data[match(unique(Data$Site),Data$Site),'Lat..DD.DDDDD.']

# display stations
plot( x=x_stations, y=y_stations)

# reformat data
Ymat = matrix(NA, nrow=n_stations, ncol=n_years)
for(YearI in 1:n_years){
for(SiteI in 1:n_stations){
  Which = which(Data$Year==unique(Data$Year)[YearI] & Data$Site==unique(Data$Site)[SiteI])
  if(length(Which)>=1) Ymat[SiteI,YearI] = Data[Which[1],Species]
}}
# vectorize data
Y = as.vector(Ymat)		  
Site = as.vector(row(Ymat))
Year = as.vector(col(Ymat))
NAind = as.integer(ifelse(is.na(Y),1,0))

# Show data
tapply( Y, FUN=mean, INDEX=Year)

# Build SPDE object using INLA (must pass mesh$idx$loc when supplying Boundary)
Cutoff = 1e-12 # 1e-1 
Boundary = NULL # inla.nonconvex.hull(points=cbind(x_stations, y_stations), convex=-0.1)	
if(MeshType=="Samples") mesh = inla.mesh.create( cbind(x_stations, y_stations), boundary=Boundary, plot.delay=NULL, cutoff=Cutoff, extend=list(n=8,offset=-0.15), refine=F )  # loc_samp
if(MeshType=="Refined") mesh = inla.mesh.create( cbind(x_stations, y_stations), boundary=Boundary, plot.delay=NULL, cutoff=Cutoff, extend=list(n=8,offset=-0.15), refine=list(min.angle=26) )  # loc_samp  ;  ,max.edge.data=0.08,max.edge.extra=0.2
spde = inla.spde2.matern(mesh,alpha=2)

# Settings
newtonOption(smartsearch=TRUE)

###################
# Parameter estimation
###################

Version = c("spatial_gompertz", "spatial_gompertz_state_as_random", "spatial_gompertz_aniso", "spatial_gompertz_parallel")[2]

# Run spatial model
if( Version %in% c("spatial_gompertz","spatial_gompertz_parallel") ){
  X = cbind( rep(1,n_data) )
  Data = list(n_data=n_stations*n_years, Y=Y, NAind=NAind, n_knots=mesh$n, n_stations=n_stations, meshidxloc=mesh$idx$loc-1, n_years=n_years, n_p=ncol(X), X=X, G0=spde$param.inla$M0, G1=spde$param.inla$M1, G2=spde$param.inla$M2)
  Parameters = list(alpha=c(0.0), phi=0.0, log_tau_E=0.0, log_tau_O=0.0, log_kappa=0.0,	rho=0.5, Epsilon_input=matrix(rnorm(mesh$n*n_years),nrow=mesh$n,ncol=n_years), Omega_input=rnorm(mesh$n))
  Random = c("Epsilon_input","Omega_input")
}
if( Version=="spatial_gompertz_state_as_random"){
  X = cbind( rep(1,mesh$n*n_years) )
  Data = list(n_data=n_stations*n_years, Y=Y, NAind=NAind, n_knots=mesh$n, n_stations=n_stations, meshidxloc=mesh$idx$loc-1, n_years=n_years, n_p=ncol(X), X=X, G0=spde$param.inla$M0, G1=spde$param.inla$M1, G2=spde$param.inla$M2)
  Parameters = list(alpha=c(0.0), phi=0.0, log_tau_U=1.0, log_tau_O=1.0, log_kappa=0.0,	rho=0.5, log_Dji=matrix(rnorm(mesh$n*n_years),nrow=mesh$n,ncol=n_years), Omega_input=rnorm(mesh$n))
  Random = c("log_Dji","Omega_input")
}

# Fix parameters (optional)
Map = list()
if( FALSE ){
  Map[["rho"]] = factor(NA)
  Parameters[["rho"]] = 0
}

# Make object
compile( paste0(Version,".cpp") )
dyn.load( dynlib(Version) )
start_time = Sys.time()
obj <- MakeADFun(data=Data, parameters=Parameters, random=Random, map=Map, hessian=FALSE)
  #obj$par <- opt0$par
obj$fn(obj$par)
obj$gr(obj$par)
#Report = obj$report()

# Run optimizer
obj$control <- list(trace=1,parscale=rep(1,13),REPORT=1,reltol=1e-12,maxit=100)
opt = nlminb(obj$par, objective=obj$fn, gradient=obj$gr, lower=c(rep(-20,2),rep(-10,3),-0.999), upper=c(rep(20,2),rep(10,3),0.999), control=list(eval.max=1e4, iter.max=1e4, trace=1))
opt[["final_gradient"]] = obj$gr( opt$par )
opt[["total_time"]] = Sys.time() - start_time
unlist( Report[c('Range','SigmaO','SigmaE','SigmaU')] )
Report = obj$report()

# Get standard errors
SD = try( sdreport(obj) )

# Parallelized experimentations
ben <- benchmark(obj,n=100,cores=1:3)
plot(ben)

# Parallel
ben <- benchmark(obj,n=1,cores=1:2,expr=expression(do.call("nlminb",obj)))
plot(ben)

