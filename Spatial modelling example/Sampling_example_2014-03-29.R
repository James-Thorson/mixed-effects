
####################
# NOTES
#
# BY DEFAULT: Alpha = 3/2
# Nu = 3/2 - D/2, where D = number of dimensions (i.e. Nu = 1/2)
# Kappa = exp(theta2)
# Range = sqrt(8*Nu)/Kappa
#
####################

#source("http://www.math.ntnu.no/inla/givemeINLA.R")
#inla.upgrade()
setwd("C:/Users/James.Thorson/Desktop/UW Hideaway/Collaborations/2013 -- Random effect review/Code/Spatial model")

require(INLA) # Download from: http://www.r-inla.org/download
require(lattice)
require(gridExtra)

n = 1000
s1 = runif(n, min=0, max=1)
s2 = runif(n, min=0, max=1)
a = 1
# spatial process                                    
ln_y_hat = 2.0*exp(-((s1-s2)/0.2)^2) + a*s1                      
## make simulated data with no spatial component     
y = rpois(n, lambda=exp(ln_y_hat))                 
                                                     
# Fit model -- INLA
coords <- cbind( s1, s2 )
data = data.frame(y=y,coords=coords,x=coords[,1])
  data = na.omit(data)
  coords = as.matrix(data[,2:3])
pl.dom <- cbind( c(min(coords[,1]),max(coords[,1]),max(coords[,1]),min(coords[,1])), c(min(coords[,2]),min(coords[,2]),max(coords[,2]),max(coords[,2])))
mesh <- inla.mesh.2d(coords, pl.dom, max.edge=1, offset=-0.5)
# Configure model
spde <- inla.spde2.matern(mesh, alpha=2)
A <- inla.spde.make.A(mesh, loc=coords)                                                                 
stk <- inla.stack(data=list(resp=data$y), A=list(A,1,1), effects=list(i=1:spde$n.spde, intercept=rep(1,nrow(data)), x=data$x), tag='est')
formula = resp ~ 0 + intercept + x + f(i, model=spde)
# Run model
result = inla(formula, family="poisson", data=inla.stack.data(stk), verbose=TRUE, control.predictor=list(A=inla.stack.A(stk), compute=TRUE), control.family=list(hyper = list(theta = list(initial = log(1/0.01^2), fixed = FALSE))), keep=FALSE)
summary(result)

# Calculate true expected surface - INLA
s1_proj = outer( seq(0,1,length=100), rep(1,100) )
s2_proj = outer( rep(1,100), seq(0,1,length=100) )
ln_y_exp = 2.0*exp(-((s1_proj-s2_proj)/0.2)^2) + a*s1_proj
# Caclualte estimated expected surface
mesh_proj = inla.mesh.projector(mesh, xlim=c(0,1), ylim=c(0,1), dims=c(100,100))
ln_y_hat = inla.mesh.project(mesh_proj, result$summary.ran$i$mean) + result$summary.fix['intercept','mean'] + result$summary.fix['x','mean']*s1_proj

# Fit model -- GLM
Glm = glm( y ~ s1 + s2 + I(s1^2) + I(s2^2) + I(s1*s2), family="poisson")
Glm = step( Glm )
ln_y_hat_glm = array(predict( Glm, type="link", newdata=data.frame("s1"=as.vector(s1_proj), "s2"=as.vector(s2_proj))), dim=dim(s1_proj))
plot( x=exp(ln_y_hat_glm), y=exp(ln_y_exp))

# Figure showing data and mesh
png(file="INLA_example--Data.png", width=6, height=3, res=200, units="in")
  # Plot data and mesh
  par(mfrow=c(1,2), mar=c(2.5,2.5,1,0), mgp=c(1.5,0.25,0), tck=-0.02)
  # data
  plot(x=s1, y=s2, cex=sqrt(y+0.1)/2, type="p", main="Data")
  # mesh
  plot(mesh, asp=1, main="")
  lines(rbind(pl.dom,pl.dom[1,]), col="black", lwd=2)
  title("Regions")
  #points(coords, col="black", cex=0.1)
dev.off()
# Figure showing expected surface, and predictions from INLA and GLM
Col = colorRampPalette(colors=c("blue","white","red"))
At = seq( min(ln_y_exp,ln_y_hat,ln_y_hat_glm),max(ln_y_exp,ln_y_hat,ln_y_hat_glm), length=100)
WhichNum = function(Vec){
  Num = rep(NA, length(Vec))
  for(i in 1:length(Num)) Num[i] = which.max( At>Vec[i] )
  Num = ifelse( Num==length(At), length(At)-1, Num)
  return(Num) 
}
png(file="INLA_example--True_vs_Pred.png", width=6, height=3, res=200, units="in")
  grid.arrange(
    levelplot(ln_y_exp, xlab='', ylab='', main='ln(True value)', col.regions=rev(topo.colors(99)), at=At, scales=list(draw=FALSE)),
    levelplot(ln_y_hat, xlab='', ylab='', main='ln(Exp. value)', col.regions=rev(topo.colors(99)), at=At, scales=list(draw=FALSE)), 
    nrow=1
  )
dev.off()
# Figure showing true, pred, and data
Col = colorRampPalette(colors=c("red","purple","blue","green","yellow"))
png(file="INLA_example--True Pred-INLA Pred-GLM.png", width=3, height=7.5, res=400, units="in")
  grid.arrange(    # trellis.par.get()
    levelplot(ln_y_exp, xlab='', ylab='', main=expression(bold("E[ ln(y) ]")), col.regions=(Col(99)), at=At, scales=list(draw=FALSE), par.settings=list(layout.heights=list(top.padding=0,bottom.padding=-1)) ), # 
    #xyplot( s2 ~ s1, cex=1, pch=20, col=rev(topo.colors(99))[WhichNum(log(y))], mar=c(5,5,5,5), xlab="", ylab="", main="y", scales=list(tck=0.05), par.settings=list(layout.heights=list(top.padding=0,bottom.padding=-1),layout.widths=list(right.padding=0,left.padding=-2), par.main.text=list(alpha=1, lineheight=0, cex=1.3)) ),
    levelplot(ln_y_hat, xlab='', ylab='', main=expression(paste(bold("INLA:ln("),bold(hat(y)),")",sep=" ")), col.regions=(Col(99)), at=At, scales=list(draw=FALSE), par.settings=list(layout.heights=list(top.padding=0,bottom.padding=-1)) ), 
    levelplot(ln_y_hat_glm, xlab='', ylab='', main=expression(paste(bold("GLM:ln("),bold(hat(y)),")",sep=" ")), col.regions=(Col(99)), at=At, scales=list(draw=FALSE), par.settings=list(layout.heights=list(top.padding=0,bottom.padding=-1)) ), 
    #lplot.xy(xy=1),
    nrow=3
  )
dev.off()

# Explore GRF parameters
Interp = inla.spde2.result(result, "i", spde)
par(mfrow=c(1,2))
plot(Interp[["marginals.range.nominal"]][[1]])
plot(Interp[["marginals.variance.nominal"]][[1]])
# Explore Thetas
theta1 = result$summary.hy[1,'mean']
theta2 = result$summary.hy[2,'mean']
sqrt(exp(theta1))            # Positively correlated with range
1/sqrt(exp(theta2))          # Negatively correlated with maximum correction (i.e., at zero distance)
# Summarize results           # family = "poisson"
post.se <- inla.tmarginal(function(x) sqrt(1/x),res$marginals.hy[[1]])   # kappa -> 1/SD for covariance;
res.field <- inla.spde2.result(res, 'i', spde, do.transf=TRUE)           # sqrt(exp(theta1)) -> range; 1/sqrt(exp(theta2)) -> variance    #v1: (-0.1, -1.2)
inla.emarginal(function(x) x, res.field$marginals.range.nominal[[1]])
# Make predictions
Pred = exp(result$summary.linear.pred[1:nrow(coords),'mean'])
Obs = data$y
matplot( x=Obs, y=Pred, type="p", pch=19 )
# Spatial prediction
result$summary.ran
result$summary.fix
result$summary.fitt[1:(nrow*ncol),'mean']
require(gridExtra)
PredMat = y
  PredMat[which(!is.na(y))] = result$summary.fitt[1:(nrow(data)),'mean']
