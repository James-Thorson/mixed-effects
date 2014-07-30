
File = "C:/Users/James.Thorson/Desktop/UW Hideaway/Collaborations/2013 -- Random effect review/Code/Individual variability/"
setwd(File)

library(abind)
library(R2jags); runif(1)

# Data from : http://datadryad.org/resource/doi:10.5061/dryad.cj36j/1
Data = read.csv( "Steelhead Size Data.csv", header=TRUE )

# Subset
Data = Data[which(Data[,'Treatment']==4 & Data[,'Strain']=="Coleman"),]
  Data = Data[which(Data[,'Ration']=="Low"),]
  Data[,'Date'] = as.Date(as.character(Data[,'Date']), "%m/%d/%y")

# Change format
Array = NULL
for(i in 1:length(unique(Data[,'Unique.ID']))){
  Temp_1 = Data[which(Data[,'Unique.ID']==unique(Data[,'Unique.ID'])[i]),]
  Temp_1 = Temp_1[which(!is.na(as.numeric(as.character(Temp_1[,'Unique.ID'])))),]
  Temp_1 = cbind( 'ID'=as.numeric(as.character(Temp_1[,'Unique.ID'])), 'Days'=c(NA,diff(Temp_1[,'Date'])), 'Date_after_start'=Temp_1[,'Date']-min(Temp_1[,'Date']), 'Growth'=c(NA,diff(Temp_1[,'Length.mm'])), 'L_t'=c(Temp_1[,'Length.mm']), 'L_prev'=c(NA,Temp_1[-nrow(Temp_1),'Length.mm']) )
  # Array representation
  if(nrow(Temp_1)>0){
    if(is.null(Array)) Array = Temp_1
    if(!is.null(Array)) Array = abind( Array, Temp_1, along=3)
  }
}

# JAGS example -- Mixed-effects model
Model_1 = "
model {
  ### Priors
  # Hyperparameters for k
  k_mean ~ dunif(0,1)     
  ln_k_mean <- log(k_mean) 
  CV_k <- pow( exp( 2*sigma_k )-1, 0.5 )
  sigma_k ~ dunif(0,1)
  tau2_k <- pow(sigma_k,-2)
  # Derived quantities for q
  q_mean <- gamma*pow(k_mean,Psi)
  # Growth parameters
  gamma ~ dunif(0,20)
  # Measurement errors
  sigma ~ dunif(0,100)
  tau2 <- pow(sigma,-2)
  # Hyperparameters for length at first sighting
  L0_mean ~ dunif(0,100)
  sigma_L0 ~ dunif(0,100)
  tau2_L0 <- pow(sigma_L0,-2)
  # States
  for(IdentI in 1:Nident){
    # k random effect for each individual
    k_ident[IdentI] ~ dlnorm(ln_k_mean,tau2_k) T(0,1)
    q_ident[IdentI] <- gamma*pow(k_ident[IdentI],Psi)
    # L0 random effect for each individual
    L_hat[1,IdentI] ~ dnorm( L0_mean, tau2_L0 )
    L_t[1,IdentI] ~ dnorm(L_hat[1,IdentI],tau2)
    # Predicting length for each observation, and computing probability of data 
    for(PeriodI in 2:Nperiod){
      L_hat[PeriodI,IdentI] <- L_hat[PeriodI-1,IdentI]*exp(-k_ident[IdentI]*Days[PeriodI,IdentI]) + gamma*pow(k_ident[IdentI],Psi-1)*(1-exp(-k_ident[IdentI]*Days[PeriodI,IdentI]))
      L_t[PeriodI,IdentI] ~ dnorm(L_hat[PeriodI,IdentI],tau2)
    }
  }
  # Predictive distribution (for plotting growth of unobserved individuals)
  k_pred ~ dlnorm(ln_k_mean,tau2_k) T(0,1) 
  L_pred[1] ~ dnorm( L0_mean, tau2_L0 )
  for(DayI in 2:MaxDaysPred){
    L_pred[DayI] <- L_pred[DayI-1]*exp(-k_pred) + pow(k_pred,Psi-1)*(1-exp(-k_pred))*gamma
  }
  q_pred <- gamma*pow(k_pred,Psi)
}
"
cat(Model_1, file="model_1.bug")
Nsim = Nburnin = 1e6
Nthin = 1e3
Data = list(Nident=dim(Array)[3], Nperiod=dim(Array)[1], MaxDays=max(Array[,'Days',],na.rm=TRUE), MaxDaysPred=max(Array[,'Date_after_start',],na.rm=TRUE), Days=Array[,'Days',], L_t=Array[,'L_t',], Psi=0.2)
# R2jags
Jags_1 <- jags(model.file="model_1.bug", working.directory=NULL, data=Data, parameters.to.save=c("gamma","k_mean","L0_mean","k_ident","sigma","sigma_k","sigma_L0","q_mean","q_ident","Linf_mean","Linf_ident","k_mean","k_ident","CV_k","L_hat","L_pred"), n.chains=3, n.thin=Nthin, n.iter=Nsim+Nburnin, n.burnin=Nburnin)
Jags_1$BUGSoutput$summary

# JAGS example -- No random effects
Model_0 = "
model {
  # Priors
  k_mean ~ dunif(0,1)
  gamma ~ dunif(0,20)
  sigma ~ dunif(0,100)
  tau2 <- pow(sigma,-2)
  L0_mean ~ dunif(0,100)
  # States
  for(IdentI in 1:Nident){
    L_hat[1,IdentI] <- L0_mean
    L_t[1,IdentI] ~ dnorm( L_hat[1,IdentI], tau2 )
    for(PeriodI in 2:Nperiod){
      L_hat[PeriodI,IdentI] <- L_hat[PeriodI-1,IdentI]*exp(-k_mean*Days[PeriodI,IdentI]) + gamma*pow(k_mean,Psi-1)*(1-exp(-k_mean*Days[PeriodI,IdentI]))
      L_t[PeriodI,IdentI] ~ dnorm(L_hat[PeriodI,IdentI],tau2)
    }
  }
  # Predictive
  L_pred[1] <- L0_mean
  for(DayI in 2:MaxDaysPred){
    L_pred[DayI] <- L_pred[DayI-1]*exp(-k_mean) + pow(k_mean,Psi-1)*(1-exp(-k_mean))*gamma
  }
}
"
cat(Model_0, file="model_0.bug")
Nsim = Nburnin = 1e5
Nthin = 1e2
Data = list(Nident=dim(Array)[3], Nperiod=dim(Array)[1], MaxDays=max(Array[,'Days',],na.rm=TRUE), MaxDaysPred=max(Array[,'Date_after_start',],na.rm=TRUE), Days=Array[,'Days',], L_t=Array[,'L_t',], Psi=0.2)
# R2jags
Jags_0 <- jags(model.file="model_0.bug", working.directory=NULL, data=Data, parameters.to.save=c("gamma","k_mean","L0_mean","k_ident","sigma","sigma_k","sigma_L0","q_mean","q_ident","Linf_mean","Linf_ident","k_mean","k_ident","CV_k","L_hat","L_pred"), n.chains=3, n.thin=Nthin, n.iter=Nsim+Nburnin, n.burnin=Nburnin)
Jags_0$BUGSoutput$summary

# Plot data vs. estimates
L_hat_1 = apply( Jags_1$BUGSoutput$sims.list$L_hat, MARGIN=2:3, FUN=median )
L_pred_1 = apply( Jags_1$BUGSoutput$sims.list$L_pred, MARGIN=2, FUN=quantile, prob=c(0.025,0.5,0.975) )
L_pred_0 = apply( Jags_0$BUGSoutput$sims.list$L_pred, MARGIN=2, FUN=quantile, prob=c(0.025,0.5,0.975) )
png(file=paste(File,"Growth.png",sep=""), width=4, height=4, res=400, units="in")
  par(mar=c(3,3,0.7,0)+0.1, mgp=c(2,0.5,0), tck=-0.02, xaxs="i", yaxs="i")
  # data
    matplot( y=Array[,'L_t',], x=Array[,'Date_after_start',], type="l", pch=20, lwd=0.5, col="black", lty="solid", xlab="", ylab="Length (mm)", ylim=c(0,max(Array[,'L_t',],na.rm=TRUE)) ) 
    mtext(side=1, text="Days after start of experiment", line=2)
  # predictions of data
    #matplot( y=L_hat_1, x=Array[,'Date_after_start',], type="l", pch=20, col="red", lty="solid", add=TRUE ) 
  # Predictive distribution (variation)
    lines( y=L_pred_1[2,], x=1:Data[['MaxDaysPred']], col="blue", lwd=3)
    polygon( y=c(L_pred_1[1,],rev(L_pred_1[3,])), x=c(1:Data[['MaxDaysPred']],Data[['MaxDaysPred']]:1), col=rgb(0,0,1,0.2), border=NA)
  # predictive distribution (conventional)
    lines( y=L_pred_0[2,], x=1:Data[['MaxDaysPred']], col="red", lwd=3)
    polygon( y=c(L_pred_0[1,],rev(L_pred_0[3,])), x=c(1:Data[['MaxDaysPred']],Data[['MaxDaysPred']]:1), col=rgb(1,0,0,0.2), border=NA)
dev.off()

# Save record
save(Jags_0, file=paste(File,"Jags_0.RData",sep=""))
save(Jags_1, file=paste(File,"Jags_1.RData",sep=""))
capture.output(Jags_0, file=paste(File,"Jags_0_summary.txt",sep=""))
capture.output(Jags_1, file=paste(File,"Jags_1_summary.txt",sep=""))
# summary
summary(Jags$BUGSoutput$sims.list$k_mean)
# Individual-specific k
mean(apply(Jags$BUGSoutput$sims.list$k_ident, MARGIN=2, FUN=mean))
sd(apply(Jags$BUGSoutput$sims.list$k_ident, MARGIN=2, FUN=mean))
# Individual-specific Linf
mean(apply(Jags$BUGSoutput$sims.list$Linf_ident, MARGIN=2, FUN=mean))
sd(apply(Jags$BUGSoutput$sims.list$Linf_ident, MARGIN=2, FUN=mean))
# Individual-specific q
mean(apply(Jags$BUGSoutput$sims.list$q_ident, MARGIN=2, FUN=mean))
sd(apply(Jags$BUGSoutput$sims.list$q_ident, MARGIN=2, FUN=mean))
