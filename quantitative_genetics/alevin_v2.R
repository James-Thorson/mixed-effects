##--------------------------------------------
## Genetic example
## Animal model code for Debes et al. 2013 alevin 
## Data from:
## http://datadryad.org/resource/doi:10.5061/dryad.9cs2v
##--------------------------------------------

## load libraries
library(reshape)
## pedigree vizualization
library(pedantics)
## this should also load MCMCglmm
## library(MCMCglmm)

## PEDIGREE DATA
pedigree.dat<-read.table("http://datadryad.org/bitstream/handle/10255/dryad.48220/Pedigree_ALEVIN.txt", header=TRUE, stringsAsFactors=FALSE)

## exploratory plots of the pedigree
## males: tan; females: green
drawPedigree(pedigree.dat, sexColours=c("#D8B365","#5AB4AC"),lwd=0.3)

## PHENOTYPIC DATA
offspring.dat<-read.table("http://datadryad.org/bitstream/handle/10255/dryad.48218/Offspring_Data.txt", header=TRUE, stringsAsFactors=FALSE)

## subset for alevins
alevin.dat<-subset(offspring.dat,  Stage=="02_Alevin")

##----------------
## MCMCglmm SETUP
##----------------

## remove animals without pedigree
idx<-alevin.dat$Ind_id%in%pedigree.dat$id
alevin.dat<-alevin.dat[idx,]

## "animal" identifier for MCMCglmm 
alevin.dat$animal <- as.factor(alevin.dat$Ind_id)

## factor variables
alevin.dat$Dam<-as.factor(alevin.dat$Dam)

## cross
alevin.dat$Reciprocal_Cross<-as.factor(alevin.dat$Reciprocal_Cross)

## tank
alevin.dat$Tank_ID<-as.factor(alevin.dat$Tank_ID)

## family
alevin.dat$ID<-as.factor(alevin.dat$ID)

##----------------------------------
## MCMCglmm RUN WITH RANDOM EFFECTS
##----------------------------------
## prior
## split the variance in 5 among the four components
bwt.var<-var(alevin.dat$body_weight, na.rm=TRUE)

prior1 <- list(G = list(G1 = list(V = bwt.var/5, nu = 1), G2 = list(V = bwt.var/5,nu = 1), G3 = list(V = bwt.var/5,nu=1), G4 = list(V = bwt.var/5, nu = 1)), R = list(V = bwt.var/5, nu = 1))

## number of mcmc iterations
nitt<-1e6

## thinning interval to obtain 1e3 posterior samples from first half burnin chain
thin<-nitt/(1e3*2)

## run the animal model
## takes a while
## compare correlation

alevin.bwt.1 <- MCMCglmm(body_weight ~ Reciprocal_Cross + degree_days_centered, random = ~animal + Dam + Tank_ID + ID, pedigree = pedigree.dat, data = alevin.dat, nitt =nitt, thin = thin, burnin=nitt/2, prior = prior1)##, start=list(QUASI=FALSE))

##save("alevin.bwt.1", file="alevin_fit1.RData")

summary(alevin.bwt.1)

autocorr.plot(alevin.bwt.1$VCV)

acf(alevin.bwt.1$VCV[,5])
## still high for animal and residual components

## plot posterior chain of variance components
par(mar=c(3,3,1,1))
plot(alevin.bwt.1$VCV)

## variance means
(variance.means<-apply(alevin.bwt.1$VCV,2,mean))

## variance proportions
round(variance.means/sum(variance.means),3)

## heritability vector
alevin.heritability <- alevin.bwt.1$VCV[, "animal"]/(alevin.bwt.1$VCV[,"animal"] + alevin.bwt.1$VCV[, "Dam"] + alevin.bwt.1$VCV[, "Tank_ID"] + alevin.bwt.1$VCV[, "ID"]+ alevin.bwt.1$VCV[, "units"])

plot(alevin.heritability)
## note very wide distribution

mean(alevin.heritability)

HPDinterval(alevin.heritability, 0.95)

##-------------------------------------
## MCMCglmm RUN WITHOUT RANDOM EFFECTS
##-------------------------------------

## improper uniform prior on variance components
prior2 <- list(R = list(V = bwt.var/5, nu = 1)) ## keep the same prior


alevin.bwt.2 <- MCMCglmm(body_weight ~ Reciprocal_Cross + degree_days_centered, data = alevin.dat, nitt =nitt, thin = thin, burnin=nitt/2, prior = prior2)

##save("alevin.bwt.2", file="alevin_fit2.RData")

summary(alevin.bwt.2)

autocorr.plot(alevin.bwt.2$VCV)

## plot posterior chain of fixed effects
par(mar=c(3,3,1,1))
plot(alevin.bwt.2$Sol, mfrow=c(9,9))

##-------
## PLOTS 
##-------

library(ggplot2)
library(reshape)

## make a dataframe with the posterior parameters by fitting method
## get cross and degree day effects
post1.0<-cbind(alevin.bwt.1$Sol[,1], alevin.bwt.1$Sol[,1]+alevin.bwt.1$Sol[,-c(1,15)], alevin.bwt.1$Sol[,15])
##
colnames(post1.0)<-c(levels(alevin.dat$Reciprocal_Cross), "Degree.days")
## change to long format
post1.1<-melt(post1.0)
## order so degree day last
post1.1$variable<-factor(as.character(post1.1$X2), levels=c(levels(post1.1$X2)[-5],levels(post1.1$X2)[5]))

## same for model without random effects
post2.0<-cbind(alevin.bwt.2$Sol[,1], alevin.bwt.2$Sol[,1]+alevin.bwt.2$Sol[,-c(1,15)], alevin.bwt.2$Sol[,15])
##
colnames(post2.0)<-c(levels(alevin.dat$Reciprocal_Cross), "Degree.days")
post2.1<-melt(post2.0)
post2.1$variable<-factor(as.character(post2.1$X2), levels=c(levels(post2.1$X2)[-5],levels(post2.1$X2)[5]))

post1.1$model<-"Random effects"
post2.1$model<-"No random effects"

post.all<-rbind(post1.1, post2.1)

## get relatedness matrix (A)
pedigree.summary<-pedigreeStats(Ped=pedigree.dat) ## takes a while
A<-pedigree.summary$Amatrix[alevin.dat$Ind_id,alevin.dat$Ind_id]

n<-100 ## number of individuals to plot
##white2red<-colorRampPalette(c("white","red"))

p3<-ggplot(melt(A[1:n,1:n]), aes(x=X1, y=X2, fill=value)) + geom_tile(colour='lightgrey')+scale_fill_gradient(low="white", high="black")+ labs(x=NULL, y=NULL) +theme(legend.position="none", axis.text.x = element_blank(), axis.ticks=element_blank(),
                    axis.text.y = element_blank(), plot.margin = unit(c(.3,2,.3,2), "cm"))+ annotate("text", x = 5, y = 96, label = "(a)")


## ggplots
## parameters
## cross
p1 <- ggplot(subset(post.all, variable!="Degree.days"), aes(variable, value, fill=model))
p1<-p1+geom_violin(scale="width")+ theme(legend.position="none", plot.margin = unit(c(.1,.1,.1,.1), "cm"), axis.text.x = element_text(angle=90, vjust=0.5))+xlab(NULL)+ylab("Coefficient value")+ scale_fill_manual(values=c("white", "darkgrey"))+ annotate("text", x = 1, y = 7.5, label = "(b)")

## degree days
p2 <- ggplot(subset(post.all, variable=="Degree.days"), aes(variable, value, fill=model))
p2<-p2+geom_violin(scale="area")+ theme(legend.position="none", plot.margin = unit(c(.1,.1,.1,.1), "cm"))+  scale_y_continuous(trans='reverse')+ylab(NULL)+xlab(NULL)+ scale_fill_manual(values=c("white", "darkgrey"))+ annotate("text", x = 0.7, y = 0.01, label = "(c)")

p_save <- list( 'p1'=p1, 'p2'=p2, 'p3'=p3, 'post.all'=post.all, 'A'=A)
  save( p_save, file="p_save.RData")

## display all the plots
jpeg("Fig_5_Alevin_plots.jpg", height=8, width=7, units="in", res=600)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(3, 2, widths=c(0.80,0.20))))
  print(p3, vp = viewport(layout.pos.row = 1:2, layout.pos.col = 1:2))
  print(p1, vp = viewport(layout.pos.row = 3, layout.pos.col = 1, height=3))
  print(p2, vp = viewport(layout.pos.row = 3, layout.pos.col = 2))
dev.off()



