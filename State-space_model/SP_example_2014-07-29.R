

File = paste0( getwd(),"/")

# Load results
INDEX = scan( paste(File,"tuna.dat",sep=""), skip=27, nlines=23)
CATCH = scan( paste(File,"tuna.dat",sep=""), skip=3, nlines=23)
# Separable -- SP_sep_v3
STD_1 = read.table( paste(File,"theta_mod_v7.std",sep=""), header=FALSE, skip=1)
  Which_1a = which(STD_1[,2]=="qU")
# Measurement error -- SP_meas_v1
STD_2 = read.table( paste(File,"theta_mod_meas_v2.std",sep=""), header=FALSE, skip=1)
  Which_2a = which(STD_2[,2]=="qU")

# Compile results
Mat = cbind('Years'=1968:1990, 'Catch'=CATCH,'Index'=INDEX)
  Mat = cbind( Mat, 'Sep_hat'=STD_1[Which_1a,3], 'Sep_upr'=(STD_1[Which_1a,3]+1.96*STD_1[Which_1a,4]), 'Sep_lwr'=(STD_1[Which_1a,3]-1.96*STD_1[Which_1a,4]) )
  Mat = cbind( Mat, 'Meas_hat'=STD_2[Which_2a,3], 'Meas_upr'=(STD_2[Which_2a,3]+1.96*STD_2[Which_2a,4]), 'Meas_lwr'=(STD_2[Which_2a,3]-1.96*STD_2[Which_2a,4]) )

# Plot results
jpeg( file=paste(File,"Fig_2_Estimates.jpeg",sep=""), width=5, height=5, res=600, units="in")
  Ratio = max(Mat[,c('Catch')])/max(Mat[,c('Sep_upr','Meas_upr','Index')])
  par( mar=c(3,3,0,3)+0.5, mgp=c(2,0.5,0), tck=-0.02, xaxs="i", yaxs="i")
  matplot( x=Mat[,'Years'], y=Mat[,c('Catch')], ylim=c(0, max(Mat[,c('Catch')])*1.05), type="l", col="black", lty="solid", lwd=1, xlab="Year", ylab="Relative abundance (kg./100 hooks)", yaxt="n" )
  polygon( x=c(Mat[,'Years'],rev(Mat[,'Years'])), y=c(Mat[,'Sep_upr'],rev(Mat[,'Sep_lwr']))*Ratio, col=rgb(0,0,0,0.4), border=NA )
  polygon( x=c(Mat[,'Years'],rev(Mat[,'Years'])), y=c(Mat[,'Meas_upr'],rev(Mat[,'Meas_lwr']))*Ratio, col=rgb(0,0,0,0.2), border=NA )
  matplot( x=Mat[,'Years'], y=Mat[,c('Index','Sep_hat','Meas_hat')]*Ratio, add=TRUE, type="l", col=c("black","black","black"), lty=c("solid","dotted","dashed"), lwd=c(2.5,2.5,2.5) )
  axis( side=4, at=pretty(c(0,max(Mat[,c('Catch')]))), labels=pretty(c(0,max(Mat[,c('Catch')]))) )
  axis( side=2, at=pretty(c(0,max(Mat[,c('Sep_upr','Meas_upr')])))*Ratio, labels=pretty(c(0,max(Mat[,c('Index','Sep_upr','Meas_upr')]))) )
  mtext( side=4, "Catch (1000 t.)", line=2)
dev.off()
