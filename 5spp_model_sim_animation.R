### EDM explanation ###
WD <- getwd()
params <- list(nu1=0.1,nu2=0.07,lambda1=3.2,lambda2=2.9,C1star=0.5,C2star=0.5,
               mu1=0.15,mu2=0.15,kappa1=2.5,kappa2=2.0,Rstar=0.3,k=1.2)
dR <- function(R,C1,C2,P1,P2) R*(1-R/params$k) -
  params$mu1*params$kappa1*(C1*R)/(R+params$Rstar) -
  params$mu2*params$kappa2*(C2*R)/(R+params$Rstar)
dC1 <- function(R,C1,C2,P1,P2) params$mu1*params$kappa1*(C1*R)/(R+params$Rstar) -
  params$nu1*params$lambda1*(P1*C1)/(C1+params$C1star) - params$mu1*C1
dC2 <- function(R,C1,C2,P1,P2) params$mu2*params$kappa2*(C2*R)/(R+params$Rstar) -
  params$nu2*params$lambda2*(P2*C2)/(C2+params$C2star) - params$mu2*C2
dP1 <- function(R,C1,C2,P1,P2) params$nu1*params$lambda1*(P1*C1)/(C1+params$C1star) -
  params$nu1*P1
dP2 <- function(R,C1,C2,P1,P2) params$nu2*params$lambda2*(P2*C2)/(C2+params$C2star) -
  params$nu2*P2
dF <- function(R,C1,C2,P1,P2) c(dR(R,C1,C2,P1,P2), dC1(R,C1,C2,P1,P2),
                                dC2(R,C1,C2,P1,P2), dP1(R,C1,C2,P1,P2), dP2(R,C1,C2,P1,P2))

tmax <- 10000; tau <- 5; dt <- 1/100; burn <- 200
d <- array(data = 0, dim = c(tmax/tau,5))
colnames(d) <-c("R","C1","C2","P1","P2")

Xi <- c(1,.5,.8,.7,.8)
#>>>RUN BURN<<<
idex <- 0
while(idex< burn*(1/dt)){
  idex <- idex+1;
  Xi <- Xi + dt*do.call(dF,as.list(Xi))
}
#>>>RUN MODEL<<<
d[1,] <- Xi
idex <- 0; tdex <- 2
while(idex<tmax*(1/dt)){
  idex <- idex+1;
  Xi <- Xi + dt*do.call(dF,as.list(Xi))
  if( (idex*dt) %% tau == 0 && tdex < tmax/tau){
    d[tdex,] <- Xi
    tdex <- tdex+1
  }
}

par(mfrow=c(1,1),mar=c(2,1,1,1)+0.1)

p<-scatterplot3d(d[500:600,c("C1","R","C2")],type="l",color="gray70",xlab="",ylab="",zlab="",scale.y=0.8,grid=F,angle=140,box=F)
p$points3d(d[600,"C1"],d[600,"R"],d[600,"C2"],pch=19,col="red")

par(mar=c(5,4,4,2)+0.1)

## Save animation
frames = 200

for(i in 1:frames){
  # creating a name for each plot file with leading zeros
  if (i < 10) {name = paste(WD,'/data/animate/000',i,'plot.png',sep='')}
  
  if (i < 100 && i >= 10) {name = paste(WD,'/data/animate/00',i,'plot.png', sep='')}
  if (i >= 100) {name = paste(WD,'/data/animate/0', i,'plot.png', sep='')}
  
  
  #saves the plot as a .png file in the working directory
  png(name,width=800,height=400)
  p<-scatterplot3d(d[500:(600+i),c("C1","R","C2")],type="l",color="gray70",xlab="",ylab="",zlab="",scale.y=0.8,grid=F,angle=140,box=F)
  p$points3d(d[(600+i),"C1"],d[(600+i),"R"],d[(600+i),"C2"],pch=19,col="red")

  dev.off()
}

for(i in 1:frames){
  # creating a name for each plot file with leading zeros
  if (i < 10) {name = paste(WD,'/data/animate/univariate/000',i,'plot.png',sep='')}
  
  if (i < 100 && i >= 10) {name = paste(WD,'/data/animate/univariate/00',i,'plot.png', sep='')}
  if (i >= 100) {name = paste(WD,'/data/animate/univariate/0', i,'plot.png', sep='')}
  
  
  #saves the plot as a .png file in the working directory
  png(name,width=800,height=400)
  plot((500+i):(600+i),d[(500+i):(600+i),"C1"],type="l",ylim=c(0,2),xaxt="n",xlab="",ylab="",bty="n",lwd=2,xaxs="i",yaxs="i",col="blue")
  lines((500+i):(600+i),d[(500+i):(600+i),"C2"],type="l",ylim=c(0,2),xlab="",ylab="",lwd=2,xaxs="i",yaxs="i",col="red")
  lines((500+i):(600+i),d[(500+i):(600+i),"R"],type="l",ylim=c(0,2),xlab="",ylab="",lwd=2,xaxs="i",yaxs="i",col="darkgreen")
  
  dev.off()
}


