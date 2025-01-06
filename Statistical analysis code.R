################################################################################
# PREPARE THE DATA                                                         #
################################################################################

# LOAD PACKAGES
library(dlnm);library(mvmeta);library(tsModel);library(lubridate);library(splines);library(INLA);library(mgcv)
source("attrdl.R")
# LOAD THE DATA
# NOTE: THE SEOUL DATA WILL BE OPEN FOR 2 YEARS.
data0<-readRDS("deathdata.RDS")

## COEFFICIENTS SUMMARY---------------------------------------------------------
n.hw<-49 #49 heatwave thresholds (i.e., 75th, 75.5th, …, 99th)
bc <- length(seq(as.Date("1997-01-01"),as.Date("1997-12-31"),1)) # Running steps 
n <- 5 #Years of training data
length <- bc*n #Length of training data
cs <- floor((length(data0$deathdate)-length)/bc)+1 #Loop times
leapyear <- c(1,0,0,0,1,0,0,0,1,0,0,0,1,0,0,0,1,0) #Is this a year after a leap year (1997-2018)
e <- 3 #df of cb
ti <- 3 #df of time
lag <- 0 #Lag days
cen <- 0 #Centering point

# CREATE THE MATRIX TO STORE THE HEALTH-BASED THRESHOLD (Th)
Th.threshold<-matrix(NA,cs,2,dimnames=list(1:cs,c("Th.percentile","Th.temperature")))

# CREATE THE MATRIX TO STORE THE RESULTS OF ASSESSMENT VALUES
as<- matrix(NA,18,6,dimnames=list(c(1:18),c("HWN","HWF","HWD","HWA","HWM","AF(%)")))

for (t in 1:cs){
  # PRINT ITERATION
  cat(t,"")
  
  # LOAD DATA
  start=(t-1)*bc+sum(leapyear[1:t]) 
  end=length+(t-1)*bc+sum(leapyear[1:t])
  data1<-data0[start:end,] 
  
  # EXTRACT DATA IN WARM SEASON
  data2 <- subset(data1,month(data1$deathdate) %in% 5:9)
  data2$dow <- format(data2$deathdate,"%w")
  data2$time <- seq(1,length(data2$deathdate),1)
  
  # CREATE THE MATRIX TO STORE THE EFFECTS (RR) of HEATWAVES (75TH-99TH)----------
  cityhw<-matrix(NA,n.hw,3,dimnames=list(1:n.hw,c(paste0("Seoul","_rr"),paste0("Seoul","_rrlow"),paste0("Seoul","_rrhigh"))))
  for (j in seq(n.hw)) {
    cb <- crossbasis(data2[,(j+7)], lag=lag, 
                     argvar=list(fun="lin"), 
                     arglag=list(fun="ns",df=3))
    model <- gam(count ~  cb + s(rh,bs="cr",k=3)+
                   s(time,bs="cr",k=ti*n)+
                   as.factor(dow)+as.factor(holiday),
                 method="GCV.Cp",family=poisson(), data=data2)
    pred <- crossreduce(cb,model,type="overall",cen=0,by=1)
    cityhw[j,]<-as.matrix(cbind(pred$RRfit,pred$RRlow,pred$RRhigh))[2,]
  }
  
  # USE BAYESIAN HIERARCHICAL MODEL TO IDENTIFY THE LOCALIZED THRESHOLD
  h<-cityhw
  yl<-min(h[,2]);yh<-max(h[,3])
  ds<-data.frame(y=seq(n.hw*3),x=c(seq(1:49),seq(50:98),seq(99:147)),group=c(rep(1,n.hw),rep(2,n.hw),rep(3,n.hw)))
  ds$y[1:49]<-h[,1]
  ds$y[50:98]<-h[,2]
  ds$y[99:147]<-h[,3]
  datadic<-data.frame(x=seq(n.hw),quantile=seq(75,99,0.5),dic=seq(n.hw))
  for (j in seq(n.hw)) {
    formula <- y ~ 1 + bs(ds$x,knots=c(j),Boundary.knots=c(2,(n.hw-1))) + f(group, model = "iid")
    output <- inla(formula,  family = "gaussian", 
                   control.compute=list(dic=TRUE),data=ds)
    datadic[j,3]<-output$dic$dic
  }
  # IDENTIFY TURNING POINT OF RR
  # NOTE: WE ASSUME THE TURNING POINT RANGES FROM 91ST TO 98TH 
  datadic<-datadic[order(datadic[,3]),]
  formula1 <- y ~ 1 + bs(ds$x,knots=datadic$x[1]) + f(group, model = "iid")
  output1 <- inla(formula1,  family = "gaussian", 
                  control.compute=list(dic=TRUE),data=ds)
  hposter0<-cbind(output1$summary.fitted.values$mean[1:49],
                  output1$summary.fitted.values$mean[50:98],
                  output1$summary.fitted.values$mean[99:147])
  hposter <- cbind(h,hposter0)
  p<-datadic[1,2]
  p<-ifelse(p>=91,p,91)
  p<-ifelse(p>=98.5,98,p)
  
  # STORE THE LOCALIZED THRESHOLD (Th)
  Th.threshold[t,1]<-p
  Th.threshold[t,2]<-round(quantile(data1$tm,p/100,na.rm =T),1)
  
  # PLOT TURNING POINT OF HEATWAVE_75TH-99TH
  # X AXIS: TEMPERATURE THRESHOLD. Y AXIS: RR 
  main<-paste0("Seoul",t)
  xlim<-seq(75,99,0.5)
  plot(xlim,h[,1],type="p",ylim=c(0.9,yh),pch=19,cex=0.9,col="black",
       las=1,xlab="",ylab="",main=main,cex.main=1.3,cex.axis=1.2)
  abline(h=1)
  for (j in seq(n.hw)) {
    lines(c(xlim[j],xlim[j]),c(h[j,2],h[j,3]),type="l")
  }
  lines(xlim,hposter[,4],type="l",col="red")
  lines(xlim,hposter[,5],type="l",col="blue")
  lines(xlim,hposter[,6],type="l",col="blue")
  abline(v=p,lty=3)

  ##############################################################################                                                              
  # ASSESSMENT OF HEHF                                                         #
  ##############################################################################

  # EXTRACT DATA IN WARM SEASON
  c1 <- subset(data1, month(data1$deathdate) %in% 5:9)
  c1$dow <- format(c1$deathdate, "%w")
  c1$time <- seq(1, length(c1$deathdate), 1)
  round(range(c1$HEHF_, na.rm=T), 1)
  
  # 1. ASSESSMENT OF HEATWAVE CHARACTERISTICS----------------------------------
  c2 <- subset(c1, c1$HEHF_>0)
  c2$diff <- 1
  c2$num <- 1:length(c2$tm)
  c2$diff[2:length(c2$deathdate)] <- diff(c2$deathdate)
  o <- as.numeric(which(c2$diff != 1))
  w <- matrix(1,(length(o)+1),1)
  w[2:(length(o)+1)] <- o
  as[t,1] <- round(length(w)/n,1)
  as[t,2] <- round(length(c2$HEHF_)/n,1)
  as[t,3] <- max(diff(w))
  maxi <- which.max(c2$HEHF_)
  as[t,4] <- round(c2$tm[maxi],1)
  as[t,5] <- round(mean(c2$tm),1)
  
  # 2. MORTALITY RISK-----------------------------------------------------------
  # coef IS THE MATRIX FOR THE OUTCOME PARAMETERS
  # vcov IS THE LIST WITH (CO)VARIANCE MATRICES
  coef <- matrix(NA,t,e)
  vcov <- vector("list",t)
  
  # DEFINE THE CROSSBASIS
  argvar <- list(fun="ns",df=e)
  cb <- crossbasis(c1$HEHF_,lag=lag,argvar=argvar,arglag=list(fun="ns",df=3))
  
  # RUN THE MODEL AND OBTAIN PREDICTIONS
  model <- gam(count ~  cb + s(time,bs="cr",k=ti*n) + s(rh,bs="cr",k=3)
               +as.factor(dow) + as.factor(holiday),
               method="GCV.Cp",family=quasipoisson(), data=c1)
  pred <- crosspred(cb, model, cen=cen, by=1, cumul=T, from=0, to=max(c1$HEHF_,na.rm=T))
  
  # X AXIS: HEHF. Y AXIS: RR
  plot(pred,"overall",col="#CC3333",lwd=3,las=1,cex.axis=1.2,bty="o",xlab="",ylab="",
       ci.arg=list(col="#FFCCCC"), main=main,cex.main=1.3,ylim=c(0.9,1.5))
  
  # REDUCTION TO OVERALL CUMULATIVE
  red <- crossreduce(cb, model, cen=cen)
  coef[t,] <- coef(red)
  vcov[[t]] <- red$vcov
  
  # 3. ATTRIBUTABLE FRACTION----------------------------------------------------
  # NUMBER OF SIMULATION RUNS FOR COMPUTING EMPIRICAL CI
  nsim <- 1000
  # CREATE THE ARRAY TO STORE THE CI OF ATTRIBUTABLE DEATHS
  arraysim <- array(NA,dim=c(1,1,nsim),dimnames=list("Seoul",c("an")))
  # HEHF > 0 REPRESENT A HEATWAVE
  ad <- subset(c1, HEHF_ > 0)
  matsim <- attrdl(ad$HEHF_,cb,ad$count,model,type="an",dir="forw",cen=0)
  arraysim[1,"an",] <- attrdl(ad$HEHF_,cb,ad$count,model,type="an",dir="forw",cen=0,sim=T,nsim=nsim)
  totdeath <- sum(ad$count,na.rm=T)
  
  # ATTRIBUTABLE NUMBERS
  ancity <- round(matsim,1)
  ancitylow <- round(apply(arraysim,c(1,2),quantile,0.025),1)
  ancityhigh <- round(apply(arraysim,c(1,2),quantile,0.975),1)
  
  # ATTRIBUTABLE FRACTIONS
  afcity <- round(ancity/totdeath*100,1)
  afcitylow <- round(ancitylow/totdeath*100,1)
  afcityhigh <- round(ancityhigh/totdeath*100,1)
  as[t,6] <- afcity
  
}

# HEALTH-BASED THRESHOLDS FOR SEOUL
Th.threshold
round(apply(Th.threshold,2,mean,na.rm=T),1)

# ASSESSMENT OF HEHF
as
round(apply(as,2,mean,na.rm=T),1)

