########################################################################
#STEP 1:                                                               #
#CONSRUCTION OF THE EHFh TO DEFINE A LOCALIZED HEATWAVE                #
########################################################################

# CREATE THE MATRIX TO STORE THE LOCALIZED THRESHOLD (Th)
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
  
  # 1. USE TIERED HEALTH RISK-BASED (THR) APPROACH TO CALCULATE THE LOCALIZED THRESHOLD (Th)
  nyear<-unique(data1$year)
  v<-vector("list",length(nyear))
  for (j in seq(length(nyear))) {
    a1<-subset(data1, data1$year == nyear[j])
    hwdata<-matrix(NA,length(a1$tm),n.hw)
    for (s in seq(n.hw)) {
      q_th<-seq(0.75,0.99,0.005)
      hwdata[,s]<-per_hw(a1$tm,2,q_th[s])
    }
    v[[j]]<-hwdata
  }
  hwall<-do.call("rbind",v)
  colnames(hwall)<-paste0("hw_",seq(75,99,0.5),"p","_2d")
  data2<-cbind(data1,hwall)
  
  # EXTRACT DATA IN WARM SEASON
  data3 <- subset(data2,month(data2$deathdate) %in% c(5,6,7,8,9))
  data3$dow <- format(data3$deathdate,"%w")
  data3$time <- seq(1,length(data3$deathdate),1)
  
  # CREATE THE MATRIX TO STORE THE EFFECTS (RR) of HEATWAVES (75TH-99TH)----------
  cityhw<-matrix(NA,n.hw,3,dimnames=list(1:n.hw,c(paste0("Seoul","_rr"),paste0("Seoul","_rrlow"),paste0("Seoul","_rrhigh"))))
  for (j in seq(n.hw)) {
    cb <- crossbasis(data3[,(j+7)], lag=lag, 
                     argvar=list(fun="lin"), 
                     arglag=list(fun="ns",df=3))
    model <- gam(count ~  cb + s(rh,bs="cr",k=3)+
                   s(time,bs="cr",k=ti*n)+
                   as.factor(dow)+as.factor(holiday),
                 method="GCV.Cp",family=poisson(), data=data3)
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
  # X AXIS: HEATWAVE THRESHOLD. Y AXIS: RR 
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
  
  # 2. CALCULARE EHI_h.sig, EHI_accl, AND EHFh----------------------------------
  nyear<-unique(data1$year)
  v1<-vector("list",length(nyear))
  qpoint<-Th.threshold[t,1]/100
  for (j in 1:length(nyear)) {
    a1<-subset(data1,year(data1$deathdate) %in% nyear[j])
    Th<-round(quantile(a1$tm,qpoint,na.rm=T),1)
    v1[[j]]<-as.matrix(calc_EHI_h.sig(a1$tm, Th))
  }
  EHI_hsig <- do.call("rbind", v1)
  c <- cbind(data1, EHI_hsig)
  c$EHI_accl <- calc_EHI_accl(c$tm, 30, 3)
  c$EHFh <- round(calc_EHFh(c$EHI_hsig, c$EHI_accl), 1)
  
  ##############################################################################                                                              #
  # STEP 2:                                                                    #
  # ASSESSMENT OF EHFh                                                         #
  ##############################################################################

  # EXTRACT DATA IN WARM SEASON
  c1 <- subset(c, month(c$deathdate) %in% c(5,6,7,8,9))
  c1$dow <- format(c1$deathdate, "%w")
  c1$time <- seq(1, length(c1$deathdate), 1)
  round(range(c1$EHFh, na.rm=T), 1)
  
  # 1. ASSESSMENT OF HEATWAVE CHARACTERISTICS----------------------------------
  c2 <- subset(c1, c1$EHFh>0)
  c2$diff <- 1
  c2$num <- 1:length(c2$tm)
  c2$diff[2:length(c2$deathdate)] <- diff(c2$deathdate)
  o <- as.numeric(which(c2$diff != 1))
  w <- matrix(1,(length(o)+1),1)
  w[2:(length(o)+1)] <- o
  as[t,1] <- round(length(w)/n,1)
  as[t,2] <- round(length(c2$EHFh)/n,1)
  as[t,3] <- max(diff(w))
  maxi <- which.max(c2$EHFh)
  as[t,4] <- round(c2$tm[maxi],1)
  as[t,5] <- round(mean(c2$tm),1)
  
  # 2. MORTALITY RISK-----------------------------------------------------------
  # coef IS THE MATRIX FOR THE OUTCOME PARAMETERS
  # vcov IS THE LIST WITH (CO)VARIANCE MATRICES
  coef <- matrix(NA,t,e)
  vcov <- vector("list",t)
  
  # DEFINE THE CROSSBASIS
  argvar <- list(fun="ns",df=e)
  cb <- crossbasis(c1$EHFh,lag=lag,argvar=argvar,arglag=list(fun="ns",df=3))
  
  # RUN THE MODEL AND OBTAIN PREDICTIONS
  model <- gam(count ~  cb + s(time,bs="cr",k=ti*n) + s(rh,bs="cr",k=3)
               +as.factor(dow) + as.factor(holiday),
               method="GCV.Cp",family=quasipoisson(), data=c1)
  pred <- crosspred(cb, model, cen=cen, by=1, cumul=T, from=0, to=max(c1$EHFh,na.rm=T))
  
  # X AXIS: EHFh. Y AXIS: RR
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
  # EHFh > 0 REPRESENT A HEATWAVE
  ad <- subset(c1, EHFh > 0)
  matsim <- attrdl(ad$EHFh,cb,ad$count,model,type="an",dir="forw",cen=0)
  arraysim[1,"an",] <- attrdl(ad$EHFh,cb,ad$count,model,type="an",dir="forw",cen=0,sim=T,nsim=nsim)
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

# LOCALIZED THRESHOLDS FOR SEOUL
Th.threshold
round(apply(Th.threshold,2,mean,na.rm=T),1)

# ASSESSMENT OF EHFh
as
round(apply(as,2,mean,na.rm=T),1)

