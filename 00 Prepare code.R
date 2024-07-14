################################################################################
# PREPARE THE FUNCTION                                                         #
################################################################################

# PERCENTILE-BASED HEATWAVE FUNCTION-------------------------------------------
per_hw <- function(v,days,q) { #v is temperature; days is heatwaves days;q is quantile;
  q.new<-quantile(v,q,na.rm =T)
  s.hw<-rep(0,length(v))
  s.hw[v>=q.new]<-1
  s.hw.lag<-Lag(s.hw,c(0:(days-1)))
  s.hw.sum<-apply(s.hw.lag,1,sum)
  s.hw.sum[s.hw.sum<days]<-0
  for (i in (days:length(s.hw.sum))){
    if (s.hw.sum[i]==days){
      s.hw.sum[c((i-days+1):(i-1))]<-1
    }
  }
  s.hw.sum[s.hw.sum>0]<-1
  hw<-s.hw.sum
  return(hw)
}

# EHFh FUNCTION-----------------------------------------------------------------
calc_EHI_h.sig <- function(temps, Th_threshold, sig_spread = 3) {
  zoo::rollapply(temps, width = sig_spread, FUN = function(x) { sum(x)/sig_spread }, align = 'right', fill = NA) - Th_threshold
}

calc_EHI_accl <- function(temps, accl_length = 30, sig_spread = 3) {
  # figures to average for acclimatisation period
  width_ <- list( -(accl_length+sig_spread-1) : -sig_spread )
  
  # first 3 days: i, i-1, i-2
  part1 <- zoo::rollapply(temps, width = sig_spread, FUN = function(x) { sum(x)/sig_spread }, align = 'right', fill = NA)
  
  # next 30 days after that: i-3, i-4, ..., i-32
  part2 <- zoo::rollapply(temps, width = width_, FUN = function(x) { sum(x)/accl_length }, align = 'right', fill = NA)
  
  part1 - part2
}

calc_EHFh <- function(EHI_h.sig, EHI_accl) {
  # EHFh = EHI_h.sig * max(1, EHI_accl)
  EHI_h.sig * ifelse(EHI_accl < 1, 1, EHI_accl)
}


################################################################################
# PREPARE THE DATA                                                         #
################################################################################

# LOAD PACKAGES
library(dlnm);library(mvmeta);library(tsModel);library(lubridate)
library(splines);library(INLA);library(mgcv)
source("attrdl.R")
# LOAD THE DATA
# NOTE: THE SEOUL DATA WILL BE OPEN FOR 2 YEARS.
data0<-readRDS("SouthKoreaSeoul.rds")

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

## LAYOUT FIGURES---------------------------------------------------------------
par(mfrow = c(6,6),mar = c(2.0,3.0,1.5,1.0), ps=13)

