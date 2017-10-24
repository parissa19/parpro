#######16.2
#a
epsilon<-rnorm(48,0,100)
epsilon

#b
Time<-1:48
dat<-100+10*sin(2*pi*Time/12)-20*cos(2*pi*Time/12)+eps
dat
dat.ts <-ts(dat , frequency=4)
dat.ts
#c
plot.ts(dat.ts,main="time series of simulated data")

#d
decompose(dat.ts)
plot(decompose(dat.ts))

newdat<-dat.ts[-c(37,38,39,40,41,42,43,44,45,46,47,48)]
newdat



DEC.ADD<-function(DATA, s, l)
{
  n <- length(DATA); W <- rep(NA, times = n); D <- rep(NA, times = n)
  s0 <- rep(0, times = n); m <- rep(0, times = s); y <- 0
  forecast <- 0; shat <- rep(0, times = n); seaseffect <- 0
  
  Time <- seq(1:n)
  q <- floor(s/2)
  
  if(s %% 2 == 0)
    for(i in (q + 1):(n - q)) {
      W[i] <- (0.5 * (DATA[i - q] + DATA[i + q]) + sum(DATA[(i - q + 1):(i + q - 1)]))/s		}
  
  else 
    
    for(i in (q + 1):(n - q)) {
      W[i] <- mean(DATA[(i - q):(i + q)])}
  
  for(i in (q + 1):(n - q)) {
    D[i] <- DATA[i] - W[i]
    if(i %% s == 0) {
      s0[s] <- s0[s] + D[i]
      m[s] <- m[s] + 1		}
    
    else 
      for(j in 1:(s - 1)) {
        if(i %% s == j) {
          s0[j] <- s0[j] + D[i]
          m[j] <- m[j] + 1	}
      }
  }
  for(j in 1:s) {
    seaseffect[j] <- s0[j]/m[j] - mean(s0[1:s]/m)	}
  
  shat <- rep(seaseffect, length.out = n)
  y <- DATA - shat
  betahat <- lm(y ~ Time)$coefficients
  yhat <- predict(lm(y ~ Time))
  forecast <-as.vector( yhat + shat)
  error <- as.vector(DATA - forecast)
  FORECAST <- 0
  for(i in 1:l) {
    if((n + i) %% s == 0) {
      FORECAST[i] <- (betahat[1] + betahat[2] * (n + i)) + seaseffect[s]
    }
    else FORECAST[i] <- (betahat[1] + betahat[2] * (n + i)) + seaseffect[(n + i) %% s]
  }
  DECAD <- data.frame(DATA, W, D, shat, y, yhat, forecast	, error)
  return(list(DECAD=DECAD, betahat=betahat, seaseffect=seaseffect,forecast=forecast, FORECAST=FORECAST,error=error))
}
result1<-DEC.ADD(newdat,4,12)
result1


####FORECASTING 12 observation
###Error of Classic  
f<-result1$FORECAST
f

d<-dat
d

Ert<-vector("numeric",12L)
for(i in 1:12)
Ert[i]<-(f[i]-d[i+36])^2
Ert

s1=sum(Ert)


vec<-c(result1$FORECAST,result1$forecast)
vec
lt<-as.list(vec)
lt

par(mfrow=c(1,1))
plot.ts(newdat,ylab="Simulation data",xlab="Time",type="o",xlim=c(1,100))
for(i in 37:48){points(i,lt[i],col="red")}
lines(37:48,lt[37:48],lwd=0.1,col="blue")
title(main="Additive TS Model")


###########################
Time<-1:36
V1<-sin(2*pi*Time/12)
V2<-cos(2*pi*Time/12)
model<-lm(newdat~V1+V2)
model
####################
Time<-37:48
V1<-sin(2*pi*Time/12)
V2<-cos(2*pi*Time/12)
reg.fcast<-model$coefficients[1]+model$coefficients[2]*V1+model$coefficients[3]*V2
reg.fcast


###Error of regression fit
EReg<-vector("numeric",12L)
for(i in 1:12)
EReg[i]<-(reg.fcast[i]-d[36+i])^2
EReg

sum(EReg)
