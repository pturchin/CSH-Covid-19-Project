### Numerical solution of the SIRD model: Korea
{#########  Solve the SIRD model for given parameter values
  dat <- read.table('CovidDat.csv', sep=",", header=TRUE, stringsAsFactors = 	FALSE)
  dur <- nrow(dat)
  
  ### Initial conditions and parameter valuse
  N <- 50000000  # total population = 50 mln (S. Korea)
  S0 <- N
  I0 <- 0.033
  D0 <- R0 <- 0
  beta0 <- 0.432
  theta <- 0.234
  gamma <- 0.0078
  delta <- 0.00044
  q1 <- -0.36
  q2 <- 0.68
  b_date <- 36.1
  q_date <- 33
  theta_q <- 0.3
  
  ### Calculate trajectory
  model <- data.frame(Time = 1:dur, S=S0, I=I0, R=R0, D=D0)
  pred <- data.frame(Time = 1:dur, I=NA, R=NA, D=NA, C=NA, delC=NA,delD=NA)
  
  for(t in 1:(dur-1)){
    beta <- beta0*(1-1/(1+exp(-theta*(t-b_date))))
    model$S[t+1] <- model$S[t] - beta*model$S[t]*model$I[t]/N
    model$I[t+1] <- model$I[t] + beta*model$S[t]*model$I[t]/N - (gamma+delta)*model$I[t]
    model$R[t+1] <- model$R[t] + gamma*model$I[t]
    model$D[t+1] <- model$D[t] + delta*model$I[t]
  }
  
  ### Data prediction
  pred <- data.frame(Time = 1:dur, I=NA, R=NA, D=NA, C=NA, delC=NA,delD=NA)
  for(t in 1:(dur)){
    q <- q1 + (q2-q1)/(1+exp(-theta_q*(t-q_date)))
    pred$I[t] <- q*model$I[t]
    pred$R[t] <- q*model$R[t]
    pred$D[t] <- 1*model$D[t]
  }
  pred$C <- pred$I + pred$R + pred$D
  pred$delC[2:dur] <- pred$C[2:dur] - pred$C[1:(dur-1)]
  pred$delD[2:dur] <- pred$D[2:dur] - pred$D[1:(dur-1)]
  
  ### Plot results and calculate the fit
  par(mfrow=c(3,2))
  plot(dat$Time,dat$I, type="p", pch=16, col="dark green", main="Infected (Active Cases)", xlab="",ylab="")
  lines(pred$Time, pred$I, lwd=3, col="dark green")
  plot(dat$Time,dat$D, type="p", pch=16, col="dark red", main="Total Deaths", xlab="",ylab="")
  lines(pred$Time, pred$D, lwd=3, col="dark red")
  plot(dat$Time,dat$delC, type="p", pch=16, col="dark green", main="New Cases", xlab="",ylab="")
  lines(pred$Time, pred$delC, lwd=3, col="dark green")
  plot(dat$Time,dat$delD, type="p", pch=16, col="dark red", main="Daily Deaths", xlab="",ylab="")
  lines(pred$Time, pred$delD, lwd=3, col="dark red")
  plot(dat$Time,dat$C, type="p", pch=16, col="blue", main="Total Confirmed Cases", xlab="Days from Jan. 22",ylab="")
  lines(pred$Time, pred$C, lwd=3, col="blue")
  plot(dat$Time,dat$R, type="p", pch=16, col="blue", main="Recovered", xlab="Days from Jan. 22",ylab="")
  lines(pred$Time, pred$R, lwd=3, col="blue")
  par(mfrow=c(1,1))
  
  predR2 <- data.frame(Var=c("I","D","delC","delD","C","R","Mean"), R2=NA)
  pred <- pred[2:dur,]
  dat <- dat[2:dur,]
  predR2$R2[1] <- 1 - sum((pred$I - dat$I)^2)/sum((dat$I - mean(dat$I))^2)
  predR2$R2[2] <- 1 - sum((pred$D - dat$D)^2)/sum((dat$D - mean(dat$D))^2)
  predR2$R2[3] <- 1 - sum((pred$delC - dat$delC)^2)/sum((dat$delC - mean(dat$delC))^2)
  predR2$R2[4] <- 1 - sum((pred$delD - dat$delD)^2)/sum((dat$delD - mean(dat$delD))^2)
  predR2$R2[5] <- 1 - sum((pred$C - dat$C)^2)/sum((dat$C - mean(dat$C))^2)
  predR2$R2[6] <- 1 - sum((pred$R - dat$R)^2)/sum((dat$R - mean(dat$R))^2)
  predR2$R2[7] <- mean(predR2$R2[1:6])
  print(predR2, digits=3)
}

{ ### Plot b(t) and q(t)
  par(mfrow=c(1,2))
  t <- 1:dur
  beta <- beta0*(1-1/(1+exp(-theta*(t-b_date))))
  plot(t,beta, type="l", lwd=3, col="dark red", main="(a)")
  abline(v=seq(0,60, by=10), h=seq(0,1,by=0.05), col="grey")
  
  q <- q1 + (q2-q1)/(1+exp(-theta_q*(t-q_date)))
  q[q < 0] <- 0
  plot(t,q, type="l", lwd=3, col="dark red", main="(b)")
  abline(v=seq(0,60, by=10), h=seq(0,1,by=0.2), col="grey")
  par(mfrow=c(1,1))
}

 
 