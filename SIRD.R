### Numerical solution of the SIRD model
### t=0 at Feb 15; t=30 at March 15
### Initial conditions
dur <- 30 # Duration of simulation in days
N <- 50000000  # total population = 50 mln 
S0 <- N
I0 <- 10
D0 <- R0 <- 0

### Parameters
beta <- seq(0.6,0, by=-0.02)  ### beta declines linearly
gamma <- 0.1
delta <- 0.01
q <- 0.5
r <- 1

### Calculate trajectory
dat <- data.frame(Time = 1:dur, S=S0, I=I0, R=R0, D=D0, delI=0,delD=0)
for(t in 1:(dur-1)){
  dat$S[t+1] <- dat$S[t] - beta[t]*dat$S[t]*dat$I[t]/N
  dat$I[t+1] <- dat$I[t] + beta[t]*dat$S[t]*dat$I[t]/N - (gamma+delta)*dat$I[t]
  dat$R[t+1] <- dat$R[t] + gamma*dat$I[t]
  dat$D[t+1] <- dat$D[t] + delta*dat$I[t]
  dat$delI[t+1] <- beta[t]*dat$S[t]*dat$I[t]/N
  dat$delD[t+1] <- dat$D[t+1] - dat$D[t]
}
dat$I <- q*dat$I
dat$delI <- q* dat$delI

### Plot results
par(mfrow=c(2,2))
plot(dat$Time,dat$I, type="l", lwd=3, col="dark green", main="Active Cases")
plot(dat$Time,dat$D, type="l", lwd=3, col="dark red", main="Total Deaths")

plot(dat$Time,dat$delI, type="l", lwd=3, col="dark green", main="New Cases")
plot(dat$Time,dat$delD, type="l", lwd=3, col="dark red", main="Daily Deaths")

par(mfrow=c(1,1))
