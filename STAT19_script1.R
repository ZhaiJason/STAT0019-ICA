# Informative prior on length of influenza episodes
## Compute the value of parameters (mulog,sigmalog) for a logNormal 
## distribution to have mean and sd (m,s)
lognPar <- function(m,s) {
  s2 <- s^2
  mulog <- log(m) - 0.5 * log(1+s2/m^2)
  s2log <- log(1+(s2/m^2))
  sigmalog <- sqrt(s2log)
  list(mulog = mulog, sigmalog = sigmalog)
}
# Defines the data
# Number of interventions (t=1: no prophylaxis; t=2: antibiotics prophylaxis) 
T <- 2                  

# Informative prior on incidence of complications in "healthy" adults
m.p1 <- 0.08         # original mean length for patients with complications
s.p1 <- 0.02         # original sd length for patients with complications
a.p1 <-  m.p1*((1-m.p1)*m.p1/s.p1^2-1)
b.p1 <-  (1-m.p1)*((1-m.p1)*m.p1/s.p1^2-1)

# Evidence synthesis on effectiveness of antibiotics prophylaxis vs status quo
r0 <- r1 <- n0 <- n1 <- numeric()   # defines observed cases & sample sizes
r0 <- 37
r1 <- 19
n0 <- 161
n1 <- 166

# Data on costs
unit.cost.comp <- 168.19  # unit (daily) cost of NI
unit.cost.no <- 113.61
unit.e.comp <- 0.0013151
unit.e.no <- 0.0025205
c.anti <- 201.47          # cost of GP visit to administer prophylactic NI

# Informative prior on length of impatient stay
m.l.comp <- 33.66         # original mean length for patients with complications
s.l.comp <- 2.1       # original sd length for patients with complications
mu.l.comp <- lognPar(m.l.comp,s.l.comp)$mulog # mean time to recovery (log scale)
sigma.l.comp <- lognPar(m.l.comp,s.l.comp)$sigmalog # sd time to recovery (log scale)
tau.l.comp <- 1/sigma.l.comp^2    # precision time to recovery (log scale)
#a.gamma.comp <- m.l.comp^2/s.l.comp^2
#b.gamma.comp <- m.l.comp/s.l.comp^2

m.l.no <- 4.36         # original mean length for patients with no complications
s.l.no <- 1.8       # original sd length for patients with no complications
mu.l.no <- lognPar(m.l.no,s.l.no)$mulog # mean time to recovery (log scale)
sigma.l.no <- lognPar(m.l.no,s.l.no)$sigmalog # sd time to recovery (log scale)
tau.l.no <- 1/sigma.l.no^2    # precision time to recovery (log scale)
#a.gamma.no <- m.l.no^2/s.l.no^2
#b.gamma.no <- m.l.no/s.l.no^2

# Parameters of unstructured effects
mean.delta <- 0
sd.delta <- sqrt(10000)
tau.delta <- 1/sd.delta^2
mean.alpha <- 0
sd.alpha <- sqrt(10000)
prec.alpha <- 1/sd.alpha^2

# Prepares to launch OpenBUGS
library(R2OpenBUGS)

# a.gamma.comp = a.gamma.comp, b.gamma.comp = b.gamma.comp, a.gamma.no = a.gamma.no, b.gamma.no = b.gamma.no,
# Creates the data list
data <- list(r0=r0,r1=r1,n0=n0,n1=n1,a.p1=a.p1,b.p1=b.p1,
             mu.l.comp=mu.l.comp,tau.l.comp=tau.l.comp,mu.l.no=mu.l.no,tau.l.no=tau.l.no,
             mean.alpha=mean.alpha,prec.alpha=prec.alpha,
             mean.delta=mean.delta,tau.delta=tau.delta)

# Points to the txt file where the OpenBUGS model is saved
filein <- "code.txt"

# Defines the parameters list
params <- c("p1","p2","l.no","l.comp","rho","alpha","delta")

# Sets the number of iterations, burnin and thinning
n.iter <- 10000
n.burnin <- 1000
n.thin <- 1

# Finally calls OpenBUGS to do the MCMC run and saves results to the object "es"
# model_bugs <- bugs(data=data, inits=NULL,parameters.to.save=params,model.file=filein,
#                   n.chains=2, n.iter, n.burnin, n.thin, DIC=TRUE)

# Displays the summary statistics
#print(model_bugs,digits=3,intervals=c(0.025, 0.975))

library(R2jags)
model_jags <- jags(data=data, inits=NULL,parameters.to.save=params,model.file=filein,
                   n.chains=2, n.iter, n.burnin, n.thin, DIC=TRUE)
print(model_jags, digits=3, intervals = c(0.025, 0.975))

# Convergence check through traceplots (example for node p1)
plot(model_jags$BUGSoutput$sims.list$p1[1:9000],t="l",col="blue",ylab="p1")
points(model_jags$BUGSoutput$sims.list$p1[9001:18000],t="l",col="red")
# Attaches the es object to the R workspace (to use the posteriors for the economic analysis)
attach.jags(model_jags)

# Runs economic analysis 
# cost of treatment
c <- e <- matrix(NA, n.sims, T)
c[,1] <- p1*unit.cost.comp*mu.l.comp + (1-p1)*unit.cost.no*mu.l.no
c[,2] <- c.anti + p2*unit.cost.comp*mu.l.comp + (1-p2)*unit.cost.no*mu.l.no
e[,1] <- p1*(365-mu.l.comp + mu.l.comp*unit.e.comp) + (1-p1)*(365-mu.l.no + mu.l.no*unit.e.no)
e[,2] <- p2*(365-mu.l.comp + mu.l.comp*unit.e.comp) + (1-p2)*(365-mu.l.no + mu.l.no*unit.e.no)

library(BCEA)
treats <- c("status quo","prophylatic antibiotics")
m <- bcea(e,c,ref=2,treats,Kmax=10000, plot = TRUE)
ceplane.plot(m,wtp=10000)

#NBt <- function(a,b,k) {
#  return NB1 <- k*mu.e1-mu.c1
#  NB2 <- k*mu.e2-mu.c2
#}