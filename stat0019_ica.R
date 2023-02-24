# Data Input ===================================================================

# Setting working directory to current file location
setwd(getwd())

# Read file and inspect
rct_data <- read.csv("rct_data.csv")
# head(rct_data)

# Break down by arms of treatment
rct.1 <- rct_data$x[rct_data$T == 0]
rct.2 <- rct_data$x[rct_data$T == 1]
rm(rct_data)

# Gibbs Sampling ===============================================================

## Model Setup information -----------------------------------------------------

# Placebo-controlled effectiveness study data
y1 <- sum(rct.1)
y2 <- sum(rct.2)
n1 <- length(rct.1)
n2 <- length(rct.2)
rm(rct.1, rct.2)

# Prior of population rate of complication p0 (Beta with specified parameters)
mu.p1 <- 0.08 # (0.04 + 0.12) / 2
sd.p1 <- 0.015 # Approximate value sd, so the specified range is within +-3 sd
shape1.p1 <- round(mu.p1 * ((1 - mu.p1) * mu.p1 / sd.p1^2 - 1)) # alpha
shape2.p1 <- round((1 - mu.p1) * ((1 - mu.p1) * mu.p1 / sd.p1^2 - 1)) # beta
rm(mu.p1, sd.p1)

# Plot for specified Beta(26,300)
# x <- seq(0,0.2,0.001)
# plot(x, dbeta(x, shape1.p1, shape2.p1), type = "l")
# abline(v = c(0.04, 0.08, 0.12))
# pbeta(0.12, shape1.p1, shape2.p1) - pbeta(0.04, shape1.p1, shape2.p1)
# pbeta(0.12, 14.64, 168.36) - pbeta(0.04, 14.64, 168.36)

# Minimally informative prior of alpha and delta (large variance normal prior)
mu.alpha <- 0
mu.delta <- 0
tau.alpha <- 1e-2
tau.delta <- 1e-2

# Length of staying (using Gamma with specified parameters)
mu.l0 <- 4.36
mu.l1 <- 33.66
sd.l0 <- 1.8
sd.l1 <- 2.1
alpha.l0 <- mu.l0^2 / sd.l0^2
alpha.l1 <- mu.l1^2 / sd.l1^2
beta.l0 <- mu.l0 / sd.l0^2
beta.l1 <- mu.l1 / sd.l1^2
rm(mu.l0, mu.l1, sd.l0, sd.l1)

# Specify data input for running BUGS/JAGS sampler
data <- list(
	y1 = y1, y2 = y2, n1 = n1, n2 = n2,
	shape1.p1 = shape1.p1, shape2.p1 = shape2.p1,
	mu.alpha = mu.alpha, tau.alpha = tau.alpha,
	mu.delta = mu.delta, tau.delta = tau.delta,
	alpha.l0 = alpha.l0, beta.l0 = beta.l0,
	alpha.l1 = alpha.l1, beta.l1 = beta.l1
)

# Specify BUGS/JAGS script
filein <- "rct_mod.txt"

# Specify parameters to sample
params <- c("p1", "p2", "rho", "alpha", "delta", "l0", "l1")

## Running Model ---------------------------------------------------------------

model.bugs <- R2OpenBUGS::bugs(
	data = data, inits = NULL, parameters.to.save = params, model.file = filein,
	n.chains = 2,
	n.iter = 10000,
	n.burnin = 100,
	n.thin = 1,
	DIC = TRUE
	# debug = TRUE
)
print(model.bugs, digits = 3)

# model.jags <- R2jags::jags(
# 	data = data, inits = NULL, parameters.to.save = params, model.file = filein,
# 	n.chains = 2,
# 	n.iter = 10000,
# 	n.burnin = 100,
# 	n.thin = 1,
# 	DIC = TRUE
# )
# print(model.jags)

## Trace plots -----------------------------------------------------------------

# par(mfrow = c(4, 2))
# for (i in params) {
# 	plot(model.bugs$sims.array[ , 1, i], type = "l", col = "red",
# 		 main = paste("Traceplot for", i), ylab = "value")
# 	lines(model.bugs$sims.array[ , 2, i], col = "blue")
# }

# CEA ==========================================================================

## Reading output --------------------------------------------------------------

# Read output from model.bugs
p1 <- model.bugs$sims.list$p1
p2 <- model.bugs$sims.list$p2
l0 <- model.bugs$sims.list$l0
l1 <- model.bugs$sims.list$l1
n.sims <- model.bugs$n.sims

# Or read output from model.jags
# p1 <- model.jags$BUGSoutput$sims.list$p1
# p2 <- model.jags$BUGSoutput$sims.list$p2
# l0 <- model.jags$BUGSoutput$sims.list$l0
# l1 <- model.jags$BUGSoutput$sims.list$l1
# n.sims <- model.jags$BUGSoutput$n.sims

## Analysis with BCEA ----------------------------------------------------------

# Population average cost & benefits
c <- e <- matrix(nrow = n.sims, ncol = 2)
c[ ,1] <- p1 * (168.19 * l1) + (1 - p1) * (113.61 * l0)
c[ ,2] <- p2 * (168.19 * l1) + (1 - p2) * (113.61 * l0) + 201.47
e[ ,1] <- (p1 * (l1 * 0.0013151 + (365 - l1)) + (1 - p1) * (l0 * 0.0025205 + (365 - l0)))
e[ ,2] <- (p2 * (l1 * 0.0013151 + (365 - l1)) + (1 - p2) * (l0 * 0.0025205 + (365 - l0)))

library(BCEA)
he <- bcea(e, c, ref = 2, c("status quo", "antibiotics prophylaxis"), Kmax = 1000)
ceplane.plot(he, wtp = 1000)
eib.plot(he, plot.cri = TRUE)
contour(he)
ceac.plot(he)
summary(he)
