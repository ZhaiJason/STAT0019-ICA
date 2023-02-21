# Data Input ===================================================================

rct_data <- read.csv("rct_data.csv")
head(rct_data)

rct.1 <- rct_data$x[rct_data$T == 0]
rct.2 <- rct_data$x[rct_data$T == 1]

# Sample data
y1 <- sum(rct.1)
y2 <- sum(rct.2)
n1 <- length(rct.1)
n2 <- length(rct.2)

rm(rct_data, rct.1, rct.2)

# Sampling Set-up ===============================================================

# Population information
mu.p1 <- 0.08 # (0.04 + 0.12) / 2
sd.p1 <- 0.015 # 0.08 / 3
shape1.p1 <- round(mu.p1 * ((1 - mu.p1) * mu.p1 / sd.p1^2 - 1)) # alpha
shape2.p1 <- round((1 - mu.p1) * ((1 - mu.p1) * mu.p1 / sd.p1^2 - 1)) # beta
rm(mu.p1, sd.p1)

x <- seq(0,0.2,0.001)
plot(x, dbeta(x, shape1.p1, shape2.p1), type = "l")
# plot(x, dbeta(x, 29, 320), type = "l") # Eric's prior
abline(v = c(0.04, 0.08, 0.12))
pbeta(0.12, shape1.p1, shape2.p1) - pbeta(0.04, shape1.p1, shape2.p1)

# Unstructured effects
mu.alpha <- 0
sd.alpha <- sqrt(10)
tau.alpha <- 1 / sd.alpha^2
mu.delta <- 0
sd.delta <- sqrt(10)
tau.delta <- 1 / sd.delta^2
rm(sd.alpha, sd.delta)

# Length of staying (1 represents complication case)
mu.l0 <- 4.36
mu.l1 <- 33.66
sd.l0 <- 1.8
sd.l1 <- 2.1
tau.l0 <- 1 / sd.l0
tau.l1 <- 1 / sd.l1
rm(sd.l0, sd.l1)

data <- list(
	y1 = y1, y2 = y2, n1 = n1, n2 = n2,
	shape1.p1 = shape1.p1, shape2.p1 = shape2.p1,
	mu.alpha = mu.alpha, tau.alpha = tau.alpha,
	mu.delta = mu.delta, tau.delta = tau.delta,
	mu.l0 = mu.l0, mu.l1 = mu.l1, tau.l0 = tau.l0, tau.l1 = tau.l1
)

filein <- "rct_mod.txt"

params <- c("p1", "p2", "rho", "alpha", "delta", "l0", "l1")

# Gibbs Sampling ===============================================================

model.bugs <- R2OpenBUGS::bugs(
	data = data, inits = NULL, parameters.to.save = params, model.file = filein,
	n.chains = 2,
	n.iter = 10000,
	n.burnin = 100,
	n.thin = 1,
	DIC = TRUE,
	# debug = TRUE
)

model.jags <- R2jags::jags(
	data = data, inits = NULL, parameters.to.save = params, model.file = filein,
	n.chains = 2,
	n.iter = 10000,
	n.burnin = 100,
	n.thin = 1,
	DIC = TRUE
)

print(model.bugs)
print(model.jags)





# print(model, digits = 3, intervals = c(0.025, 0.975)) # Interval changes the output quantiles, this may require load bmhe package