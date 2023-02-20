# Playground ===================================================================






# Analysis Tools ===============================================================



# dataJags <- list("N", "y", "X", "h", "k")
# filein <- "modeNormal.txt" # The file name for the above BUGS/JAGS code
# params <- c("alpha", "beta", "sigma") # Parameters to monitor/return from sampler
# inits <- function() { # Auto generate starting point for MCMC
# 	list(
# 		alpha = rnorm(1),
# 		beta = rnorm(1),
# 		lsigma = rnorm(1)
# 	)
# }
# # If used blocking, can consider the following inits:
# inits <- function() {
# 	list(
# 		lsigma = rnorm(1),
# 		coef = rnorm(2, 0, 1) # Notece for simplicity, initialise coef using two independent draws from N(0,1)
# 	)
# }

# model <- jags(
# 	data = data, inits = inits, parameters.to.save = params, model.file = filein,
# 	n.chains = 2,
# 	n.iter = 10000,
# 	n.burin = 4500,
# 	n.thin = 1,
# 	DIC = TRUE
)

# print(model, digits = 3, intervals = c(0.025, 0.975)) # Interval changes the output quantiles, this may require load bmhe package
# names(model) # Inspect elements contained in the output
# attach.bugs(model$BUGSoutput) # Access from R the actual MCMC output
# 
# model$BUGSoutput$sims.array # A 3-dimensional array [iterations stored, number of chains, number of parameters monitored] that can be used to create a trace plot
# chain1 <- model$BUGSoutput$sims.array[1:500, 1, "theta"]
# chain2 <- model$BUGSoutput$sims.array[1:500, 2, "theta"]
# chains <- model$BUGSoutput$sims.array[1:500, 1:2, "theta"]
# plot(chain1, type = "l", col = "blue", xlab = "Iteration", ylab = "", ylim = range(chains))
# points(chain2, type = "l", col = "red")