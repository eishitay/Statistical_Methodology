# Load the required libraries 
library(ggplot2)
library(coda)

# Data and prior
y = rnorm(1000)
m = 0
C = 10
v = 1
s = 1

# Initial values 
mu = 0
sigma2 = 1

# Save structures 
n_iter = 10000
mu_keep = rep(NA, n_iter)
sigma_keep = rep(NA, n_iter)

# Pre-calculations
n = length(y)
sum_y = sum(y)
vp = v + n 
vs2 = v*s^2

# Gibbs Sampler 
for (i in 1:n_iter) {
  # Sample mu 
  Cp = 1/(n/sigma2 + 1/C)
  mp = Cp*(sum_y/sigma2 + m/C)
  mu = rnorm(1,mp,sqrt(Cp))
  
  # Sample sigma 
  vpsp2 = vs2 + sum((y-mu)^2)
  sigma2 = 1/rgamma(1,vp/2,vpsp2/2)
  
  # Save the iterations 
  mu_keep[i] = mu 
  sigma_keep[i] = sqrt(sigma2)
}


# Traceplots for mu and sigma 
traceplot(as.mcmc(mu_keep), main=expression("Traceplot for"~ mu), xlab = "Iterations", ylab = "value")
traceplot(as.mcmc(sigma_keep), main=expression("Traceplot for"~ sigma^2), xlab = "Iterations", ylab = "value")

plot(density(mu_keep), main = "", xlab = expression(~ mu ~" values"))
plot(density(sigma_keep), main = "", xlab = expression(~ sigma^2 ~" values"))

acf(as.mcmc(mu_keep) , main = "", ylab = expression(~mu))
acf(as.mcmc(sigma_keep) , main = "", ylab = expression(~ sigma^2))
















