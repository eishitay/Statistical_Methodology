# Import the required libraries 
library(ggplot2)


# Exponential Gamma conjugate model 
set.seed(121)
n <- 100
lambda <- 1
yis <- rexp(n,rate = lambda)
a <- 1
b <- 1
a1 <- a + n 
b1 <- b + n*mean(yis)
posterior_mean <- a1/b1

# Simulate samples from gamma distribution
sim1000 <- rgamma(n = 1000, shape = a1, rate = b1)

# Plot the simulated data with an overlayed curve for n = 1000 samples  
hist(sim1000,freq = F, main = "N = 1000", xlim = c(0,2), xlab = expression(lambda))
curve(dgamma(x,a1,b1), from = 0, to = 2, add = T, col = "red", lwd = 2)
legend(1.5,3,legend = c("True posterior"), fill = "red")

# Plot the simulated data with an overlayed curve for n = 10 samples 
sim10 <- rgamma(n = 10, shape = a1, rate = b1)

# Plot the simulated data with an overlayed curve for n = 10 samples 
hist(sim10, freq = F, main = "N = 10", xlim =c(0,2), xlab = expression(lambda))
curve(dgamma(x,a1,b1), from = 0, to = 2, add = T, col = "red", lwd = 2)
legend(1.5,4,legend = c("True posterior"), fill = "red")


