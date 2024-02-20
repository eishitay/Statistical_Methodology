# Loading the required libraries 
library(R2jags)
library(runjags)
library(tidyverse)
library(tidybayes)
library(fields)
library(mvtnorm)
library(invgamma)
library(MASS)


set.seed(1)
year <- sample(1:500,100) %>% sort
x <- year/100

# Create a matrix of distances between every combination
d <- fields::rdist(x)
phi <- 5
K <- exp(-(phi^2)*(d^2))
sigma_g <- 2
f <- rmvnorm(1, sigma=(sigma_g^2)*K)

#sigma_y <- matrix(c(0.3,0,0,0.3), ncol=2)

tau <- 0.7
eps <- rnorm(100,0,(tau)^2) # random noise
y <- c(f) + eps 

dat = tibble(year = year,
             x = x,
             y = y, 
             gp = c(y))

ggplot(dat, aes(x = year, y = y)) +
  geom_point(aes(colour = "Simulated Data")) +
  geom_line(aes(x = year, y = gp, colour = "Gaussian Process")) +
  labs(colour = "")

# The likelihood value 
likelihood_norm <- function(tau_value,sigma_value,phi_value,y)
{
  K <- sigma_value^2*exp(-(phi_value^2)*(d^2))
  H <- K + (tau_value^2)*diag(nrow(K))
  cat('H = ', dim(H), '\n')
  dmvnorm(y,sigma = H,log = TRUE)
}
#likelihood_norm(sqrt(0.3),2,5,y)

prior_tau_value <- function(tau2_value)
{
  
  dnorm(tau2_value^2,mean = 0.5,sd = 1, log = TRUE)
  
}

prior_sigma_value <- function(sigma2_value)
{
  dnorm(sigma2_value^2,mean = 4,sd = 1, log = TRUE)
}

prior_phi_value <- function(phi_value)
{
  dnorm(phi_value^2,mean = 25,sd = 1, log = TRUE)
}


metropolis_hasting <- function(iterations)
{
  # Initializations
  tau_initial <- 1
  tau_val <- tau_initial
  sigma_initial <- 1
  sigma_val <- sigma_initial
  phi_initial <- 1
  phi_val <- phi_initial
  
  phi_samples <- numeric(iterations)
  tau_samples <- numeric(iterations)
  sigma_samples <- numeric(iterations)
  
  for (i in 1:iterations)
  {
    #phi_proposed <- runif(1,min = 0.01, max = 30)
    phi_proposed <- rnorm(1,mean = 5, sd = 1)
    cat('phi prop:', phi_proposed,'\n')
    sigma_proposed <- rnorm(1,mean = 2, sd = 1)
    #sigma_proposed <- runif(1,min = 0.01, max = 10)
    cat('sigma prop:', sigma_proposed, '\n')
    #tau_proposed <- runif(1,min = 0.01, max = 1)
    tau_proposed <- rnorm(1,mean = 0.7, sd = 1)
    cat('tau proposed:', tau_proposed,'\n')
    #density_proposed <- likelihood_norm(tau_proposed,sigma_proposed,phi_proposed,y)*prior_phi_value(phi_proposed)*prior_sigma_value(sigma_proposed)*prior_tau_value(tau_proposed)
    density_current <-  likelihood_norm(tau_val,sigma_val,phi_val,y)+prior_phi_value(phi_val)+prior_sigma_value(sigma_val)+prior_tau_value(tau_val)
    density_proposed_tau <- likelihood_norm(tau_proposed,sigma_val,phi_val,y)+prior_phi_value(phi_val)+prior_sigma_value(sigma_val)+prior_tau_value(tau_proposed)
    density_proposed_phi <- likelihood_norm(tau_val,sigma_val,phi_proposed,y)+prior_phi_value(phi_proposed)+prior_sigma_value(sigma_val)+prior_tau_value(tau_val)
    density_proposed_sigma <- likelihood_norm(tau_val,sigma_proposed,phi_val,y)+prior_phi_value(phi_val)+prior_sigma_value(sigma_proposed)+prior_tau_value(tau_val)
    
    # Calculate alphas for different parameters 
    alpha_tau <- min(1,exp(density_proposed_tau/density_current))
    alpha_phi <- min(1,exp(density_proposed_phi/density_current))
    alpha_sigma <- min(1,exp(density_proposed_sigma/density_current))
    cat('tau density:', density_proposed_tau,'\n')
    cat('Alpha tau', alpha_tau,'\n')
    # Accept or reject samples for tau
    if (runif(1) < alpha_tau) {
      tau_val <- tau_proposed
      tau_samples[i] <- tau_val
    }
    else
    {
      tau_samples[i] <- tau_val
    }
    
    # Accept or reject samples for phi 
    if (runif(1) < alpha_phi) {
      phi_val <- phi_proposed
      phi_samples[i] <- phi_val
    }
    else
    {
      phi_samples[i] <- phi_val
    }
    
    # Accept or reject samples for sigma 
    if (runif(1) < alpha_sigma) {
      sigma_val <- sigma_proposed
      sigma_samples[i] <- sigma_val
    }
    else
    {
      sigma_samples[i] <- sigma_val
    }
    
    
    
  }
  sample_list <- list("tau" = tau_samples, "phi" = phi_samples, "sigma" = sigma_samples)
  return(sample_list)
  
}

s = metropolis_hasting(10000)
sigma2_simulation <- (mean(s$sigma))
tau2_simulation <- (mean(s$tau))
phi2_simulation <- (mean(s$phi))


cov_matrix <-  sigma2_simulation^2*exp(-(phi2_simulation^2)*(d^2))
new_noise <- tau2_simulation^2*diag(nrow(K))
fnew <- cov_matrix%*%solve(cov_matrix + new_noise)%*%y


f_mu_new = fnew
f_s_new = (tau2_simulation*diag(nrow(cov_matrix)))%*%solve(tau2_simulation*diag(nrow(cov_matrix))+cov_matrix)%*%cov_matrix
f_sims = mvrnorm(10000,mu = f_mu_new, Sigma = f_s_new)
f_simulated = colMeans(f_sims)

plot(x,y)
lines(x,f,col="red", lwd = 3)
lines(x,f_simulated,col = "blue",lwd = 2)
legend(4,3, legend = c("true","simulated"),fill = c("red","blue"))
plot(f,f_simulated, main = "", xlab = "True F", ylab = "Simulated F")

tau2_samples_lst = s$tau[1001:length(s$tau)]
sigma2_samples_lst = s$sigma[1001:length(s$sigma)]
phi2_samples_lst = s$phi[1001:length(s$phi)]

traceplot(mcmc(tau2_samples_lst), main = "", ylab = expression(~tau ~"values"))
traceplot(mcmc(sigma2_samples_lst), main = "", ylab = expression(~sigma ~ "values"))
traceplot(mcmc(phi2_samples_lst), main = "", ylab = expression(~phi ~ "values"))


# effective sample sizes 
effectiveSize(mcmc(tau2_samples_lst))
effectiveSize(mcmc(sigma2_samples_lst))
effectiveSize(mcmc(phi2_samples_lst))

# ACF plots 
acf(tau2_samples_lst,main = "", ylab = expression(~tau))
acf(sigma2_samples_lst,main = "", ylab = expression(~sigma ) )
acf(phi2_samples_lst,main = "", ylab = expression(~phi))

