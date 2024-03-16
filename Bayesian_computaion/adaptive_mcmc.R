# Load the required libraries 
library(MASS)
library(mvtnorm)
library(ggplot2)
library(coda)

# The target distribution to sample from is a bi-variate normal distribution
mean_val <- c(10,10)
cov_matrix <- matrix(c(2,1,1,2),nrow = 2, ncol = 2)
target_distribution <- function(x,mean_val = mean_val,cov_matrix = cov_matrix)
{
  return(dmvnorm(x,mean = c(10,10), sigma = matrix(c(2,1,1,2),nrow = 2, ncol = 2)))
}

# Initial value of x0,m0,c0, sigma_d and beta values 
initial_x0 <- c(0,0)
initial_m0 <- c(0,0)
initial_c0 <- matrix(c(1,0,0,1),nrow = 2, ncol = 2)
sigma_d <- 0.1
beta <- 0.05 


# Adaptive MCMC algorithm 
adaptive_mcmc <- function(initial_x0,initial_m0,initial_c0,sigma_d,beta, num_samples)
{
  d <- 2 
  # Initialize a matrix to collect the samples 
  samples <- matrix(0, nrow = num_samples, ncol = d)
  sigma_opt <- sigma_d
  x_n <- initial_x0
  c_n <- initial_c0
  m_n <- initial_m0
  # Calculate the total samples that were accepted 
  acceptance_no <- 0 
  for (i in 1:num_samples) 
  {
    # Propose new values 
    if (num_samples <= 2*d)
    {
      x_prop <- mvrnorm(mu = x_n,Sigma = initial_c0)
    }
    else
    {
      sigma_opt <- (2.38**2)/d
      sigma_init <- (0.1**2)/d
      I_d <- matrix(c(1,0,0,1),nrow = 2, ncol = 2)
      x_prop <- (1-beta)*mvrnorm(1,mu = x_n, Sigma = sigma_opt*c_n)+beta*mvrnorm(1,mu = x_n, Sigma = sigma_init*I_d)
    }
    print(x_prop)
    # Calculate acceptance ratio 
    alpha <- min(1,target_distribution(x_prop,mean_val = m_n,cov_matrix = c_n)/target_distribution(x_n,mean_val = m_n,cov_matrix = c_n))
    
    if(runif(1) < alpha)
    {
      x_n <- x_prop
      acceptance_no <- acceptance_no + 1
      # Update the mean and covraiance matrix 
      
    }
    else
    {
      x_n <- x_n
    }
    
    # Update the values of mean and covraiance parameters 
    diff <- x_n - m_n
    c_n <- c_n + (1/(num_samples+1))*((outer(diff,diff)) - c_n)
    m_n <- (num_samples/(num_samples+1))*m_n + (1/(num_samples+1))*(x_n - m_n)
    samples[i, ] <- x_n
    
  }
  # Calculate the acceptance ratio 
  acceptance_rate <- (acceptance_no / num_samples)*100
  cat('Acceptance Rate: ', acceptance_rate, '\n')
  # Return a list of sampled values and acceptance rates 
  return(list(Samples = samples, AR = acceptance_rate))
}

num_samples <- 100000
# Collect samples from the target density using adaptive mcmc algorithm 
samples <- adaptive_mcmc(initial_x0, initial_m0, initial_c0, sigma_d, beta =0.05, num_samples=10000)
samples

# Plot the Traceplots for sampled values to check for convergence 
plot_traceplot <- function(samples, true_value) {
  num_iterations <- length(samples)
  num_dimensions <- 1
  
  iterations <- rep(1:num_iterations, each = num_dimensions)
  dimension <- rep(1:num_dimensions, times = num_iterations)
  
  df <- data.frame(Iteration = iterations, Dimension = dimension, Sample = as.vector(samples))
  m_val <- mean(df$Sample)
  
  # Create data frame for legend
  legend_df <- data.frame(yintercept = c(true_value, m_val),
                          label = c("True Value", "Mean Value"),
                          linetype = c("dashed", "solid"),
                          color = c("red", "green"))
  
  ggplot(df, aes(x = Iteration, y = Sample, color = factor(Dimension))) + 
    geom_line() +
    geom_hline(data = legend_df, aes(yintercept = yintercept, linetype = linetype, color = color)) + # Thickening the green line
    geom_text(data = legend_df, aes(x = Inf, y = yintercept, label = label), color = "black", hjust = -0.2, vjust = -0.5) +
    scale_color_manual("", 
                       values = c("grey", "red", "green"),
                       labels = c("Sampled Value", "True Value", "Mean Value")) +
    scale_linetype_manual("", 
                          values = c("dashed", "solid"),
                          labels = c("True Value", "Mean Value")) +
    labs(title = "Traceplots",
         x = "Iteration",
         y = "Sample") + 
    theme_minimal() +
    theme(legend.position = "bottom")
}

plot_traceplot(samples$Samples[,1], 10)  
plot_traceplot(samples$Samples[,2],10)

# Function to plot contour plot of the target density with sampled points
plot_contourplot <- function(target_distribution, samples) {
  # Generate grid of x, y values
  x <- seq(min(samples[, 1]) - 1, max(samples[, 1]) + 1, length.out = 100)
  y <- seq(min(samples[, 2]) - 1, max(samples[, 2]) + 1, length.out = 100)
  grid <- expand.grid(x = x, y = y)
  
  # Evaluate target distribution on grid
  grid$z <- apply(grid, 1, function(xy) target_distribution(c(xy[1], xy[2])))
  
  # Plot contour plot
  contour_plot <- ggplot(grid, aes(x = x, y = y)) +
    geom_contour(aes(z = z, color = "Target Distribution"), bins = 20) +
    geom_point(data = as.data.frame(samples), aes(x = samples[, 1], y = samples[, 2], color = "Sampled Points")) + # Add points for samples
    scale_color_manual(name = "Legend", 
                       values = c("Target Distribution" = "black", "Sampled Points" = "lightgreen"), 
                       labels = c("Sampled Points", "Target Contours")) + # Add legend for both contour and sampled points
    theme_minimal() +
    labs(title = "Contour Plot with Sampled Points",
         x = "X",
         y = "Y") +
    guides(color = guide_legend(reverse = TRUE)) # Reverse the order of the legends
  
  return(contour_plot)
}

# Overlay contour plot of the target distribution with the sampled data points 
contour_plot <- plot_contourplot(target_distribution, samples$Samples)
print(contour_plot)

# Function to plot ACF plots for samples
plot_acf <- function(samples) {
  par(mfrow = c(1, 1))  # Set the layout to 1x1
  acf_result <- acf(samples, main = "Autocorrelation Function (ACF) Plot", ylab = "Autocorrelation", col = "blue")  # Plot ACF
  return(acf_result)
}

# Check for auto-correlation in the samples 
acf_plot <- plot_acf(samples$Samples[, 1])  # Pass the samples to the function
acf_plot <- plot_acf(samples$Samples[,2])

# Thin the mcmc samples chain 

thinned_samples <- function(thinning_factor,samples)
{
  thinned_samples <- samples[seq(1,length(samples),by = thinning_factor)]
  return(thinned_samples)
}

# Thin samples collected by the adaptive MCMC chain by a factor of 20 
new_samples <- thinned_samples(20,samples$Samples[,1])

# Plot the ACF plot of the newly thinned samples 
plot_acf(new_samples)


# Calculate the effective sample sizes 
ess <- effectiveSize(samples$Samples[,1])

# Run 5 chains from different start points and calculate the psrf value to check for convergence 

run_chains <- function(num_chains)
{
  chains_total <- list()
  for(i in seq(1, num_chains))
  {
    chain_num <- paste("chain", i)  # Creating the chain number string
    # Define the range for random initialization
    lower_bound <- 1
    upper_bound <- 20
    
    # Randomly instantiate initial_x0 values
    initial_x0 <- runif(length(initial_x0), min = lower_bound, max = upper_bound)
    
    # Print the randomly instantiated values
    cat("Initial X0:", initial_x0, "\n")
    
    c <- adaptive_mcmc(initial_x0, initial_m0, initial_c0, sigma_d, beta = 0.05, num_samples = 50000)
    chains_total[[i]] <- c$Samples  # Append chain number and samples
  }
  
  return(chains_total)
}


# Calculate the PSRF values 
calculate_metrics <- function(total_chains_list, d)
{
  chain_variances <- c()
  all_samples <- c()
  chain_means <- c()
  n <- length(total_chains_list[[1]][,1])
  m <- length(total_chains_list)
  # Calculate within chain variance 
  for (i in seq(1,m))
  {
    chain_variances <- c(chain_variances, sd(total_chains_list[[i]][,d])**2)
    all_samples <- c(all_samples,total_chains_list[[i]][,d])
    chain_means <- c(chain_means, mean(total_chains_list[[i]][,d]))
    
  }
  within_chain_variance <- sum(chain_variances)/length(chain_variances)
  total_mean <- mean(all_samples)
  pairwise_mean_diff_squared <- (total_mean - chain_means)**2
  between_chain_variance <- (n/(m-1))*sum(pairwise_mean_diff_squared)
  pooled_variance <- ((n-1)/n)*within_chain_variance + (1/m)*between_chain_variance
  psrf <- sqrt(pooled_variance/within_chain_variance)
  
  return(psrf)
  
}

# Run 3 chains and sample from them 
chains <- run_chains(3)
# Calculate the PSRF value for 3 chains 
calculate_metrics(chains,1)



