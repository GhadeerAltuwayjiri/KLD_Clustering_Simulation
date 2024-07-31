library(MASS)
library(flowClust)

# Define the KLD function
KLD <- function(fit){
  KL_total <- 0
  
  for(k in 1:9) {
    mu_p <- Memory.mgd7.PS80.15.2_9@mu[k,]
    mu_q <- fit@mu[k,]
    
    sigma_p <- Memory.mgd7.PS80.15.2_9@sigma[k, , ]
    sigma_q <- fit@sigma[k, , ]
    
    d <- length(mu_p)
    
    KL1 <- (log(det(sigma_q)) - log(det(sigma_p)) - d) 
    KL2 <- t(mu_p - mu_q) %*% solve(sigma_q, mu_p - mu_q)
    KL3 <- sum(diag(solve(sigma_q, sigma_p)))
    
    KL <- 0.5 * (KL1 + KL2 + KL3)
    KL_total <- KL_total + KL
  }
  KL_total
}

# Define the cov_sim variations
cov_sim_values <- c(0, 0.1, 0.2, 0.5, 0.7)

# Initialize a list to hold the simulated data for each cov_sim value
simulated_data_lists <- list()

# Iterate over the cov_sim values
for (cov_index in seq_along(cov_sim_values)) {
  cov_sim <- cov_sim_values[cov_index]
  
  # Initialize lists for this cov_sim
  simulated_sample_list <- list()
  flow_sim_list <- list()
  MM_list <- list()
  Zhat_list1 <- list()
  Ztoy_list1 <- list()
  pi1 <- list()
  
  # Adjust the number of simulations if necessary
  num_simulations <- 5
  
  # Loop to generate data, perform flowClust, and store in lists
  for (j in 1:num_simulations) {
    set.seed(123 + j)
    
    # Calculate standard deviation for channels
    sd_channel <- apply(mu, 1, sd)
    
    # Generate modified mu for each sample
    mu_samp_list <- list()
    for(k in 1:5) {
      mu_samp_list[[k]] <- mu + matrix(rnorm(6 * 9, sd = cov_sim * sd_channel), nrow = 6, ncol = 9)
    }
    
    # Simulate datasets
    Ztoy <- sample(1:G, size = n, prob = pi, replace = TRUE)
    ytoy <- t(sapply(Ztoy, function(z) MASS::mvrnorm(n = 1, mu = mu_samp_list[[j]][, z], Sigma = sigma[, , z])))
    colnames(ytoy) <- Memory.mgd7.PS80.15.2[[9]]@varNames
    
    simulated_sample <- new("flowFrame", exprs = as.matrix(ytoy), parameters = parameters(Memory.mgd7.PS80.15.S1$lymphocyte))
    flow_sim <- flowClust(simulated_sample, K = 9, trans = 0, nu = Inf)
    
    if (!is.null(flow_sim@z)) {
      simulated_sample_list[[j]] <- simulated_sample
      flow_sim_list[[j]] <- flow_sim
      pi1[[j]] <- flow_sim@w
      # Store the simulated data and clustering results
      MM_list[[j]] <- ytoy
      Zhat <- apply(flow_sim@z, 1, which.max)
      Zhat_list1[[j]] <- Zhat
      Ztoy_list1[[j]] <- Ztoy
    } else {
      warning(paste("Clustering failed for simulation", j))
    }
  }
  
  # Store the lists for this cov_sim value
  simulated_data_lists[[paste("cov_sim", cov_sim, sep = "_")]] <- list(
    simulated_samples = simulated_sample_list, 
    flow_sims = flow_sim_list, 
    MM = MM_list, 
    Zhat = Zhat_list1,
    pi = pi1,
    Ztoy = Ztoy_list1
  )
}

# Initialize lists to store fit models, KLD values, and Z values for each kappa and cov_sim combination
fit_list2 <- list()
KLD_mat2 <- list()
Z_mat2 <- list()
Z_mat5 <- list()
pi2 <- list()

# Loop over each cov_sim value in simulated_data_lists
for (cov_sim_name in names(simulated_data_lists)) {
  # Extract the list of simulated samples and their corresponding flow sims for the current cov_sim
  simulated_samples_for_cov_sim <- simulated_data_lists[[cov_sim_name]]$simulated_samples
  flow_sims_for_cov_sim <- simulated_data_lists[[cov_sim_name]]$flow_sims
  
  # Determine the number of simulations for the current cov_sim
  num_simulations_for_cov_sim <- length(simulated_samples_for_cov_sim)
  
  # Prepare sub-lists for storing results corresponding to the current cov_sim
  fit_list2[[cov_sim_name]] <- list()
  KLD_mat2[[cov_sim_name]] <- matrix(0, nrow = num_simulations_for_cov_sim, ncol = length(Kappa_vals))
  Z_mat2[[cov_sim_name]] <- list()
  Z_mat5[[cov_sim_name]] <- list()
  pi2[[cov_sim_name]] <- list()  # Initialize the list for storing mixture weights for each kappa
  
  # Prepare the reference fit for the current cov_sim
  reference_fit <- flow_sims_for_cov_sim[[1]]  # Assuming the first simulation can serve as a reference
  toy_ref <- reference_fit
  toy_ref@sigma <- aperm(sigma, c(3, 2, 1))  # Adjust sigma dimensions if necessary
  toy_ref@mu <- t(mu)  # Transpose mu if necessary
  toy_ref@w <- pi  # Use the mixture weights from the first simulation as an example
  reference_fit <- toy_ref
  
  # Loop over kappa values
  for (i in 1:length(Kappa_vals)) {
    # Adjust the prior for each kappa value using the reference fit
    prior.sd1 <- Ghadeer2Prior(reference_fit, kappa = Kappa_vals[i])
    
    # Prepare sub-lists for storing results for the current kappa within the current cov_sim
    fit_list2[[cov_sim_name]][[i]] <- list()
    Z_mat2[[cov_sim_name]][[i]] <- list()
    Z_mat5[[cov_sim_name]][[i]] <- list()
    pi2[[cov_sim_name]][[i]] <- list()  # Prepare a sub-list for mixture weights for the current kappa
    
    for (j in 1:num_simulations_for_cov_sim) {
      # Print the current status
      print(paste("Processing cov_sim =", cov_sim_name, ", kappa =", Kappa_vals[i], "and simulation =", j))
      
      # Fit the model to the simulated sample
      secondary_fit <- flowClust(simulated_samples_for_cov_sim[[j]], K = 9, prior = prior.sd1, lambda = 1, trans = 0, usePrior = "yes")
      
      # Store the fit object and other metrics for each kappa, cov_sim, and sample combination
      fit_list2[[cov_sim_name]][[i]][[j]] <- secondary_fit
      KLD_mat2[[cov_sim_name]][j, i] <- KLD(secondary_fit)
      Z_mat2[[cov_sim_name]][[i]][[j]] <- table(apply(secondary_fit@z, 1, which.max))
      Z_mat5[[cov_sim_name]][[i]][[j]] <- apply(secondary_fit@z, 1, which.max)
      pi2[[cov_sim_name]][[i]][[j]] <- secondary_fit@w  # Store mixture weights
    }
  }
}
