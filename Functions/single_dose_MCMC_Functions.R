# Prior Function 
prior_function <- function(parameters_vector, informative_prior) {
  
  k_abs <- parameters_vector[1]
  bioavailability <- parameters_vector[2]
  sigma <- parameters_vector[3]
  k_alb_so <- parameters_vector[4]
  k_alb <- parameters_vector[5]
  
  if (informative_prior == TRUE) {
    k_abs_prior_value <- dtmvnorm(k_abs, mean = 5, sigma = 7, lower = 1, log = TRUE)
    bioavailability_prior_value <- dtmvnorm(bioavailability, mean = 0.005, sigma = 0.01, lower = 0, log = TRUE)
    sigma_prior_value <- dtmvnorm(sigma, mean = 15, sigma = 7, lower = 5, log = TRUE)
    k_alb_so_prior_value <- dtmvnorm(k_alb_so, mean = 0.125, sigma = 0.04, lower = 0, log = TRUE)
    k_alb_prior_value <- dtmvnorm(k_alb, mean = 0.125, sigma = 0.04, lower = 0, log = TRUE)
  } else {
    k_abs_prior_value <- dtmvnorm(k_abs, mean = 5, sigma = 7, lower = 1, log = TRUE)
    bioavailability_prior_value <- dtmvnorm(bioavailability, mean = 0.005, sigma = 0.01, lower = 0, log = TRUE)
    sigma_prior_value <- dtmvnorm(sigma, mean = 15, sigma = 7, lower = 0.5, log = TRUE)
    k_alb_so_prior_value <- dtmvnorm(k_alb_so, mean = 0.125, sigma = 0.10, lower = 0, log = TRUE)
    k_alb_prior_value <- dtmvnorm(k_alb, mean = 0.125, sigma = 0.10, lower = 0, log = TRUE)
  }

  prior_output <- k_abs_prior_value + bioavailability_prior_value + sigma_prior_value + k_alb_so_prior_value + k_alb_prior_value 
  return(unname(prior_output))
  
}

#Likelihood Function
loglikelihood_function <- function(parameters_vector, Alb_data, Alb_SO_data, metabolite_availability, model_instance, end_time, dose) {

  # Running the model and extracting relevant output
  model_runner <- Albendazole_PK_Model(k_abs = parameters_vector[1], bioavailability = parameters_vector[2], 
                                       sigma = parameters_vector[3], k_alb_so = parameters_vector[4],
                                       k_alb = parameters_vector[5], dose = dose)

  times <- seq(0, end_time, 0.01)
  model_output <- model_runner$run(times)
  Times <- model_output[, 1]
  Alb <- model_output[, 5]
  Alb_SO <- model_output[, 6]
  if (metabolite_availability == "Both") {
    Alb_time_index <- match.closest(Alb_data$Time, Times)
    Alb_subset <- Alb[Alb_time_index]
    
    if(sum(is.na(Alb_subset)) > 0) {
      print(Alb_subset)
      print(parameters_vector)
    }
    for (i in 1:length(Alb_subset)) {
      if (Alb_subset[i] < 0) {
        Alb_subset[i] <- 0
      }
    }
    Alb_SO_time_index <- match.closest(Alb_SO_data$Time, Times)
    Alb_SO_subset <- Alb_SO[Alb_SO_time_index]
    for (i in 1:length(Alb_SO_subset)) {
      if (Alb_SO_subset[i] < 0) {
        Alb_SO_subset[i] <- 0
      }
    }
  } else if (metabolite_availability == "Sulfoxide_Only") {
    Alb_SO_time_index <- match.closest(Alb_SO_data$Time, Times)
    Alb_SO_subset <- Alb_SO[Alb_SO_time_index]
    # print(Alb_SO_subset)
    for (i in 1:length(Alb_SO_subset)) {
      if (Alb_SO_subset[i] < 0) {
        Alb_SO_subset[i] <- 0
      }
    }
  }

  # Calculate the loglikelihoods 
  if (metabolite_availability == "Both") {
    loglik_Alb <- sum(dpois(round(Alb_data$Alb), lambda = Alb_subset, log = TRUE))
    loglik_Alb_SO <- sum(dpois(round(Alb_SO_data$Alb_SO), lambda = Alb_SO_subset, log = TRUE))
    overall_loglikelihood <- loglik_Alb + loglik_Alb_SO 
    if(sum(is.nan(dpois(round(Alb_data$Alb), lambda = Alb_subset, log = TRUE))) > 0 | sum(is.na(dpois(round(Alb_data$Alb), lambda = Alb_subset, log = TRUE))) > 0) {
      print(round(Alb_data$Alb))
      print(Alb_subset)
      print(parameters_vector)
      print(overall_loglikelihood)
    }
    if(sum(is.nan(dpois(round(Alb_SO_data$Alb_SO), lambda = Alb_SO_subset, log = TRUE))) > 0 | sum(is.nan(dpois(round(Alb_SO_data$Alb_SO), lambda = Alb_SO_subset, log = TRUE))) > 0) {
      print(round(Alb_SO_data$Alb_SO))
      print(Alb_SO_subset)
      print(parameters_vector)
      print(overall_loglikelihood)
    }
  } else if (metabolite_availability == "Sulfoxide_Only") {
    loglik_Alb_SO <- sum(dpois(round(Alb_SO_data$Alb_SO), lambda = Alb_SO_subset, log = TRUE))
    overall_loglikelihood <- loglik_Alb_SO 
    if(sum(is.nan(dpois(round(Alb_SO_data$Alb_SO), lambda = Alb_SO_subset, log = TRUE))) > 0) {
      print(round(Alb_SO_data$Alb_SO))
      print(Alb_SO_subset)
      print(parameters_vector)
      print(overall_loglikelihood)
    }
  }
  return(overall_loglikelihood)
}

# Proposal Function
proposal_function <- function(current_parameters, covariance_matrix, scaling_factor) {
  number_parameters <- length(current_parameters)
  proposed_parameter_values <- rtmvnorm(1, mean = current_parameters, sigma = scaling_factor * covariance_matrix, lower = rep(0, length = number_parameters), 
                                        upper = rep(Inf, length = number_parameters), algorithm = c("rejection"))
  return(proposed_parameter_values)
}  

# Posterior Function
posterior_function <- function(parameters_vector, Alb_data, Alb_SO_data, metabolite_availability, model_instance, end_time, dose, informative_prior){
  posterior <- loglikelihood_function(parameters_vector, Alb_data, Alb_SO_data, metabolite_availability, model_instance, end_time, dose) + prior_function(parameters_vector, informative_prior)
  return(posterior)
}

# MCMC Function
MCMC_running <- function(number_of_iterations, parameters_vector, sd_proposals, start_covariance_adaptation, Alb_data, Alb_SO_data, metabolite_availability, model_instance, end_time, dose, output, informative_prior, reparam) {
  
  # Scaling Factor Adaptation Parameters 
  number_of_parameters <- length(parameters_vector)
  required_acceptance_probability <- 0.06 # required acceptance probability
  alpha_r <- -qnorm(required_acceptance_probability/2) # same as -qtrunc(p=p_r.accept/2, spec="norm", a = -Inf, b = Inf)
  r <- array(NA, dim = c(number_of_iterations + 1, 1))
  c_r <- array(NA, dim = c(number_of_iterations, 1))
  r[1, ] <- 1 # rexp(n = 1, rate = 1) # generating a random value to begin with
  acceptance_tracker <- array(0, dim = c(number_of_iterations + 1, 1))
  covariance_matrix_storage <- list()
  
  # Storage for Output
  prop_param <- matrix(nrow = number_of_iterations + 1, ncol = number_of_parameters)
  MCMC_output <- matrix(nrow = number_of_iterations + 1, ncol = number_of_parameters)
  MCMC_output[1, ] <- c(parameters_vector["k_abs"], parameters_vector["bioavailability"], parameters_vector["sigma"], 
                        parameters_vector["k_alb_so"], parameters_vector["k_alb"])
  colnames(MCMC_output) = c("k_abs", "bioavailability", "sigma", "k_alb_so", "k_alb")
  MCMC_output[1, ] <- MCMC_output[1, ] * reparam ## reparameterisation step
  Acceptances <- vector(length = number_of_iterations + 1)
  Acceptance_Ratio <- vector(length = number_of_iterations + 1)
  Logposterior_storage <- vector(length = number_of_iterations + 1)
  
  # Generating the Covariance Matrix
  sigma <- matrix(0, nrow = number_of_parameters, ncol = number_of_parameters)
  for (i in 1:number_of_parameters) {
    for (j in 1:number_of_parameters) {
      if (i == j) {
        sigma[i, j] <- sd_proposals[i]
      }
    }
  }
  
  # Calculating the Posterior Probability for the Initial Values
  current_parameter_values <- MCMC_output[1, ]
  reparam_current_parameter_values <- current_parameter_values / reparam ## reverse reparameterisation for model running 
  current_posterior <- posterior_function(reparam_current_parameter_values, Alb_data, Alb_SO_data, metabolite_availability, model_instance, end_time, dose, informative_prior)
  
  # Running the Actual MCMC
  for (i in 1:number_of_iterations){
    
    current_parameter_values <- MCMC_output[i, ]
    proposed_parameter_values <- proposal_function(MCMC_output[i, ], sigma, r[i, ])
    reparam_proposed_parameter_values <-  proposed_parameter_values / reparam ## reverse reparameterisation for model running 
    proposed_posterior <- posterior_function(reparam_proposed_parameter_values, Alb_data, Alb_SO_data, metabolite_availability, model_instance, end_time, dose, informative_prior)
    
    uncorrected_likelihood_ratio <- exp(proposed_posterior - current_posterior)
    prob_current_params_given_proposed <- dtmvnorm(current_parameter_values, mean = c(proposed_parameter_values), sigma = r[i, ] * sigma, lower = rep(0, length = number_of_parameters), upper = rep(Inf, length = number_of_parameters))
    prob_proposed_params_given_current <- dtmvnorm(proposed_parameter_values, mean = c(current_parameter_values), sigma = r[i, ] * sigma, lower = rep(0, length = number_of_parameters), upper = rep(Inf, length = number_of_parameters))
    correction_factor <- prob_current_params_given_proposed/prob_proposed_params_given_current
    corrected_likelihood_ratio <- uncorrected_likelihood_ratio * correction_factor
    
    if(runif(1) < corrected_likelihood_ratio) {
      MCMC_output[i + 1, ] <- proposed_parameter_values
      prop_param[i + 1, ] <- reparam_proposed_parameter_values
      Acceptances[i] <- 1
      Logposterior_storage[i] <- proposed_posterior
      current_posterior <- proposed_posterior
      acceptance_tracker[i, ] <- 1
    } 
    else{
      MCMC_output[i + 1, ] <- MCMC_output[i, ]
      prop_param[i + 1, ] <- prop_param[i, ]
      Acceptances[i] <- 0
      Logposterior_storage[i] <- current_posterior
      acceptance_tracker[i, ] <- 0
    }

    # Adaptation Step Involves the Robins-Munro Adaptation of r, the SD of the Proposal Distribution
    c_r[i, ] <- (1 - (1/number_of_parameters)) * ((2 * pi ^ 0.5 * exp((alpha_r ^ 2)/2)) / (2 * alpha_r)) + 1/(number_of_parameters * required_acceptance_probability * (1 - required_acceptance_probability))
    r[i + 1, ] <- abs(ifelse((acceptance_tracker[i, ] == 0), r[i, ] - (c_r[i, ] * required_acceptance_probability)/max(1000, i/number_of_parameters), 
                                                             r[i, ] + (c_r[i, ] * (1 - required_acceptance_probability))/max(1000, i/number_of_parameters))) 
    # Covariance Matrix Adaptation
    if (i > start_covariance_adaptation) {
      if (i == start_covariance_adaptation + 1) {
        previous_mean <- as.matrix(apply(MCMC_output[1:start_covariance_adaptation, ], 2, mean, na.rm = TRUE), nrow = number_of_parameters, ncol = 1)
        sigma <- cov(MCMC_output[1:start_covariance_adaptation, ])
      }
      new_mean <- as.matrix((1/i) * ((i-1) * previous_mean + MCMC_output[i + 1, ]), nrow = number_of_parameters, ncol = 1)
      sigma <- (((i-2)/(i-1)) * sigma) + (previous_mean %*% t(previous_mean)) - ((i/(i-1)) * (new_mean %*% t(new_mean))) + ((1/(i-1)) * (as.matrix(MCMC_output[i + 1, ]) %*% t(as.matrix(MCMC_output[i + 1, ])))) 
      if (i %% 1000 == 0 & output == TRUE) {
        print(c("scaling factor is ", r[i, ]))
        print(r[i, ] * sigma)
        print(c("running mean is ", apply(MCMC_output, 2, mean, na.rm = TRUE)))
      }
      previous_mean <- new_mean
      #print(reparam_proposed_parameter_values)
    }
    covariance_matrix_storage[[i]] <- r[i, ]  * sigma
    Acceptance_Ratio[i] <- sum(Acceptances)/i
    if(i %% 1000 == 0 & output == TRUE) {
      print(c("The iteration number is", i))
      print(c("The acceptance ratio is", Acceptance_Ratio[i]))
    }
  }
  
  list <- list()
  list[["MCMC_Output"]] <- MCMC_output
  list[["Acceptances"]] <- Acceptances
  list[["Posterior_Output"]] <- Logposterior_storage
  list[["Covariance_Matrices"]] <- covariance_matrix_storage
  list[["rm_r"]] <- r
  list[["rm_c"]] <- c_r
  return(list)
  
}

