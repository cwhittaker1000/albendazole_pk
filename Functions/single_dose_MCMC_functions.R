# Prior Function 
single_dose_prior_function <- function(parameters_vector, informative_prior) {
  
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
    k_abs_prior_value <- dtmvnorm(k_abs, mean = 5, sigma = 7, lower = 0.00001, log = TRUE)
    bioavailability_prior_value <- dtmvnorm(bioavailability, mean = 0.00001, sigma = 1, lower = 0.00001, log = TRUE)
    sigma_prior_value <- dtmvnorm(sigma, mean = 15, sigma = 7, lower = 0.5, log = TRUE)
    k_alb_so_prior_value <- dtmvnorm(k_alb_so, mean = 0, sigma = 1, lower = 0.00001, log = TRUE)
    k_alb_prior_value <- dtmvnorm(k_alb, mean = 0, sigma = 1, lower = 0.00001, log = TRUE)
  }

  prior_output <- k_abs_prior_value + bioavailability_prior_value + sigma_prior_value + k_alb_so_prior_value + k_alb_prior_value 
  return(unname(prior_output))
  
}

#Likelihood Function
single_dose_single_dose_loglikelihood_function <- function(parameters_vector, Alb_data, Alb_SO_data, metabolite_availability, model_instance, dose_info, time_increment) {
  
  # Extracting relevant dose information
  dose <- dose_info$amount 

  # Running the model and extracting relevant output
  model_runner <- Albendazole_PK_Model(k_abs = parameters_vector[, "k_abs"], sigma = parameters_vector[, "sigma"], 
                                       k_alb_so = parameters_vector[, "k_alb_so"], k_alb = parameters_vector[, "k_alb"], 
                                       gut_1 = (parameters_vector[, "bioavailability"]/1e+5) * (dose * 1e+6), 
                                       gut_2 = 0, liver = 0, blood_alb = 0, blood_alb_so = 0)

  times <- seq(0, max(Alb_SO_data$Time), time_increment)
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
        Alb_subset[i] <- 0.00000001
      }
    }
    Alb_SO_time_index <- match.closest(Alb_SO_data$Time, Times)
    Alb_SO_subset <- Alb_SO[Alb_SO_time_index]
    for (i in 1:length(Alb_SO_subset)) {
      if (Alb_SO_subset[i] < 0) {
        Alb_SO_subset[i] <- 0.00000001
      }
    }
  } else if (metabolite_availability == "Sulfoxide_Only") {
    Alb_SO_time_index <- match.closest(Alb_SO_data$Time, Times)
    Alb_SO_subset <- Alb_SO[Alb_SO_time_index]
    for (i in 1:length(Alb_SO_subset)) {
      if (Alb_SO_subset[i] < 0) {
        Alb_SO_subset[i] <- 0.00000001
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
  return(list(likelihood = overall_loglikelihood, Times = model_output[, 1], Alb_output = model_output[, 5], Alb_SO_output = model_output[, 6]))
}

#Likelihood Function
single_dose_multi_dose_loglikelihood_function <- function(parameters_vector, Alb_data, Alb_SO_data, metabolite_availability, model_instance, dose_info, time_increment) {
  
  # Extracting relevant dose information
  dose_amounts <- dose_info$amount
  dose_times <- dose_info$times
  
  for (i in 1:length(dose_times)) {
    if (i == 1) {
      mod <- Albendazole_PK_Model(k_abs = parameters_vector[, "k_abs"], sigma = parameters_vector[, "sigma"], 
                                  k_alb_so = parameters_vector[, "k_alb_so"], k_alb = parameters_vector[, "k_alb"], 
                                  gut_1 = (parameters_vector[, "bioavailability"]/1e+5) * (dose_amounts[i] * 1e+6), 
                                  gut_2 = 0, liver = 0, blood_alb = 0, blood_alb_so = 0)
    } else {
      mod <- Albendazole_PK_Model(k_abs = parameters_vector[, "k_abs"], sigma = parameters_vector[, "sigma"], 
                                  k_alb_so = parameters_vector[, "k_alb_so"], k_alb = parameters_vector[, "k_alb"], 
                                  gut_1 = final_vals["Gut_1"] + (parameters_vector[, "bioavailability"]/1e+5) * (dose_amounts[i] * 1e+6), 
                                  gut_2 = final_vals["Gut_2"], liver = final_vals["Liver"], blood_alb = final_vals["Blood_Alb"], 
                                  blood_alb_so = final_vals["Blood_Alb_SO"])
    }
    if (i != length(dose_times)) {
      times <- seq(0, dose_times[i + 1] - dose_times[i], time_increment)
    } else {
      times <- seq(0, max(Alb_SO_data$Time) - dose_times[i], time_increment)
    }
    
    model_output <- mod$run(times)
    if (i == 1) {
      final <- model_output
    } else {
      model_output[, "t"] <- model_output[, "t"] + max(final[, "t"]) 
      final <- rbind(final[-dim(final)[1], ], model_output) 
    }
    final_vals <- single_dose_extract_final_values(model_output)
  }
  
  Times <- final[, "t"]
  Alb <- final[, "Blood_Alb"]      
  Alb_SO <- final[, "Blood_Alb_SO"]
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
  return(list(likelihood = overall_loglikelihood, Times = final[, "t"], Alb_output = final[, "Blood_Alb"], Alb_SO_output = final[, "Blood_Alb_SO"]))
}

# Proposal Function
proposal_function <- function(current_parameters, covariance_matrix, scaling_factor) {
  number_parameters <- length(current_parameters)
  proposed_parameter_values <- rtmvnorm(1, mean = current_parameters, sigma = scaling_factor * covariance_matrix, lower = rep(0, length = number_parameters), 
                                        upper = rep(Inf, length = number_parameters), algorithm = c("rejection"))
  return(proposed_parameter_values)
}  

# Posterior Function
single_dose_posterior_function <- function(parameters_vector, Alb_data, Alb_SO_data, metabolite_availability, model_instance, dose_info, informative_prior, dosage, time_increment){
  if (dosage == "Single") {
    loglik <- single_dose_single_dose_loglikelihood_function(parameters_vector, Alb_data, Alb_SO_data, metabolite_availability, model_instance, dose_info, time_increment)
    posterior <- loglik$likelihood + single_dose_prior_function(parameters_vector, informative_prior)
    alb <- loglik$Alb_output
    alb_so <- loglik$Alb_SO_output
    times <- loglik$Times
  } else if (dosage == "Multiple") {
    loglik <- single_dose_multi_dose_loglikelihood_function(parameters_vector, Alb_data, Alb_SO_data, metabolite_availability, model_instance, dose_info, time_increment)
    posterior <- loglik$likelihood + single_dose_prior_function(parameters_vector, informative_prior)
    alb <- loglik$Alb_output
    alb_so <- loglik$Alb_SO_output
    times <- loglik$Times
  } else {
    stop("Error - specify single or multiple")
  }
  return(list(posterior = posterior, alb = alb, alb_so = alb_so, times = times))
}

# Extract Final Values
single_dose_extract_final_values <- function(model_output) {
  end <- dim(model_output)[1]
  Gut_1 <- model_output[end, "Gut_1"]
  Gut_2 <- model_output[end, "Gut_2"]
  Liver <- model_output[end, "Liver"]
  Blood_Alb <- model_output[end, "Blood_Alb"]
  Blood_Alb_SO <- model_output[end, "Blood_Alb_SO"]
  return(c(Gut_1, Gut_2, Liver, Blood_Alb, Blood_Alb_SO))
}

# Function up Updating Proposals Based on Johnstone-Chang Algorithm updating the scaling factor
jc_prop_update <- function(accepted, i, current_sf, previous_mu, current_parameters,
                           current_covariance_matrix, required_acceptance_ratio) {
  
  cooldown <- (i + 1) ^ -0.6
  new_covariance_matrix <- ((1 - cooldown) * current_covariance_matrix) +
    (cooldown * (t(current_parameters - previous_mu) %*% (current_parameters - previous_mu)))
  new_mu <- ((1 - cooldown) * previous_mu) + (cooldown * current_parameters)
  log_new_scaling_factor <- log(current_sf) + cooldown * (accepted - required_acceptance_ratio)
  new_scaling_factor = exp(log_new_scaling_factor);
  
  return(list("covariance_matrix" = new_covariance_matrix,
              "mu" = new_mu,
              "scaling_factor" = new_scaling_factor))
}

# MCMC Function
single_dose_MCMC_running <- function(number_of_iterations, parameters_vector, sd_proposals, start_covariance_adaptation, 
                                Alb_data, Alb_SO_data, metabolite_availability, model_instance, dose_info, dosage,
                                output, informative_prior, reparam, refresh, time_increment) {
  
  # Scaling Factor Adaptation Parameters 
  scaling_factor <- 1
  required_acceptance_ratio <- 0.1 # required acceptance probability
  number_of_parameters <- length(parameters_vector)
  alb_output_storage <- matrix(nrow = number_of_iterations, ncol = length(seq(0, max(Alb_SO_data$Time), time_increment)))
  alb_SO_output_storage <- matrix(nrow = number_of_iterations, ncol = length(seq(0, max(Alb_SO_data$Time), time_increment)))
  
  # Storage for Output
  prop_param <- matrix(nrow = number_of_iterations + 1, ncol = number_of_parameters)
  MCMC_output <- matrix(nrow = number_of_iterations + 1, ncol = number_of_parameters)
  MCMC_output[1, ] <- c(parameters_vector["k_abs"], parameters_vector["bioavailability"], parameters_vector["sigma"], 
                        parameters_vector["k_alb"], parameters_vector["k_alb_so"])
  colnames(MCMC_output) = c("k_abs", "bioavailability", "sigma", "k_alb", "k_alb_so")
  MCMC_output[1, ] <- MCMC_output[1, ] * reparam ## reparameterisation step
  Acceptances <- vector(length = number_of_iterations + 1)
  Acceptance_Ratio <- vector(length = number_of_iterations + 1)
  Logposterior_storage <- vector(length = number_of_iterations + 1)
  Correction_problems <- vector(mode = "numeric", length = number_iterations)
  
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
  reparam_current_parameter_values <- t(as.matrix(current_parameter_values / reparam)) ## reverse reparameterisation for model running 
  raw_current_posterior <- single_dose_posterior_function(reparam_current_parameter_values, Alb_data, Alb_SO_data, metabolite_availability, 
                                                          model_instance, dose_info, informative_prior, dosage, time_increment)
  current_posterior <- raw_current_posterior$posterior
  
  # Running the Actual MCMC
  for (i in 1:number_of_iterations){
    
    current_parameter_values <- MCMC_output[i, ]
    proposed_parameter_values <- proposal_function(MCMC_output[i, ], sigma, scaling_factor)
    reparam_proposed_parameter_values <-  proposed_parameter_values / reparam ## reverse reparameterisation for model running 
    colnames(reparam_proposed_parameter_values) <- c("k_abs", "bioavailability", "sigma", "k_alb", "k_alb_so")
    raw_proposed_posterior <- single_dose_posterior_function(reparam_proposed_parameter_values, Alb_data, Alb_SO_data, metabolite_availability, 
                                                             model_instance, dose_info, informative_prior, dosage, time_increment)
    proposed_posterior <- raw_proposed_posterior$posterior
    
    uncorrected_likelihood_ratio <- exp(proposed_posterior - current_posterior)
    if (sum(eigen(sigma)$values < 0) > 0) {
      sigma <- sigma + diag(0.000001, nrow = number_of_parameters)
    }
    prob_current_params_given_proposed <- dtmvnorm(current_parameter_values, mean = c(proposed_parameter_values), sigma = scaling_factor * sigma, lower = rep(0, length = number_of_parameters), upper = rep(Inf, length = number_of_parameters))
    prob_proposed_params_given_current <- dtmvnorm(proposed_parameter_values, mean = c(current_parameter_values), sigma = scaling_factor * sigma, lower = rep(0, length = number_of_parameters), upper = rep(Inf, length = number_of_parameters))
    correction_factor <- prob_current_params_given_proposed/prob_proposed_params_given_current
    corrected_likelihood_ratio <- uncorrected_likelihood_ratio * correction_factor
    if (is.na(correction_factor)) {
      correction_factor <- 1
      corrected_likelihood_ratio <- uncorrected_likelihood_ratio * correction_factor
      Correction_problems[i] <- 1
    }
    
    # print(c(i, proposed_posterior, current_posterior, uncorrected_likelihood_ratio, scaling_factor, 
    #         prob_current_params_given_proposed, prob_proposed_params_given_current,
    #         correction_factor, corrected_likelihood_ratio, proposed_parameter_values, current_parameter_values))
    
    if(runif(1) < corrected_likelihood_ratio) {
      MCMC_output[i + 1, ] <- proposed_parameter_values
      prop_param[i + 1, ] <- reparam_proposed_parameter_values
      Acceptances[i] <- 1
      Logposterior_storage[i] <- proposed_posterior
      current_posterior <- proposed_posterior
      raw_current_posterior <- raw_proposed_posterior
      alb_output_storage[i, ] <- raw_proposed_posterior$alb
      alb_SO_output_storage[i, ] <- raw_proposed_posterior$alb_so
    } 
    else{
      MCMC_output[i + 1, ] <- MCMC_output[i, ]
      prop_param[i + 1, ] <- prop_param[i, ]
      Acceptances[i] <- 0
      Logposterior_storage[i] <- current_posterior
      alb_output_storage[i, ] <- raw_current_posterior$alb
      alb_SO_output_storage[i, ] <- raw_current_posterior$alb_so
    }
    
    # Adaptation Step Involves the Johnstone Chang Algorithm
    if (i >= start_covariance_adaptation) {
      timing_cov <- i - start_covariance_adaptation + 1 # iteration relative to when covariance adaptation started
      if (i == start_covariance_adaptation) {
        previous_mu <- matrix(colMeans(MCMC_output[1:(start_covariance_adaptation+1), ]), nrow = 1) # double check whether the +1 should be there
        current_parameters <- matrix(MCMC_output[start_covariance_adaptation+1, ], nrow = 1)
        temp <- jc_prop_update(Acceptances[i], timing_cov, scaling_factor, previous_mu,
                               current_parameters, sigma, required_acceptance_ratio)
        scaling_factor <- temp$scaling_factor
        sigma <- temp$covariance_matrix
        previous_mu <- temp$mu
      } else {
        current_parameters <- matrix(MCMC_output[i+1, ], nrow = 1)
        temp <- jc_prop_update(Acceptances[i], timing_cov, scaling_factor, previous_mu,
                               current_parameters, sigma, required_acceptance_ratio)
        scaling_factor <- temp$scaling_factor
        sigma <- temp$covariance_matrix
        previous_mu <- temp$mu
      }
    }
    Acceptance_Ratio[i] <- sum(Acceptances)/i
    if(i %% refresh == 0 & output == TRUE) {
      print(c("The iteration number is", i))
      print(c("The acceptance ratio is", Acceptance_Ratio[i]))
      print(c("Mean is ", round(apply(MCMC_output, 2, mean, na.rm = TRUE), 4) / reparam))
    }
  }
  list <- list()
  list[["MCMC_Output"]] <- MCMC_output
  list[["Acceptances"]] <- Acceptances
  list[["Posterior_Output"]] <- Logposterior_storage
  list[["Alb"]] <- alb_output_storage
  list[["Alb_SO"]] <- alb_SO_output_storage
  list[["Times"]] <- raw_current_posterior$times
  list[["Corection_problems"]] <- Correction_problems
  return(list)
}
