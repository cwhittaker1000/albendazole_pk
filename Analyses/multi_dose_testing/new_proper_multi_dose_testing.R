# Load required libraries
library(odin); library(tidyverse); library(MALDIquant); library(tmvtnorm)

# Generate model instance and source relevant functiosn
source("multi_dose_testing/new_test_multi_dose_model_ODIN.R")
source("multi_dose_testing/new_multi_dose_MCMC_Functions.R")

# Generate fake data 
dose_times <- c(0, 10, 20, 30)
dose_amounts <- c(400000, 400000, 400000, 400000)
parameters_vector <- c(0.5, 0.001, 15, 0.2, 0.2, 0.05, 3, 5)

for (i in 1:length(dose_times)) {
  if (i == 1) {
    mod <- Albendazole_PK_Model(k_abs = parameters_vector[1], sigma = parameters_vector[3], 
                                k_alb = parameters_vector[4], k_alb_so_init = parameters_vector[5], 
                                k_env = parameters_vector[6], EC50 = parameters_vector[7], Emax = parameters_vector[8],
                                gut_1 = (parameters_vector[2]/1e+5) * (dose_amounts[i] * 1e+6), # remember equivalent change in likelihood function required
                                gut_2 = 0, liver = 0, blood_alb = 0, blood_alb_so = 0, enzyme = 0)
  } else {
    mod <- Albendazole_PK_Model(k_abs = parameters_vector[1], sigma = parameters_vector[3], 
                                k_alb = parameters_vector[4], k_alb_so_init = parameters_vector[5], 
                                k_env = parameters_vector[6], EC50 = parameters_vector[7], Emax = parameters_vector[8],
                                gut_1 = final_vals["Gut_1"] + (parameters_vector[2]/1e+5) * (dose_amounts[i] * 1e+6), 
                                gut_2 = final_vals["Gut_2"], liver = final_vals["Liver"], blood_alb = final_vals["Blood_Alb"], 
                                blood_alb_so = final_vals["Blood_Alb_SO"], enzyme = final_vals["Enzyme_conc"])
  }
  if (i != length(dose_times)) {
    times <- seq(0, dose_times[i + 1] - dose_times[i], 1)
  } else {
    times <- seq(0, dose_times[i] - dose_times[i-1] + 15, 1)
  }
  model_output <- mod$run(times)
  if (i == 1) {
    final <- model_output
  } else {
    model_output[, "t"] <- model_output[, "t"] + max(final[, "t"]) 
    final <- rbind(final[-dim(final)[1], ], model_output) 
  }
  final_vals <- extract_final_values(model_output)
}

noise <- rnorm(length(final[, "Blood_Alb_SO"]), mean = 0, sd = final[, "Blood_Alb_SO"]/10)
noisy_Alb_SO <- final[, "Blood_Alb_SO"] + noise
noisy_Alb_SO[noisy_Alb_SO < 0] <- 0
Alb_SO_data <- data.frame(Time = seq(0, length(noisy_Alb_SO)-1), Alb_SO = noisy_Alb_SO)
plot(Alb_SO_data$Time, Alb_SO_data$Alb_SO)
lines(Alb_SO_data$Time, final[, "Blood_Alb_SO"])

n_mcmc <- 1000
names(parameters_vector) <- c("k_abs", "bioavailability", "sigma", "k_alb", "k_alb_so_init", "k_env", "EC50", "Emax")
reparam <- c(10, 1000, 0.1, 10, 100, 100, 0.1, 0.1)

out <- MCMC_running(number_of_iterations = n_mcmc,
                    parameters_vector = parameters_vector,
                    sd_proposals = rep(0.1, length(parameters_vector)),
                    start_covariance_adaptation = min(1000, round(n_mcmc/2, 0)),
                    Alb_data = NA, 
                    Alb_SO_data = Alb_SO_data,
                    metabolite_availability = "Sulfoxide_Only",
                    model_instance = Albendazole_PK_Model,
                    dose_times = dose_times, 
                    dose_amounts = dose_amounts,
                    output = TRUE,
                    informative_prior = FALSE,
                    reparam = reparam,
                    refresh = 1000)

mcmc_output <- out$MCMC_Output

output <- matrix(nrow = n_mcmc, ncol = 56)
output3 <- matrix(nrow = n_mcmc, ncol = 56)
for (i in 1:n_mcmc) {
  parameters_vector <- mcmc_output[i, ] / reparam
  for (j in 1:length(dose_times)) {
    if (j == 1) {
      mod <- Albendazole_PK_Model(k_abs = parameters_vector[1], sigma = parameters_vector[3], 
                                  k_alb = parameters_vector[4], k_alb_so_init = parameters_vector[5], 
                                  k_env = parameters_vector[6], EC50 = parameters_vector[7], Emax = parameters_vector[8],
                                  gut_1 = (parameters_vector[2]/1e+5) * (dose_amounts[j] * 1e+6), # remember equivalent change in likelihood function required
                                  gut_2 = 0, liver = 0, blood_alb = 0, blood_alb_so = 0, enzyme = 0)
    } else {
      mod <- Albendazole_PK_Model(k_abs = parameters_vector[1], sigma = parameters_vector[3], 
                                  k_alb = parameters_vector[4], k_alb_so_init = parameters_vector[5], 
                                  k_env = parameters_vector[6], EC50 = parameters_vector[7], Emax = parameters_vector[8],
                                  gut_1 = final_vals["Gut_1"] + (parameters_vector[2]/1e+5) * (dose_amounts[j] * 1e+6), 
                                  gut_2 = final_vals["Gut_2"], liver = final_vals["Liver"], blood_alb = final_vals["Blood_Alb"], 
                                  blood_alb_so = final_vals["Blood_Alb_SO"], enzyme = final_vals["Enzyme_conc"])
    }
    if (j != length(dose_times)) {
      times <- seq(0, dose_times[j + 1] - dose_times[j], 1)
    } else {
      times <- seq(0, dose_times[j] - dose_times[j-1] + 15, 1)
    }
    model_output <- mod$run(times)
    if (j == 1) {
      final <- model_output
    } else {
      model_output[, "t"] <- model_output[, "t"] + max(final[, "t"]) 
      final <- rbind(final[-dim(final)[1], ], model_output) 
    }
    final_vals <- extract_final_values(model_output)
  }
  output[i, ] <- final[, "Blood_Alb_SO"]
  output3[i, ] <- final[, "k_alb_so"]
  if (i %% 1000 == 0) {
    print(i)
  }
}

mean <- apply(output, 2, mean)
lower <- apply(output, 2, quantile, 0.025)
upper <- apply(output, 2, quantile, 0.975)
plot(mean, ylim = c(0, max(upper)), type = "l")
lines(lower, col = "red")
lines(upper, col = "red")
points(Alb_SO_data$Alb_SO, pch = 20)

mean3 <- apply(output3, 2, mean)
lower3 <- apply(output3, 2, quantile, 0.025)
upper3 <- apply(output3, 2, quantile, 0.975)
plot(mean3, ylim = c(0, max(upper3)), type = "l")
lines(lower3, col = "red")
lines(upper3, col = "red")

plot(mean3/parameters_vector["k_alb_so_init"], ylim = c(0, max(upper3/parameters_vector["k_alb_so_init"])), type = "l")
lines(lower3/parameters_vector["k_alb_so_init"], col = "red")
lines(upper3/parameters_vector["k_alb_so_init"], col = "red")

hist(mcmc_output[(n_mcmc/2):n_mcmc, 1]/reparam[1])
abline(v = parameters_vector[1], col = "red", lwd = 5)

hist(mcmc_output[(n_mcmc/2):n_mcmc, 2]/reparam[2])
abline(v = parameters_vector[2], col = "red", lwd = 5)

hist(mcmc_output[(n_mcmc/2):n_mcmc, 3]/reparam[3])
abline(v = parameters_vector[3], col = "red", lwd = 5)

hist(mcmc_output[(n_mcmc/2):n_mcmc, 4]/reparam[4])
abline(v = parameters_vector[4], col = "red", lwd = 5)

hist(mcmc_output[(n_mcmc/2):n_mcmc, 5]/reparam[5])
abline(v = parameters_vector[5], col = "red", lwd = 5)

hist(mcmc_output[(n_mcmc/2):n_mcmc, 6]/reparam[6])
abline(v = parameters_vector[6], col = "red", lwd = 5)

hist(mcmc_output[(n_mcmc/2):n_mcmc, 7]/reparam[7])
abline(v = parameters_vector[7], col = "red", lwd = 5)

hist(mcmc_output[(n_mcmc/2):n_mcmc, 8]/reparam[8])
abline(v = parameters_vector[8], col = "red", lwd = 5)

x <- mcmc_output[, c("EC50", "Emax")]
plot(x[, 1]/x[, 2])
abline(h = 1, col = "red", lwd = 3)

hist(x[, 1]/x[, 2])

x <- seq(0, 2, 0.01)
output2 <- matrix(nrow = n_mcmc, ncol = length(x)) 
for (i in 1:n_mcmc) {
  Emax = mcmc_output[i, "Emax"]
  EC50 = mcmc_output[i, "EC50"]
  output2[i, ] <- 1 + (Emax * x)/(EC50 + x)
}

mean2 <- apply(output2, 2, mean)
lower2 <- apply(output2, 2, quantile, 0.025)
upper2 <- apply(output2, 2, quantile, 0.975)
plot(x, mean2, ylim = c(1, max(upper2)), type = "l")
lines(x, lower2, col = "red")
lines(x, upper2, col = "red")

output3 <- matrix(nrow = n_mcmc, ncol = length(x)) 
for (i in 1:n_mcmc) {
  Emax = mcmc_output[i, "Emax"]
  EC50 = mcmc_output[i, "EC50"]
  output2[i, ] <- 1 + (Emax * x)/(EC50 + x)
}

mean2 <- apply(output2, 2, mean)
lower2 <- apply(output2, 2, quantile, 0.025)
upper2 <- apply(output2, 2, quantile, 0.975)
plot(x, mean2, ylim = c(1, max(upper2)), type = "l")
lines(x, lower2, col = "red")
lines(x, upper2, col = "red")

x[4000, ]

params <- mcmc_output[4000, , drop = FALSE]
prior_function(params, informative_prior = FALSE)

cor(mcmc_output)

apply(mcmc_output[10000:20000, ], 2, mean)/reparam
parameters_vector
