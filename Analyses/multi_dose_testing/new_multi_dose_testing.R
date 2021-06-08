library(odin); library(tidyverse)

#source("multi_dose_testing/multi_dose_model_ODIN.R")
source("multi_dose_testing/new_test_multi_dose_model_ODIN.R")

extract_final_values <- function(model_output) {
  end <- dim(model_output)[1]
  Gut_1 <- model_output[end, "Gut_1"]
  Gut_2 <- model_output[end, "Gut_2"]
  Liver <- model_output[end, "Liver"]
  Blood_Alb <- model_output[end, "Blood_Alb"]
  Blood_Alb_SO <- model_output[end, "Blood_Alb_SO"]
  Enzyme_conc <- model_output[end, "Enzyme_conc"]
  return(c(Gut_1, Gut_2, Liver, Blood_Alb, Blood_Alb_SO, Enzyme_conc))
}

# Generating Dummy Data 
bioavailability <- 0.001
dose <- 400
mod <- Albendazole_PK_Model(k_abs = 0.5, sigma = 15, k_alb = 0.2, k_alb_so_init = 0.13, 
                            k_env = 0.05, Emax = 5, EC50 = 3,
                            gut_1 = (bioavailability/1e+5) * (dose * 1e+6), gut_2 = 0, liver = 0, 
                            blood_alb = 0, blood_alb_so = 0, enzyme = 0)
output1 <- mod$run(t = seq(0, 20, 1))
final1 <- extract_final_values(output1)
mod2 <- Albendazole_PK_Model(k_abs = 0.5, sigma = 15, k_alb = 0.2, k_alb_so_init = 0.13, 
                             k_env = 0.05, Emax = 5, EC50 = 3,
                             gut_1 = final1["Gut_1"] + (bioavailability/1e+5) * (dose * 1e+6), gut_2 = final1["Gut_2"], liver = final1["Liver"], 
                             blood_alb = final1["Blood_Alb"], blood_alb_so = final1["Blood_Alb_SO"], enzyme = final1["Enzyme_conc"])
output2 <- mod2$run(t = seq(0, 20, 1))
final <- rbind(output1, output2[-1, ])

plot(final[, "Enzyme_conc"], type = "l")
plot(final[, "Blood_Alb_SO"], type = "l")
plot(final[, "k_alb_so"], type = "l", ylim = c(0, max(final[, "k_alb_so"])))
plot(final[, "k_alb_so"]/0.13, type = "l")
plot(final[, "x"], type = "l")

plot(final[, "Blood_Alb_SO"], type = "l")
plot(final[, "Blood_Alb"], type = "l")
plot(final[, "Gut_1"], type = "l")
plot(final[, "Gut_2"], type = "l")
plot(final[, "Liver"], type = "l")


Emax = 0
EC50 = 3
x <- seq(0, 10, 0.1)
y <- 1 + (Emax * x)/(EC50 + x)
plot(x, y)

final_noise_inc <- cbind(final, final[, "Blood_Alb_SO"] + rnorm(n = length(final[, "Blood_Alb_SO"]), 
                                                      mean = 0,
                                                      sd = final[, "Blood_Alb_SO"])/4)
colnames(final_noise_inc)[length(colnames(final_noise_inc))] <- "noisy"
final_noise_inc[, "noisy"][final_noise_inc[, "noisy"] <= 0] <- 0

plot(final_noise_inc[, "noisy"])
lines(final_noise_inc[, "Blood_Alb_SO"], type = "l")




parameters_vector <- c(0.5, 0.001, 15, 0.2, 0.13, 0.2, 0.05)
dose_amounts <- c(400, 400, 400)
dose_times <- c(0, 25, 50)
Alb_SO_data <- data.frame(temp = 0, Time = 75)
# final <- matrix(NA, ncol = 7, nrow = 0)
# colnames(final) <- c("t", "Gut_1", "Gut_2", "Liver", "Blood_Alb", "Blood_Alb_SO", "Enzyme_conc")

for (i in 1:length(dose_times)) {
  
  if (i == 1) {
    mod <- Albendazole_PK_Model(k_abs = parameters_vector[1], sigma = parameters_vector[3], 
                                k_alb = parameters_vector[4], k_alb_so_init = parameters_vector[5], 
                                auto_ind = parameters_vector[6], k_env = parameters_vector[7],
                                gut_1 = (parameters_vector[2]/1e+5) * (dose_amounts[i] * 1e+6), 
                                gut_2 = 0, liver = 0, 
                                blood_alb = 0, blood_alb_so = 0, enzyme = 0)
  } else {
    mod <- Albendazole_PK_Model(k_abs = parameters_vector[1], sigma = parameters_vector[3], 
                                k_alb = parameters_vector[4], k_alb_so_init = parameters_vector[5], 
                                auto_ind = parameters_vector[6], k_env = parameters_vector[7],
                                gut_1 = final_vals["Gut_1"] + (parameters_vector[2]/1e+5) * (dose_amounts[i] * 1e+6), 
                                gut_2 = final_vals["Gut_2"], liver = final_vals["Liver"], 
                                blood_alb = final_vals["Blood_Alb"], blood_alb_so = final_vals["Blood_Alb_SO"], enzyme = final_vals["Enzyme_conc"])
  }
  
  if (i != length(dose_times)) {
    times <- seq(0, dose_times[i + 1] - dose_times[i], 1)
  } else {
    times <- seq(0, max(Alb_SO_data$Time) - dose_times[i], 1)
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

plot(final[, "t"], final[, "Blood_Alb_SO"])
plot(final[, "t"], final[, "Enzyme_conc"])
