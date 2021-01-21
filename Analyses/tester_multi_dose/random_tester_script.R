par(mfrow = c(1, 1))
k <- 10

initial_parameters <- c(k_abs = 0.7, bioavailability = 0.001694134, sigma = 2.622298, k_alb_so = 0.1, k_alb = 0.1)
time_series_number <- k
Single_PK_Dataset <- filter(All_PK_Data, Temporal_ID == k) # filter by time_series_number
dose <- Single_PK_Dataset$Dose_mg[1] * 1000 # need to alter the model so that it takes dose as an input, probably have it accessed as part of Single_PK_Dataset within the MCMC
end_time <- max(Single_PK_Dataset$Time)

Albendazole_Time <- Single_PK_Dataset$Time[Single_PK_Dataset$Metabolite == "Alb"]  
Albendazole_Conc <- Single_PK_Dataset$Converted_Concentration[Single_PK_Dataset$Metabolite == "Alb"]  
Alb_data <- data.frame(Time = Albendazole_Time, Alb = Albendazole_Conc)

Albendazole_Sulfoxide_Time <- Single_PK_Dataset$Time[Single_PK_Dataset$Metabolite == "AlbSO"]  
Albendazole_Sulfoxide_Conc <- Single_PK_Dataset$Converted_Concentration[Single_PK_Dataset$Metabolite == "AlbSO"]  
Alb_SO_data <- data.frame(Time = Albendazole_Sulfoxide_Time, Alb_SO = Albendazole_Sulfoxide_Conc)

if (sum(Single_PK_Dataset$Metabolite == "Alb") >= 1) {
  metabolite_availability <- "Both" 
} else {
  metabolite_availability <- "Sulfoxide_Only"
}

bloop <- MCMC_running(280000, initial_parameters, sd_proposals, start_covariance_adaptation, Alb_data, Alb_SO_data, metabolite_availability, Albendazole_PK_Model, end_time, dose, TRUE, informative_prior)
burn_in <- 0.25
iterations <- number_iterations
chain <- bloop$MCMC_Output[(burn_in*number_iterations):number_iterations, ]

median_k_abs <- median(chain[, 1])
median_bioavailability <- median(chain[, 2])
median_sigma <- median(chain[, 3])
median_k_alb_so <- median(chain[, 4])
median_k_alb <- median(chain[, 5])
model_runner <- Albendazole_PK_Model(k_abs = median_k_abs, bioavailability = median_bioavailability, sigma = median_sigma, k_alb_so = median_k_alb_so, k_alb = median_k_alb, dose = dose) 
times <- seq(0, 100, length.out = 1000)
median_out <- model_runner$run(times)

max <- max(c(as.vector(median_out[, 5:6]), Alb_SO_data$Alb_SO))
plot(median_out[, 1], median_out[, 5], type = "l", col = "#E5005F",  ylim = c(0, max), xlim = c(0, 48), ylab = "Concentration (ng/ml)", xlab = "Time (Hours)", las = 1, lwd = 2)
lines(median_out[, 1], median_out[, 6], type = "l", col = "#7609BA", lwd = 2)
points(Alb_data$Time, Alb_data$Alb, pch = 20, col = "#E5005F", cex = 2)
points(Alb_SO_data$Time, Alb_SO_data$Alb_SO, pch = 20, col = "#7609BA", cex = 2)

parameters_vector <- c(median_k_abs, median_bioavailability, median_sigma, median_k_alb_so, median_k_alb)
loglikelihood_function(parameters_vector, Alb_data, Alb_SO_data, metabolite_availability, Albendazole_PK_Model, end_time, dose)
  
median_k_abs #0.2430469
median_bioavailability # 0.001694134 
median_sigma #2.622298
median_k_alb_so #0.09346704
median_k_alb #0.07739172

k_abs <- 2.255504
bioavailability <- 0.000659322
sigma <- 3.694552
k_alb_so <- 0.03760597
k_alb <- 0.03168715
parameters_vector <- c(k_abs, bioavailability, sigma, k_alb_so, k_alb)
loglikelihood_function(parameters_vector, Alb_data, Alb_SO_data, metabolite_availability, Albendazole_PK_Model, end_time, dose)

par(mfrow = c(1, 1))
model_runner <- Albendazole_PK_Model(k_abs = k_abs, bioavailability = bioavailability, sigma = sigma, k_alb_so = k_alb_so, k_alb = k_alb, dose = 400000)
times <- seq(0, 100, length.out = 1000)
out <- model_runner$run(times)
max <- max(c(as.vector(out[, 5:6]), Alb_SO_data$Alb_SO))
plot(out[, 1], out[, 5], type = "l", col = "#E5005F",  ylim = c(0, max), xlim = c(0, 48), ylab = "Concentration (ng/ml)", xlab = "Time (Hours)", las = 1, lwd = 2)
lines(out[, 1], out[, 6], type = "l", col = "#7609BA", lwd = 2)
points(Alb_data$Time, Alb_data$Alb, pch = 20, col = "#E5005F", cex = 2)
points(Alb_SO_data$Time, Alb_SO_data$Alb_SO, pch = 20, col = "#7609BA", cex = 2)






