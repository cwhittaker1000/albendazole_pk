# Loading In Required Libraries 
library(tictoc); library(deSolve); library(MALDIquant); library(mvtnorm); library(MASS); 
library(tmvtnorm); library(odin); library(dplyr); library(naniar)

# Loading in ODIN Model Instance
source("Models/Alb_PK_Single_Dose_Odin_Model.R")
source("Functions/single_dose_MCMC_functions.R")

# Importing and Processing the Data - 60 (1-60) single dose time-series; 
#                                     13 (61-73) multiple dose, single time-series; 
#                                      9 (74-82) multiple dose, multiple time-series
single_PK_Data <- read.csv("Data/single_dose_AlbPK_data.csv", stringsAsFactors = FALSE)
multiple_PK_Data <- read.csv("Data/multiple_dose_AlbPK_data.csv", stringsAsFactors = FALSE)
overall <- rbind(single_PK_Data, multiple_PK_Data) 

# Defining the MCMC Inputs
number_iterations <- 40000
initial_parameters <- c(k_abs = 2, bioavailability = 0.01, sigma = 15, k_alb = 0.2, k_alb_so = 0.1300)
sd_proposals <- c(0.002, 0.002, 0.002, 0.002, 0.002)
reparam <- c(1, 100, 0.1, 10, 10) 
start_covariance_adaptation <- 5000
informative_prior <- FALSE
burnin <- min(number_iterations/2, start_covariance_adaptation * 2)

# Dataframe for Summary Model Outputs
num_TS <- max(overall$Temporal_ID)
results <- data.frame(median_k_alb_so = rep(NA, num_TS), median_bioavailability = rep(NA, num_TS), median_AUC = rep(NA, num_TS), median_Cmax = rep(NA, num_TS),
                      mean_k_alb_so = rep(NA, num_TS), mean_bioavailability = rep(NA, num_TS), mean_AUC = rep(NA, num_TS), mean_Cmax = rep(NA, num_TS),
                      sex = rep(NA_character_, num_TS), sex_ratio = rep(NA, num_TS), state = rep(NA_character_, num_TS), dose = rep(NA, num_TS),
                      dose_single = rep(NA, num_TS), dose_binary = rep(NA_character_, num_TS), infection = rep(NA_character_, num_TS), 
                      age = rep(NA, num_TS), age_binary = rep(NA_character_, num_TS), weight = rep(NA, num_TS), drug = rep(NA_character_, num_TS), 
                      dosing = character(num_TS), number = rep(NA, num_TS), dose_data_availabilty = rep(NA, num_TS), stringsAsFactors = FALSE)

# Looping Over the Dataset and Fitting Each of the Datasets
for (k in 1:num_TS) {
  
  # Selecting individual time series for fitting 
  time_series_number <- k
  Single_PK_Dataset <- filter(overall, Temporal_ID == k) # filter by time_series_number
  
  # Extracting dose and time period dataset spans
  dosing <- unique(Single_PK_Dataset$Dosing)
  if (dosing == "Single") {
    dose_info <- data.frame(amount = unique(Single_PK_Dataset$Dose_Amount), times = 0)
  } else if (dosing == "Multiple") {
    dose_info <- data.frame(amount = unique(Single_PK_Dataset$Dose_Amount), times = eval(parse(text = unique(Single_PK_Dataset$Dose_Timing))))
  } else {
    stop("Error in single vs multiple")
  }
  if (unique(dose_info$amount) != length(dose_info$amount)) {
    amount <- as.numeric(unique(dose_info$amount) * unlist(unname(table(dose_info$times))))
    dose_info <- data.frame(amount = amount, times = unique(dose_info$times))
  }
  
  # Converting Multiple Dose Time-Series Where We have Starting Dose Time-Series Available to Single-Dose Time-Series
  data_available <- unique(Single_PK_Dataset$Data)
  if (dosing == "Multiple" & data_available == "Multiple") {
    time_2nd_dose <- dose_info$times[2]
    Single_PK_Dataset <- Single_PK_Dataset %>%
      filter(Time < time_2nd_dose) ## go back and check this is correct for each
    dose_info <- dose_info[1, ]
    amount <- dose_info$amount[1]
    dosing <- "Single"
  } 
  
  # Extracting Albendazole and Albendazole Sulfoxide data if present
  Albendazole_Time <- Single_PK_Dataset$Time[Single_PK_Dataset$Metabolite == "Alb"]  
  Albendazole_Conc <- Single_PK_Dataset$Converted_Concentration[Single_PK_Dataset$Metabolite == "Alb"]  
  Alb_data <- data.frame(Time = Albendazole_Time, Alb = Albendazole_Conc)
  
  Albendazole_Sulfoxide_Time <- Single_PK_Dataset$Time[Single_PK_Dataset$Metabolite == "AlbSO"]  
  Albendazole_Sulfoxide_Conc <- Single_PK_Dataset$Converted_Concentration[Single_PK_Dataset$Metabolite == "AlbSO"]  
  Alb_SO_data <- data.frame(Time = Albendazole_Sulfoxide_Time, Alb_SO = Albendazole_Sulfoxide_Conc)
  
  # Defining which dataset(s) to fit to
  if (sum(Single_PK_Dataset$Metabolite == "Alb") >= 1) {
    metabolite_availability <- "Both" 
    time_increment <- min(c(1, diff(Alb_data$Time), diff(Alb_SO_data$Time)))
  } else {
    metabolite_availability <- "Sulfoxide_Only"
    time_increment <- min(c(1, diff(Alb_SO_data$Time)))
  }
  
  # Running the MCMC and saving the output
  run_MCMC <- single_dose_MCMC_running (number_of_iterations = number_iterations, 
                                        parameters_vector = initial_parameters, 
                                        sd_proposals = sd_proposals, 
                                        start_covariance_adaptation = start_covariance_adaptation, 
                                        Alb_data = Alb_data, 
                                        Alb_SO_data = Alb_SO_data, 
                                        metabolite_availability = metabolite_availability, 
                                        model_instance = Albendazole_PK_Model, 
                                        dose_info = dose_info, 
                                        dosage = dosing, 
                                        output = FALSE, 
                                        informative_prior = informative_prior, 
                                        reparam = reparam, 
                                        refresh = 1000,
                                        time_increment = time_increment)
  
  # Saving Output
  if (informative_prior) {
    prior_string <- "Inf"
  } else {
    prior_string <- "Uninf"
  }
  saveRDS(list(number_iterations = number_iterations, info = Single_PK_Dataset, mcmc_output = run_MCMC, Alb_data = Alb_data, Alb_SO_data = Alb_SO_data), 
          paste0("Outputs/SingleDose_AlbPK_ModelFitting_Outputs/TS", k, "_", prior_string, "Prior", "_", dosing, "Dose.rds"))
  
  # Plotting Output
  alb_so <- run_MCMC$Alb_SO[burnin:number_iterations, ]
  mean <- apply(alb_so, 2, mean)
  lower <- apply(alb_so, 2, quantile, 0.025)
  upper <- apply(alb_so, 2, quantile, 0.975)
  plot(run_MCMC$Times, mean, type = "l", col = "#7609BA", ylab = "Concentration (ng/ml)", xlab = "Time (Hours)", las = 1, lwd = 2, main = paste0("Time Series ", k), ylim = c(0, max(c(upper, Alb_SO_data$Alb_SO))))
  points(Alb_SO_data$Time, Alb_SO_data$Alb_SO, pch = 20, col = "#7609BA", cex = 2)
  polygon(c(run_MCMC$Times, rev(run_MCMC$Times)), c(lower, rev(upper)), col = adjustcolor("#7609BA", alpha.f = 0.2), border = NA)
  
  # Generating and Saving Summary Statistics for Each Fitted Output 
  times <- seq(0, 100, 0.01)
  dose <- dose_info$amount
  
  ## Extracting and Running Model With Median Outputs
  median_k_abs <- median(run_MCMC$MCMC_Output[burnin:number_iterations, "k_abs"]) / reparam[1]
  median_bioavailability <- median(run_MCMC$MCMC_Output[burnin:number_iterations, "bioavailability"]) / reparam[2]
  median_sigma <- median(run_MCMC$MCMC_Output[burnin:number_iterations, "sigma"]) / reparam[3]
  median_k_alb <- median(run_MCMC$MCMC_Output[burnin:number_iterations, "k_alb"]) / reparam[4]
  median_k_alb_so <- median(run_MCMC$MCMC_Output[burnin:number_iterations, "k_alb_so"]) / reparam[5]
  median_model <- Albendazole_PK_Model(k_abs = median_k_abs, sigma = median_sigma, k_alb_so = median_k_alb_so, k_alb = median_k_alb, 
                                       gut_1 = (median_bioavailability/1e+5) * (dose * 1e+6), gut_2 = 0, liver = 0, blood_alb = 0, blood_alb_so = 0)
  median_out <- median_model$run(times)
  median_alb_so_PK <- median_out[, 6]
  median_AUC <- DescTools::AUC(times, median_alb_so_PK)
  median_Cmax <- max(median_alb_so_PK)
  results$median_k_alb_so[k] <- median_k_alb_so
  results$median_bioavailability[k] <- median_bioavailability
  results$median_AUC[k] <- median_AUC
  results$median_Cmax[k] <- median_Cmax
  
  ## Mean Outputs
  mean_k_abs <- mean(run_MCMC$MCMC_Output[burnin:number_iterations, "k_abs"]) / reparam[1]
  mean_bioavailability <- mean(run_MCMC$MCMC_Output[burnin:number_iterations, "bioavailability"]) / reparam[2]
  mean_sigma <- mean(run_MCMC$MCMC_Output[burnin:number_iterations, "sigma"]) / reparam[3]
  mean_k_alb <- mean(run_MCMC$MCMC_Output[burnin:number_iterations, "k_alb"]) / reparam[4]
  mean_k_alb_so <- mean(run_MCMC$MCMC_Output[burnin:number_iterations, "k_alb_so"]) / reparam[5]
  mean_model <- Albendazole_PK_Model(k_abs = mean_k_abs, sigma = mean_sigma, k_alb_so = mean_k_alb_so, k_alb = mean_k_alb, 
                                     gut_1 = (mean_bioavailability/1e+5) * (dose * 1e+6), gut_2 = 0, liver = 0, blood_alb = 0, blood_alb_so = 0)
  mean_out <- mean_model$run(times)
  mean_alb_so_PK <- mean_out[, 6]
  mean_AUC <- DescTools::AUC(times, mean_alb_so_PK)
  mean_Cmax <- max(mean_alb_so_PK)
  results$mean_k_alb_so[k] <- mean_k_alb_so # change to 1 over this quantity i.e. half life
  results$mean_bioavailability[k] <- mean_bioavailability
  results$mean_AUC[k] <- mean_AUC
  results$mean_Cmax[k] <- mean_Cmax  
  
  # Get metadata for each timeseries
  meta <- Single_PK_Dataset
  results$sex[k] <- unique(meta$Sex_Binary)
  results$sex_ratio[k] <- unique(meta$Sex_Ratio)
  results$state[k] <- unique(meta$Feeding_State)
  results$dose[k] <- unique(meta$Dose_mg_daily)
  results$dose_single[k] <- unique(meta$Dose_Amount)
  results$dose_binary[k] <- unique(meta$Dose_Binary)
  results$infection[k] <- unique(meta$Disease_Binary)
  results$age[k] <- unique(meta$Age_Cont)
  results$age_binary[k] <- unique(meta$Age_Binary)
  results$weight[k] <- unique(meta$Weight)
  results$dosing[k] <- unique(meta$Dosing)
  results$dose_data_availability[k] <- unique(meta$Data)
  results$number[k] <- unique(meta$Number)
  results$drug[k] <- unique(meta$Drug_Binary)
  
  print(k) 
  
}


# Plotting all the raw model fitting results
m <- matrix(c(1, 1, 2, 3, 1, 1, 4, 5, 1, 1, 6, 7), nrow = 3, ncol = 4, byrow = TRUE)
layout(m)
for (i in 1:max(overall$Temporal_ID)) {
  
  # Loading in results
  Single_PK_Dataset <- filter(overall, Temporal_ID == i) # filter by time_series_number
  dosing <- unique(Single_PK_Dataset$Dosing)
  temp <- readRDS(paste0("Outputs/SingleDose_SimpleModel_MCMC_Outputs/TS", i, "_prior", prior_string, "_", tolower(dosing), "Dose_simpleModel_.rds"))
  alb_so <- temp$mcmc_output$Alb_SO[max(temp$number_iterations/2, start_covariance_adaptation):temp$number_iterations, ]
  mean_alb_so <- apply(alb_so, 2, mean)
  lower_alb_so <- apply(alb_so, 2, quantile, 0.025)
  upper_alb_so <- apply(alb_so, 2, quantile, 0.975)
  
  alb <- temp$mcmc_output$Alb[max(temp$number_iterations/2, start_covariance_adaptation):temp$number_iterations, ]
  mean_alb <- apply(alb, 2, mean)
  lower_alb <- apply(alb, 2, quantile, 0.025)
  upper_alb <- apply(alb, 2, quantile, 0.975)
  
  # Plotting model outputs
  plot(temp$mcmc_output$Times, mean_alb_so, type = "l", col = "#7609BA", ylab = "Concentration (ng/ml)", xlab = "Time (Hours)", las = 1, lwd = 2, 
       main = paste0("Time Series ", i), ylim = c(0, max(upper_alb_so, temp$info$Converted_Concentration[temp$info$Metabolite == "AlbSO"])))
  points(temp$info$Time[temp$info$Metabolite == "AlbSO"], temp$info$Converted_Concentration[temp$info$Metabolite == "AlbSO"], pch = 20, col = "#7609BA", cex = 2)
  polygon(c(temp$mcmc_output$Times, rev(temp$mcmc_output$Times)), c(lower_alb_so, rev(upper_alb_so)), col = adjustcolor("#7609BA", alpha.f = 0.2), border = NA)
  if("Alb" %in% temp$info$Metabolite) {
    lines(temp$mcmc_output$Times, mean_alb, type = "l", col = "#E5005F", ylab = "Concentration (ng/ml)", xlab = "Time (Hours)", las = 1, lwd = 2, main = paste0("Time Series ", i))
    points(temp$info$Time[temp$info$Metabolite == "Alb"], temp$info$Converted_Concentration[temp$info$Metabolite == "Alb"], pch = 20, col = "#E5005F", cex = 2)
    polygon(c(temp$mcmc_output$Times, rev(temp$mcmc_output$Times)), c(lower_alb, rev(upper_alb)), col = adjustcolor("#E5005F", alpha.f = 0.2), border = NA)
  }
  
  colours <- rainbow(5)
  naming <- colnames(temp$mcmc_output$MCMC_Output)
  for (j in 1:dim(temp$mcmc_output$MCMC_Output)[2]) {
    plot(temp$mcmc_output$MCMC_Output[max(temp$number_iterations/2, start_covariance_adaptation):temp$number_iterations, j],
         type = "l", col = colours[j], main = naming[j], ylab = "")
  }
  plot.new()
  
}

# Regression Results - need to decide how and what to include here given data sparsity
miss_var_summary(results)

# Sex
sex <- results[!is.na(results$sex) & results$sex != "Unclear" & results$sex != "Female", ]
summary(lm(median_bioavailability ~ sex, data = sex))
summary(lm(median_k_alb_so ~ sex, data = sex))
summary(lm(median_Cmax ~ sex, data = sex))

sex_ratio <- results[!is.na(results$sex_ratio), ]
summary(lm(median_bioavailability ~ sex_ratio, data = sex_ratio))
summary(lm(median_k_alb_so ~ sex_ratio, data = sex_ratio))
summary(lm(median_Cmax ~ sex_ratio, data = sex_ratio))

# Fasting State
state <- results[!is.na(results$state) & results$state != "Mixture" & results$state != "Unclear", ]
summary(lm(median_bioavailability ~ state, data = state))
summary(lm(median_k_alb_so ~ state, data = state))
summary(lm(median_Cmax ~ state, data = state))

# Dose
dose_single <- results[!is.na(results$dose_single), ]
summary(lm(median_bioavailability ~ dose_single, data = dose_single))
summary(lm(median_k_alb_so ~ dose_single, data = dose_single))
summary(lm(median_Cmax ~ dose_single, data = dose_single))

dose_binary <- results[!is.na(results$dose_binary), ]
summary(lm(median_bioavailability ~ dose_binary, data = dose_binary))
summary(lm(median_k_alb_so ~ dose_binary, data = dose_binary))
summary(lm(median_Cmax ~ dose_binary, data = dose_binary))

# Infection
infection <- results[!is.na(results$infection) & results$infection != "Mixture" & results$infection != "None", ]
summary(lm(median_bioavailability ~ infection, data = infection))
summary(lm(median_k_alb_so ~ infection, data = infection))
summary(lm(median_Cmax ~ infection, data = infection))

# Age
age <- results[!is.na(results$age), ]
summary(lm(median_bioavailability ~ age, data = age))
summary(lm(median_k_alb_so ~ age, data = age))
summary(lm(median_Cmax ~ age, data = age))

age_binary <- results[!is.na(results$age_binary) & results$age_binary != "Unclear", ]
summary(lm(median_bioavailability ~ age_binary, data = age_binary))
summary(lm(median_k_alb_so ~ age_binary, data = age_binary))
summary(lm(median_Cmax ~ age_binary, data = age_binary))

# Weight
weight <- results[!is.na(results$weight), ]
summary(lm(median_bioavailability ~ weight, data = weight))
summary(lm(median_k_alb_so ~ weight, data = weight))
summary(lm(median_Cmax ~ weight, data = weight))

# Drug
drug <- results[!is.na(results$drug), ]
summary(lm(median_bioavailability ~ drug, data = drug))
summary(lm(median_k_alb_so ~ drug, data = drug))
summary(lm(median_Cmax ~ drug, data = drug))

# Dosing
dosing <- results[!is.na(results$dosing), ]
summary(lm(median_bioavailability ~ dosing, data = dosing))
summary(lm(median_k_alb_so ~ dosing, data = dosing))
summary(lm(median_Cmax ~ dosing, data = dosing))

dosing <- results[!is.na(results$dosing), ]
summary(lm(median_bioavailability ~ dosing, data = dosing))
summary(lm(median_k_alb_so ~ dosing, data = dosing))
summary(lm(median_Cmax ~ dosing, data = dosing))



colnames(results)
summary(lm(median_AUC ~ dosing + sex + state + infection, data = results))

summary(lm(median_bioavailability ~ dose_binary + dosing + state + infection, data = results))
summary(lm(median_k_alb_so ~ dose_binary + dosing + state + infection, data = results))
summary(lm(median_Cmax ~ dose_binary + dosing + state + infection, data = results))


summary(lm(median_bioavailability ~ dose_binary + dosing + sex + state + infection, data = results))
summary(lm(median_k_alb_so ~ dose_binary + dosing + sex + state + infection, data = results))
summary(lm(median_Cmax ~ dose_binary + dosing + sex + state + infection, data = results))

summary(lm(median_Cmax ~ dosing + sex + state, data = results))


summary(lm(median_bioavailability ~ sex + state + dose_binary + drug + age_binary + infection, data = results))
summary(lm(median_k_alb_so ~ sex + state + dose_binary + drug + age_binary + infection, data = results))
summary(lm(median_AUC ~ sex + state + dose_binary + drug + age_binary + infection, data = results))
summary(lm(median_Cmax ~ sex + state + dose_binary + drug + age_binary + infection, data = results))

x <- results[!is.na(results$sex_ratio) & !is.na(results$state)  & !is.na(results$dose_binary) &
               !is.na(results$drug) & !is.na(results$age_binary) & !is.na(results$infection), ]

sum(!is.na(results$sex))
sum(!is.na(results$sex_ratio))
sum(!is.na(results$state))
sum(!is.na(results$dose_binary))
sum(!is.na(results$drug))
sum(!is.na(results$age))
sum(!is.na(results$age_binary))
sum(!is.na(results$weight))
sum(!is.na(results$infection))


summary(lm(mean_bioavailability ~ sex_ratio + state + dose_binary + drug + age + weight + infection, data = results))
summary(lm(mean_k_alb_so ~ sex_ratio + state + dose_binary + drug + age + weight + infection, data = results))
summary(lm(mean_AUC ~ sex_ratio + state + dose_binary + drug + age + weight + infection, data = results))
summary(lm(mean_Cmax ~ sex_ratio + state + dose_binary + drug + age + weight + infection, data = results))

# Generating median outputs for plotting
times <- seq(0, 50, 0.1)
median_output <- matrix(nrow = max(overall$Temporal_ID), ncol = length(times))
for (i in 1:max(overall$Temporal_ID)) {
  
  # Loading in results
  Single_PK_Dataset <- filter(overall, Temporal_ID == i) # filter by time_series_number
  dosing <- unique(Single_PK_Dataset$Dosing)
  temp <- readRDS(paste0("Outputs/SingleDose_SimpleModel_MCMC_Outputs/TS", i, "_prior", prior_string, "_", tolower(dosing), "Dose_simpleModel_.rds"))

  # Running model with median parameter estimates
  median_k_abs <- median(temp$mcmc_output$MCMC_Output[max(temp$number_iterations/2, start_covariance_adaptation):temp$number_iterations, "k_abs"]) / reparam[1]
  median_bioavailability <- median(temp$mcmc_output$MCMC_Output[max(temp$number_iterations/2, start_covariance_adaptation):temp$number_iterations, "bioavailability"]) / reparam[2]
  median_sigma <- median(temp$mcmc_output$MCMC_Output[max(temp$number_iterations/2, start_covariance_adaptation):temp$number_iterations, "sigma"]) / reparam[3]
  median_k_alb <- median(temp$mcmc_output$MCMC_Output[max(temp$number_iterations/2, start_covariance_adaptation):temp$number_iterations, "k_alb"]) / reparam[4]
  median_k_alb_so <- median(temp$mcmc_output$MCMC_Output[max(temp$number_iterations/2, start_covariance_adaptation):temp$number_iterations, "k_alb_so"]) / reparam[5]
  median_model <- Albendazole_PK_Model(k_abs = median_k_abs, sigma = median_sigma, 
                                       k_alb_so = median_k_alb_so, k_alb = median_k_alb, 
                                       gut_1 = (median_bioavailability/1e+5) * (dose * 1e+6), 
                                       gut_2 = 0, liver = 0, blood_alb = 0, blood_alb_so = 0)
  median_out <- median_model$run(times)
  median_alb_so_PK <- median_out[, 6]
  median_output[i, ] <- median_alb_so_PK
  
  print(i)
  
}

# Plotting dosing results
par(mfrow = c(1, 1))
dosing <- factor(results$dosing)
colours <- c("#4E90CE", "#CE4E61")
indices <- which(!is.na(dosing))
mean_single <- apply(median_output[dosing == "Single", ], 2, mean)
mean_multiple <- apply(median_output[dosing == "Multiple", ], 2, mean)
for (i in indices) {
  if (i == 1) {
    plot(times, median_output[i, ], type = "l", ylim = c(0, 1000), las = 1, xlab = "Time Since Treatment (Hours)", ylab = "Concentration (ng/ml)",
         col = ifelse(dosing[i] == "Single", adjustcolor(colours[1], alpha = 0.2), adjustcolor(colours[2], alpha = 0.2)))
  } else {
    lines(times, median_output[i, ], col = ifelse(dosing[i] == "Single", adjustcolor(colours[1], alpha = 0.2), adjustcolor(colours[2], alpha = 0.2)))
  }
}
lines(times, mean_single, lwd = 3, col = colours[1])
lines(times, mean_multiple, lwd = 3, col = colours[2])
legend("topright", inset = 0.02, legend = c(paste0("Single (n = ", sum(dosing == "Single"), ")"), 
                                            paste0("Multiple (n = ", sum(dosing == "Multiple"), ")")), 
       title = "Dosing", col = c(colours[1], colours[2]), lty = 1, lwd = 5, cex = 1, title.adj = 0.5)


# Plotting sex results
par(mfrow = c(1, 1))
sex <- factor(results$sex)
colours <- c("#4E90CE", "#CE4E61")
indices <- which(!is.na(sex) & sex != "Unclear")
mean_men <- apply(median_output[sex == "Male", ], 2, mean)
mean_mixture <- apply(median_output[sex == "Mixture", ], 2, mean)
for (i in indices) {
  if (i == 1) {
    plot(times, median_output[i, ], type = "l", ylim = c(0, 1000), las = 1, xlab = "Time Since Treatment (Hours)", ylab = "Concentration (ng/ml)",
         col = ifelse(sex[i] == "Male", adjustcolor(colours[1], alpha = 0.2), adjustcolor(colours[2], alpha = 0.2)))
  } else {
    lines(times, median_output[i, ], col = ifelse(sex[i] == "Male", adjustcolor(colours[1], alpha = 0.2), adjustcolor(colours[2], alpha = 0.2)))
  }
}
lines(times, mean_men, lwd = 3, col = colours[1])
lines(times, mean_mixture, lwd = 3, col = colours[2])
legend("topright", inset = 0.02, legend = c(paste0("Males (n = ", sum(sex == "Male"), ")"), paste0("Mixture (n = ", sum(sex == "Mixture"), ")")), 
       title = "Sex", col = c(colours[1], colours[2]), lty = 1, lwd = 5, cex = 1, title.adj = 0.5)

# Plotting state results
state <- factor(results$state)
colours <- c("#4BBC00","#F28123")
indices <- which(!is.na(state) & state != "Unclear" & state != "Mixture")
mean_fatty_meal <- apply(median_output[state == "Fatty_meal", ], 2, mean)
mean_fasted <- apply(median_output[state == "Fasted", ], 2, mean)
for (i in indices) {
  if (i == 1) {
    plot(times, median_output[i, ], type = "l", ylim = c(0, 1000), las = 1, xlab = "Time Since Treatment (Hours)", ylab = "Concentration (ng/ml)",
         col = ifelse(state[i] == "Fatty_meal", adjustcolor(colours[1], alpha = 0.2), adjustcolor(colours[2], alpha = 0.2)))
  } else {
    lines(times, median_output[i, ], col = ifelse(state[i] == "Fatty_meal", adjustcolor(colours[1], alpha = 0.2), adjustcolor(colours[2], alpha = 0.2)))
  }
}
lines(times, mean_fatty_meal, lwd = 3, col = colours[1])
lines(times, mean_fasted, lwd = 3, col = colours[2])
legend("topright", inset = 0.02, legend = c(paste0("Fatty Meal (n = ", sum(state == "Fatty_meal"), ")"), paste0("Fasted (n = ", sum(state == "Fasted"), ")")), 
       title = "State", col = c(colours[1], colours[2]), lty = 1, lwd = 5, cex = 1, title.adj = 0.5)

# Plotting dose results
dose <- factor(results$dose_binary)
colours <- c("#3ED6BA", "#725FAD")
indices <- which(!is.na(dose))
mean_high <- apply(median_output[dose == "High", ], 2, mean)
mean_low <- apply(median_output[dose == "Low", ], 2, mean)
for (i in indices) {
  if (i == 1) {
    plot(times, median_output[i, ], type = "l", ylim = c(0, 1000), las = 1, xlab = "Time Since Treatment (Hours)", ylab = "Concentration (ng/ml)",
         col = ifelse(dose[i] == "High", adjustcolor(colours[1], alpha = 0.2), adjustcolor(colours[2], alpha = 0.2)))
  } else {
    lines(times, median_output[i, ], col = ifelse(dose[i] == "High", adjustcolor(colours[1], alpha = 0.2), adjustcolor(colours[2], alpha = 0.2)))
  }
}
lines(times, mean_high, lwd = 3, col = colours[1])
lines(times, mean_low, lwd = 3, col = colours[2])
legend("topright", inset = 0.02, legend = c(paste0("High (n = ", sum(dose == "High"), ")"), paste0("Low (n = ", sum(dose == "Low"), ")")), 
       title = "State", col = c(colours[1], colours[2]), lty = 1, lwd = 5, cex = 1, title.adj = 0.5)

# Plotting infection results
infection <- factor(results$infection)
colours <- c("#F7C911", "#CE5A90")
indices <- which(!is.na(infection) & infection != "Mixture")
mean_infected <- apply(median_output[infection == "Infected", ], 2, mean)
mean_healthy <- apply(median_output[infection == "Healthy", ], 2, mean)
for (i in indices) {
  if (i == 1) {
    plot(times, median_output[i, ], type = "l", ylim = c(0, 1000), las = 1, xlab = "Time Since Treatment (Hours)", ylab = "Concentration (ng/ml)",
         col = ifelse(infection[i] == "Healthy", adjustcolor(colours[1], alpha = 0.2), adjustcolor(colours[2], alpha = 0.2)))
  } else {
    lines(times, median_output[i, ], col = ifelse(infection[i] == "Healthy", adjustcolor(colours[1], alpha = 0.2), adjustcolor(colours[2], alpha = 0.2)))
  }
}
lines(times, mean_healthy, lwd = 3, col = colours[1])
lines(times, mean_infected, lwd = 3, col = colours[2])
legend("topright", inset = 0.02, legend = c(paste0("Healthy (n = ", sum(infection == "Healthy"), ")"), paste0("Infected (n = ", sum(infection == "Infected"), ")")), 
       title = "Infection Status", col = c(colours[1], colours[2]), lty = 1, lwd = 5, cex = 1, title.adj = 0.5)

# Plotting drug results
drug <- factor(results$drug)
colours <- c("#ADA9A9", "#D62C2C")
indices <- which(!is.na(drug))
mean_none <- apply(median_output[drug == "None", ], 2, mean)
mean_drug <- apply(median_output[drug == "Yes", ], 2, mean)
for (i in indices) {
  if (i == 1) {
    plot(times, median_output[i, ], type = "l", ylim = c(0, 1000), las = 1, xlab = "Time Since Treatment (Hours)", ylab = "Concentration (ng/ml)",
         col = ifelse(drug[i] == "None", adjustcolor(colours[1], alpha = 0.2), adjustcolor(colours[2], alpha = 0.2)))
  } else {
    lines(times, median_output[i, ], col = ifelse(drug[i] == "None", adjustcolor(colours[1], alpha = 0.2), adjustcolor(colours[2], alpha = 0.2)))
  }
}
lines(times, mean_none, lwd = 3, col = colours[1])
lines(times, mean_drug, lwd = 3, col = colours[2])
legend("topright", inset = 0.02, legend = c(paste0("None (n = ", sum(drug == "None"), ")"), paste0("Drug (n = ", sum(drug == "Yes"), ")")), 
       title = "Drug Co-Administration", col = c(colours[1], colours[2]), lty = 1, lwd = 5, cex = 1, title.adj = 0.5)

# Plotting age results
age <- factor(results$age_binary)
colours <- c("#E0A02A", "#ED7BEB")
indices <- which(!is.na(age) & age != "Unclear")
mean_adults <- apply(median_output[age == "Adult", ], 2, mean)
mean_children <- apply(median_output[age == "Children", ], 2, mean)
for (i in indices) {
  if (i == 1) {
    plot(times, median_output[i, ], type = "l", ylim = c(0, 1000), las = 1, xlab = "Time Since Treatment (Hours)", ylab = "Concentration (ng/ml)",
         col = ifelse(age[i] == "Adult", adjustcolor(colours[1], alpha = 0.2), adjustcolor(colours[2], alpha = 0.2)))
  } else {
    lines(times, median_output[i, ], col = ifelse(age[i] == "Adult", adjustcolor(colours[1], alpha = 0.2), adjustcolor(colours[2], alpha = 0.2)))
  }
}
lines(times, mean_adults, lwd = 3, col = colours[1])
lines(times, mean_children, lwd = 3, col = colours[2])
legend("topright", inset = 0.02, legend = c(paste0("Adults (n = ", sum(age == "Adult"), ")"), paste0("Children (n = ", sum(age == "Children"), ")")), 
       title = "Age", col = c(colours[1], colours[2]), lty = 1, lwd = 5, cex = 1, title.adj = 0.5)

# Plotting weight results
weight <- results$weight > 60
colours <- c("#32720C", "#CCB5FF")
indices <- which(!is.na(weight))
mean_heavy <- apply(median_output[weight[indices], ], 2, mean)
mean_light <- apply(median_output[!weight[indices], ], 2, mean)
for (i in indices) {
  if (i == 1) {
    plot(times, median_output[i, ], type = "l", ylim = c(0, 1000), las = 1, xlab = "Time Since Treatment (Hours)", ylab = "Concentration (ng/ml)",
         col = ifelse(weight[indices[i]], adjustcolor(colours[1], alpha = 0.2), adjustcolor(colours[2], alpha = 0.2)))
  } else {
    lines(times, median_output[i, ], col = ifelse(weight[indices[i]], adjustcolor(colours[1], alpha = 0.2), adjustcolor(colours[2], alpha = 0.2)))
  }
}
lines(times, mean_heavy, lwd = 3, col = colours[1])
lines(times, mean_light, lwd = 3, col = colours[2])
legend("topright", inset = 0.02, legend = c(paste0(">60kg (n = ", sum(weight[indices]), ")"), paste0("<60kg (n = ", sum(sum(!weight[indices])), ")")), 
       title = "Age", col = c(colours[1], colours[2]), lty = 1, lwd = 5, cex = 1, title.adj = 0.5)

# Compare these results with old results
# old <- readRDS("Outputs/MCMC_Output_26th_December.rds")
# old_bio_vector <- vector(length = 55)
# new_bio_vector <- vector(length = 55)
# old_kalbso_vector <- vector(length = 55)
# new_kalbso_vector <- vector(length = 55)
# for (i in 1:55) {
#   
#   temp_new <- readRDS(paste0("Outputs/SingleDose_SimpleModel_MCMC_Outputs/TS", i, "_prior", prior_string, "_singleDose_simpleModel_.rds"))
#   new_median_bioavailability <- median(temp_new$mcmc_output$MCMC_Output[max(temp_new$number_iterations/2, start_covariance_adaptation):temp_new$number_iterations, "bioavailability"]) / reparam[2]  
#   new_median_k_alb_so <- median(temp_new$mcmc_output$MCMC_Output[max(temp_new$number_iterations/2, start_covariance_adaptation):temp_new$number_iterations, "k_alb_so"]) / reparam[5]
#   
#   temp_old <- old[[i]]
#   old_median_bioavailability <- median(temp_old[25000:40001, "bioavailability"])   
#   old_median_k_alb_so <- median(temp_old[25000:40001, "k_alb_so"])
#   
#   old_bio_vector[i] <- old_median_bioavailability
#   new_bio_vector[i] <- new_median_bioavailability
#   old_kalbso_vector[i] <- old_median_k_alb_so
#   new_kalbso_vector[i] <- new_median_k_alb_so
#   
#   print(i)
#   
# }
