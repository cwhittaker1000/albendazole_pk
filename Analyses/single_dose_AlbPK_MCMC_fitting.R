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

# tst <- overall %>%
#   group_by(Temporal_ID) %>%
#   filter(row_number()==1) %>%
#   select(Temporal_ID, Dose_Amount)
# for (k in 1:max(overall$Temporal_ID)) {
#   time_series_number <- k
#   Single_PK_Dataset <- filter(overall, Temporal_ID == k) # filter by time_series_number
#   Albendazole_Sulfoxide_Time <- Single_PK_Dataset$Time[Single_PK_Dataset$Metabolite == "AlbSO"]  
#   Albendazole_Sulfoxide_Conc <- Single_PK_Dataset$Converted_Concentration[Single_PK_Dataset$Metabolite == "AlbSO"]  
#   Alb_SO_data <- data.frame(Time = Albendazole_Sulfoxide_Time, Alb_SO = Albendazole_Sulfoxide_Conc)
#   plot(Alb_SO_data$Time, Alb_SO_data$Alb_SO, pch = 20, col = "#7609BA", cex = 2, main = paste0("Time Series ", k))
#   browser()
# }

# Defining the MCMC Inputs
number_iterations <- 25000
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
                      dosing = character(num_TS), number = rep(NA, num_TS), dose_data_availability = rep(NA, num_TS), stringsAsFactors = FALSE)

# Looping Over the Dataset and Fitting Each of the Datasets
fresh_run <- FALSE
for (k in 1:num_TS) {
  
  # Set Seed
  set.seed(140194)
  
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
    if (k == 96 | k == 98) {
      Single_PK_Dataset <- Single_PK_Dataset %>% # for these two, ambiguous whether measured after giving drug - fairly certain measured after 2nd dose so remove 
        filter(Time < time_2nd_dose) 
    } else {
      Single_PK_Dataset <- Single_PK_Dataset %>%
        filter(Time <= time_2nd_dose)  
    }
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
  if (informative_prior) {
    prior_string <- "Inf"
  } else {
    prior_string <- "Uninf"
  }
  if (fresh_run) {
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
                                          output = TRUE, 
                                          informative_prior = informative_prior, 
                                          reparam = reparam, 
                                          refresh = 1000,
                                          time_increment = time_increment)
    
    # Saving Output
    saveRDS(list(number_iterations = number_iterations, info = Single_PK_Dataset, mcmc_output = run_MCMC, Alb_data = Alb_data, Alb_SO_data = Alb_SO_data), 
            paste0("Outputs/SingleDose_AlbPK_ModelFitting_Outputs/TS", k, "_", prior_string, "Prior", "_", dosing, "Dose.rds"))
    # saveRDS(list(number_iterations = number_iterations, info = Single_PK_Dataset, mcmc_output = run_MCMC, Alb_data = Alb_data, Alb_SO_data = Alb_SO_data), 
    #         paste0("Outputs/SingleDose_AlbPK_ModelFitting_Outputs/TS", k, "_", "prior", prior_string, "_singleDose_simpleModel_.rds"))
  } else {
    run_MCMC <- readRDS(paste0("Outputs/SingleDose_AlbPK_ModelFitting_Outputs/TS", k, "_", prior_string, "Prior", "_", dosing, "Dose.rds"))
    # run_MCMC <- readRDS(paste0("Outputs/SingleDose_AlbPK_ModelFitting_Outputs/TS", k, "_", "prior", prior_string, "_singleDose_simpleModel_.rds"))
    run_MCMC <- run_MCMC$mcmc_output
  }
  
  # Plotting Output
  alb_so <- run_MCMC$Alb_SO[burnin:(dim(run_MCMC$MCMC_Output)[1] - 1), ]
  mean <- apply(alb_so, 2, mean)
  lower <- apply(alb_so, 2, quantile, 0.025)
  upper <- apply(alb_so, 2, quantile, 0.975)
  plot(run_MCMC$Times, mean, type = "l", col = "#7609BA", ylab = "Concentration (ng/ml)", xlab = "Time (Hours)", las = 1, lwd = 2, main = paste0("Time Series ", k), ylim = c(0, max(c(upper, Alb_SO_data$Alb_SO))))
  points(Alb_SO_data$Time, Alb_SO_data$Alb_SO, pch = 20, col = "#7609BA", cex = 2)
  polygon(c(run_MCMC$Times, rev(run_MCMC$Times)), c(lower, rev(upper)), col = adjustcolor("#7609BA", alpha.f = 0.2), border = NA)
  
  if("Alb" %in% Single_PK_Dataset$Metabolite) {
    alb <- run_MCMC$Alb[burnin:(dim(run_MCMC$MCMC_Output)[1] - 1), ]
    mean_alb <- apply(alb, 2, mean)
    lower_alb <- apply(alb, 2, quantile, 0.025)
    upper_alb <- apply(alb, 2, quantile, 0.975)
    lines(run_MCMC$Times, mean_alb, type = "l", col = "#E5005F", ylab = "Concentration (ng/ml)", xlab = "Time (Hours)", las = 1, lwd = 2)
    points(Alb_data$Time, Alb_data$Alb, pch = 20, col = "#E5005F", cex = 2)
    polygon(c(run_MCMC$Times, rev(run_MCMC$Times)), c(lower_alb, rev(upper_alb)), col = adjustcolor("#E5005F", alpha.f = 0.2), border = NA)
  }
  
  # Generating and Saving Summary Statistics for Each Fitted Output 
  times <- seq(0, 100, 0.01)
  dose <- dose_info$amount[1]
  
  ## Extracting and Running Model With Median Outputs
  median_k_abs <- median(run_MCMC$MCMC_Output[burnin:(dim(run_MCMC$MCMC_Output)[1] - 1), "k_abs"]) / reparam[1]
  median_bioavailability <- median(run_MCMC$MCMC_Output[burnin:(dim(run_MCMC$MCMC_Output)[1] - 1), "bioavailability"]) / reparam[2]
  median_sigma <- median(run_MCMC$MCMC_Output[burnin:(dim(run_MCMC$MCMC_Output)[1] - 1), "sigma"]) / reparam[3]
  median_k_alb <- median(run_MCMC$MCMC_Output[burnin:(dim(run_MCMC$MCMC_Output)[1] - 1), "k_alb"]) / reparam[4]
  median_k_alb_so <- median(run_MCMC$MCMC_Output[burnin:(dim(run_MCMC$MCMC_Output)[1] - 1), "k_alb_so"]) / reparam[5]
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
  mean_k_abs <- mean(run_MCMC$MCMC_Output[burnin:(dim(run_MCMC$MCMC_Output)[1] - 1), "k_abs"]) / reparam[1]
  mean_bioavailability <- mean(run_MCMC$MCMC_Output[burnin:(dim(run_MCMC$MCMC_Output)[1] - 1), "bioavailability"]) / reparam[2]
  mean_sigma <- mean(run_MCMC$MCMC_Output[burnin:(dim(run_MCMC$MCMC_Output)[1] - 1), "sigma"]) / reparam[3]
  mean_k_alb <- mean(run_MCMC$MCMC_Output[burnin:(dim(run_MCMC$MCMC_Output)[1] - 1), "k_alb"]) / reparam[4]
  mean_k_alb_so <- mean(run_MCMC$MCMC_Output[burnin:(dim(run_MCMC$MCMC_Output)[1] - 1), "k_alb_so"]) / reparam[5]
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
saveRDS(results, paste0("Outputs/SingleDose_AlbPK_ModelFitting_Outputs/AllTS_Results.rds"))


# Plotting all the raw model fitting results
m <- matrix(c(1, 1, 2, 3, 1, 1, 4, 5, 1, 1, 6, 7), nrow = 3, ncol = 4, byrow = TRUE)
prior_string <- "Uninf"
layout(m)
for (i in 1:max(overall$Temporal_ID)) {
  
  # Loading in results
  Single_PK_Dataset <- filter(overall, Temporal_ID == i) # filter by time_series_number
  dosing <- unique(Single_PK_Dataset$Dosing)
  
  temp <- readRDS(paste0("Outputs/SingleDose_AlbPK_ModelFitting_Outputs/TS", i, "_", prior_string, "Prior", "_", dosing, "Dose.rds"))

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

# Loading in MCMC results table 
results <- readRDS("outputs/SingleDose_AlbPK_ModelFitting_Outputs/AllTS_Results.rds")
results$id <- seq(1:dim(results)[1])
results <- results %>% 
  relocate(id, .before = median_k_alb_so) %>%
  select(id, median_k_alb_so, median_bioavailability, median_AUC, median_Cmax)

# Loading in the metadata about each time-series
study_metadata <- read.csv("data/metadata_for_import.csv") %>%
  select(TS_Num, Dosing_Regimen, Data_Availability, 
         Single_Dose_Mg, Dose_Per_Kg, Number_People, 
         Age_Used, Age_Group, Weight_Used,
         Sex, Sex_Ratio, State, 
         Disease_Status, Disease_Specific, Onchocerciasis, Echinococcosis, Neurocysticercosis,
         Co_Drugs_Binary, Ivermectin, DEC, Anti_Epileptics) 

overall_results <- results %>%
  left_join(study_metadata, by = c("id" = "TS_Num")) %>%
  mutate(Sex = ifelse(is.na(Sex), "Unclear", Sex)) %>%
  mutate(State = ifelse(is.na(State), "Unclear", State))
miss_var_summary(overall_results)

# Sex
table(overall_results$Sex)
summary(lm(median_bioavailability ~ Sex, data = overall_results))
summary(lm(median_k_alb_so ~ Sex, data = overall_results))
summary(lm(median_Cmax ~ Sex, data = overall_results))
summary(lm(median_AUC ~ Sex, data = overall_results))
summary(lm(median_bioavailability ~ Sex_Ratio, data = overall_results))
summary(lm(median_k_alb_so ~ Sex_Ratio, data = overall_results))
summary(lm(median_Cmax ~ Sex_Ratio, data = overall_results))
summary(lm(median_AUC ~ Sex_Ratio, data = overall_results))

# Fasting State
table(overall_results$State)
summary(lm(median_bioavailability ~ State, data = overall_results))
summary(lm(median_k_alb_so ~ State, data = overall_results))
summary(lm(median_Cmax ~ State, data = overall_results))
summary(lm(median_AUC ~ State, data = overall_results))

# Dose
table(overall_results$Dose_Per_Kg)
table(overall_results$Single_Dose_Mg)
summary(lm(median_bioavailability ~ Dose_Per_Kg, data = overall_results))
summary(lm(median_k_alb_so ~ Dose_Per_Kg, data = overall_results))
summary(lm(median_Cmax ~ Dose_Per_Kg, data = overall_results))
summary(lm(median_AUC ~ Dose_Per_Kg, data = overall_results))
summary(lm(median_bioavailability ~ Single_Dose_Mg, data = overall_results))
summary(lm(median_k_alb_so ~ Single_Dose_Mg, data = overall_results))
summary(lm(median_Cmax ~ Single_Dose_Mg, data = overall_results))
summary(lm(median_AUC ~ Single_Dose_Mg, data = overall_results))

# Weight
table(overall_results$Weight_Used)
summary(lm(median_bioavailability ~ Weight_Used, data = overall_results))
summary(lm(median_k_alb_so ~ Weight_Used, data = overall_results))
summary(lm(median_Cmax ~ Weight_Used, data = overall_results))
summary(lm(median_AUC ~ Weight_Used, data = overall_results))

# Age
table(overall_results$Age_Used)
table(overall_results$Age_Group)
summary(lm(median_bioavailability ~ Age_Used, data = overall_results))
summary(lm(median_k_alb_so ~ Age_Used, data = overall_results))
summary(lm(median_Cmax ~ Age_Used, data = overall_results))
summary(lm(median_AUC ~ Age_Used, data = overall_results))
summary(lm(median_bioavailability ~ Age_Group, data = overall_results))
summary(lm(median_k_alb_so ~ Age_Group, data = overall_results))
summary(lm(median_Cmax ~ Age_Group, data = overall_results))
summary(lm(median_AUC ~ Age_Group, data = overall_results))

# Infection
table(overall_results$Disease_Status)
table(overall_results$Onchocerciasis)
table(overall_results$Echinococcosis)
table(overall_results$Neurocysticercosis)

summary(lm(median_bioavailability ~ Disease_Status + Dose_Per_Kg, data = overall_results))
summary(lm(median_k_alb_so ~ Disease_Status + Dose_Per_Kg, data = overall_results))
summary(lm(median_Cmax ~ Disease_Status + Dose_Per_Kg, data = overall_results))
summary(lm(median_AUC ~ Disease_Status + Dose_Per_Kg, data = overall_results))

summary(lm(median_bioavailability ~ Onchocerciasis + Dose_Per_Kg, data = overall_results))
summary(lm(median_k_alb_so ~ Onchocerciasis + Dose_Per_Kg, data = overall_results))
summary(lm(median_Cmax ~ Onchocerciasis + Dose_Per_Kg, data = overall_results))
summary(lm(median_AUC ~ Onchocerciasis + Dose_Per_Kg, data = overall_results))

summary(lm(median_bioavailability ~ Echinococcosis + Dose_Per_Kg, data = overall_results))
summary(lm(median_k_alb_so ~ Echinococcosis + Dose_Per_Kg, data = overall_results))
summary(lm(median_Cmax ~ Echinococcosis + Dose_Per_Kg, data = overall_results))
summary(lm(median_AUC ~ Echinococcosis + Dose_Per_Kg, data = overall_results))

summary(lm(median_bioavailability ~ Neurocysticercosis + Dose_Per_Kg, data = overall_results))
summary(lm(median_k_alb_so ~ Neurocysticercosis + Dose_Per_Kg, data = overall_results))
summary(lm(median_Cmax ~ Neurocysticercosis + Dose_Per_Kg, data = overall_results))
summary(lm(median_AUC ~ Neurocysticercosis + Dose_Per_Kg, data = overall_results))

# Drug
table(overall_results$Co_Drugs_Binary)
table(overall_results$DEC)
table(overall_results$Anti_Epileptics)
table(overall_results$Ivermectin)

summary(lm(median_bioavailability ~ Co_Drugs_Binary + Single_Dose_Mg, data = overall_results))
summary(lm(median_k_alb_so ~ Co_Drugs_Binary + Single_Dose_Mg, data = overall_results))
summary(lm(median_Cmax ~ Co_Drugs_Binary + Single_Dose_Mg, data = overall_results))
summary(lm(median_AUC ~ Co_Drugs_Binary + Single_Dose_Mg, data = overall_results))

summary(lm(median_bioavailability ~ DEC + Single_Dose_Mg, data = overall_results))
summary(lm(median_k_alb_so ~ DEC + Single_Dose_Mg, data = overall_results))
summary(lm(median_Cmax ~ DEC + Single_Dose_Mg, data = overall_results))
summary(lm(median_AUC ~ DEC + Single_Dose_Mg, data = overall_results))

summary(lm(median_bioavailability ~ Anti_Epileptics + Single_Dose_Mg, data = overall_results))
summary(lm(median_k_alb_so ~ Anti_Epileptics + Single_Dose_Mg, data = overall_results))
summary(lm(median_Cmax ~ Anti_Epileptics + Single_Dose_Mg, data = overall_results))
summary(lm(median_AUC ~ Anti_Epileptics + Single_Dose_Mg, data = overall_results))

summary(lm(median_bioavailability ~ Ivermectin + Single_Dose_Mg, data = overall_results))
summary(lm(median_k_alb_so ~ Ivermectin + Single_Dose_Mg, data = overall_results))
summary(lm(median_Cmax ~ Ivermectin + Single_Dose_Mg, data = overall_results))
summary(lm(median_AUC ~ Ivermectin + Single_Dose_Mg, data = overall_results))

# Overall
miss_var_summary(overall_results)

summary(lm(median_bioavailability ~ Sex + State + Age_Group + Single_Dose_Mg + Co_Drugs_Binary + Disease_Status, data = overall_results))
summary(lm(median_k_alb_so ~ Sex + State + Age_Group + Single_Dose_Mg + Co_Drugs_Binary + Disease_Status, data = overall_results))
summary(lm(median_Cmax ~ Sex + State + Age_Group + Single_Dose_Mg + Co_Drugs_Binary + Disease_Status, data = overall_results))
summary(lm(median_AUC ~ Sex + State + Age_Group + Single_Dose_Mg + Co_Drugs_Binary + Disease_Status, data = overall_results))

summary(lm(median_bioavailability ~ Sex + State + Age_Group + Dose_Per_Kg + Co_Drugs_Binary + Disease_Status, data = overall_results))
summary(lm(median_k_alb_so ~ Sex + State + Age_Group + Dose_Per_Kg + Co_Drugs_Binary + Disease_Status, data = overall_results))
summary(lm(median_Cmax ~ Sex + State + Age_Group + Dose_Per_Kg + Co_Drugs_Binary + Disease_Status, data = overall_results))
summary(lm(median_AUC ~ Sex + State + Age_Group + Dose_Per_Kg + Co_Drugs_Binary + Disease_Status, data = overall_results))

# Generating median outputs for plotting
times <- seq(0, 50, 0.1)
median_output <- matrix(nrow = max(overall$Temporal_ID), ncol = length(times))
for (i in 1:max(overall$Temporal_ID)) {
  
  # Loading in results
  Single_PK_Dataset <- filter(overall_results, id == i) # filter by time_series_number
  dosing <- unique(Single_PK_Dataset$Dosing)
  dosing <- ifelse(dosing == "Single_Dose", "Single", "Multiple")
  if (i %in% c(95, 96, 97, 98, 99, 100, 107, 108)) {
    dosing <- "Single"
  }
  temp <- readRDS(paste0("Outputs/SingleDose_AlbPK_ModelFitting_Outputs/TS", i, "_", prior_string, "Prior", "_", dosing, "Dose.rds"))
  
  # Running model with median parameter estimates
  median_k_abs <- median(temp$mcmc_output$MCMC_Output[max(temp$number_iterations/2, start_covariance_adaptation):temp$number_iterations, "k_abs"]) / reparam[1]
  median_bioavailability <- median(temp$mcmc_output$MCMC_Output[max(temp$number_iterations/2, start_covariance_adaptation):temp$number_iterations, "bioavailability"]) / reparam[2]
  median_sigma <- median(temp$mcmc_output$MCMC_Output[max(temp$number_iterations/2, start_covariance_adaptation):temp$number_iterations, "sigma"]) / reparam[3]
  median_k_alb <- median(temp$mcmc_output$MCMC_Output[max(temp$number_iterations/2, start_covariance_adaptation):temp$number_iterations, "k_alb"]) / reparam[4]
  median_k_alb_so <- median(temp$mcmc_output$MCMC_Output[max(temp$number_iterations/2, start_covariance_adaptation):temp$number_iterations, "k_alb_so"]) / reparam[5]
  median_model <- Albendazole_PK_Model(k_abs = median_k_abs, sigma = median_sigma, 
                                       k_alb_so = median_k_alb_so, k_alb = median_k_alb, 
                                       gut_1 = (median_bioavailability/1e+5) * (400 * 1e+6), 
                                       gut_2 = 0, liver = 0, blood_alb = 0, blood_alb_so = 0)
  median_out <- median_model$run(times)
  median_alb_so_PK <- median_out[, 6]
  median_output[i, ] <- median_alb_so_PK
  
  print(i)
  
}

# Plotting dosing results
# par(mfrow = c(1, 1))
# dosing <- factor(results$dosing)
# colours <- c("#4E90CE", "#CE4E61")
# indices <- which(!is.na(dosing))
# mean_single <- apply(median_output[dosing == "Single", ], 2, mean)
# mean_multiple <- apply(median_output[dosing == "Multiple", ], 2, mean)
# for (i in indices) {
#   if (i == 1) {
#     plot(times, median_output[i, ], type = "l", ylim = c(0, 1000), las = 1, xlab = "Time Since Treatment (Hours)", ylab = "Concentration (ng/ml)",
#          col = ifelse(dosing[i] == "Single", adjustcolor(colours[1], alpha = 0.2), adjustcolor(colours[2], alpha = 0.2)))
#   } else {
#     lines(times, median_output[i, ], col = ifelse(dosing[i] == "Single", adjustcolor(colours[1], alpha = 0.2), adjustcolor(colours[2], alpha = 0.2)))
#   }
# }
# lines(times, mean_single, lwd = 3, col = colours[1])
# lines(times, mean_multiple, lwd = 3, col = colours[2])
# legend("topright", inset = 0.02, legend = c(paste0("Single (n = ", sum(dosing == "Single"), ")"), 
#                                             paste0("Multiple (n = ", sum(dosing == "Multiple"), ")")), 
#        title = "Dosing", col = c(colours[1], colours[2]), lty = 1, lwd = 5, cex = 1, title.adj = 0.5)

# Plotting sex results
par(mfrow = c(2, 4))
sex <- factor(overall_results$Sex)
colours <- c("#4E90CE", "#CE4E61")
indices <- which(!is.na(sex) & sex != "Unclear" & sex != "FeMale")
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
state <- factor(overall_results$State)
colours <- c("#4BBC00","#F28123")
indices <- which(!is.na(state) & state != "Unclear" & state != "Mixture")
mean_fatty_meal <- apply(median_output[state == "Fatty Meal", ], 2, mean)
mean_fasted <- apply(median_output[state == "Fasted", ], 2, mean)
for (i in indices) {
  if (i == 1) {
    plot(times, median_output[i, ], type = "l", ylim = c(0, 1000), las = 1, xlab = "Time Since Treatment (Hours)", ylab = "Concentration (ng/ml)",
         col = ifelse(state[i] == "Fatty Meal", adjustcolor(colours[1], alpha = 0.2), adjustcolor(colours[2], alpha = 0.2)))
  } else {
    lines(times, median_output[i, ], col = ifelse(state[i] == "Fatty Meal", adjustcolor(colours[1], alpha = 0.2), adjustcolor(colours[2], alpha = 0.2)))
  }
}
lines(times, mean_fatty_meal, lwd = 3, col = colours[1])
lines(times, mean_fasted, lwd = 3, col = colours[2])
legend("topright", inset = 0.02, legend = c(paste0("Fatty Meal (n = ", sum(state == "Fatty Meal"), ")"), paste0("Fasted (n = ", sum(state == "Fasted"), ")")), 
       title = "State", col = c(colours[1], colours[2]), lty = 1, lwd = 5, cex = 1, title.adj = 0.5)

# Plotting dose results
overall_results$Dose_Binary <- ifelse(overall_results$Single_Dose_Mg <= 400, "Low", "High")
dose <- factor(overall_results$Dose_Binary)
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
infection <- factor(overall_results$Disease_Status)
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
drug <- factor(overall_results$Co_Drugs_Binary)
colours <- c("#ADA9A9", "#D62C2C")
indices <- which(!is.na(drug))
mean_none <- apply(median_output[drug == "No", ], 2, mean)
mean_drug <- apply(median_output[drug == "Yes", ], 2, mean)
for (i in indices) {
  if (i == 1) {
    plot(times, median_output[i, ], type = "l", ylim = c(0, 1000), las = 1, xlab = "Time Since Treatment (Hours)", ylab = "Concentration (ng/ml)",
         col = ifelse(drug[i] == "None", adjustcolor(colours[1], alpha = 0.2), adjustcolor(colours[2], alpha = 0.2)))
  } else {
    lines(times, median_output[i, ], col = ifelse(drug[i] == "No", adjustcolor(colours[1], alpha = 0.2), adjustcolor(colours[2], alpha = 0.2)))
  }
}
lines(times, mean_none, lwd = 3, col = colours[1])
lines(times, mean_drug, lwd = 3, col = colours[2])
legend("topright", inset = 0.02, legend = c(paste0("None (n = ", sum(drug == "No"), ")"), paste0("Drug (n = ", sum(drug == "Yes"), ")")), 
       title = "Drug Co-Administration", col = c(colours[1], colours[2]), lty = 1, lwd = 5, cex = 1, title.adj = 0.5)

# Plotting age results
age <- factor(overall_results$Age_Group)
colours <- c("#E0A02A", "#ED7BEB")
indices <- which(!is.na(age) & age != "Unclear" & age != "Mixture")
mean_adults <- apply(median_output[!is.na(age) & age == "Adults", ], 2, mean)
mean_children <- apply(median_output[!is.na(age) & age == "Children", ], 2, mean)
for (i in indices) {
  if (i == 1) {
    plot(times, median_output[i, ], type = "l", ylim = c(0, 1000), las = 1, xlab = "Time Since Treatment (Hours)", ylab = "Concentration (ng/ml)",
         col = ifelse(age[i] == "Adults", adjustcolor(colours[1], alpha = 0.2), adjustcolor(colours[2], alpha = 0.2)))
  } else {
    lines(times, median_output[i, ], col = ifelse(age[i] == "Adults", adjustcolor(colours[1], alpha = 0.2), adjustcolor(colours[2], alpha = 0.2)))
  }
}
lines(times, mean_adults, lwd = 3, col = colours[1])
lines(times, mean_children, lwd = 3, col = colours[2])
legend("topright", inset = 0.02, legend = c(paste0("Adults (n = ", sum(!is.na(age) & age == "Adults"), ")"), paste0("Children (n = ", sum(!is.na(age) & age == "Children"), ")")), 
       title = "Age", col = c(colours[1], colours[2]), lty = 1, lwd = 5, cex = 1, title.adj = 0.5)

# Plotting weight results
weight <- overall_results$Weight_Used > 60
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
       title = "Weight", col = c(colours[1], colours[2]), lty = 1, lwd = 5, cex = 1, title.adj = 0.5)

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
