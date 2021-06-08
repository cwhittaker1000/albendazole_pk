# Loading In Required Libraries 
library(tictoc); library(deSolve); library(MALDIquant); library(mvtnorm); library(MASS); library(tmvtnorm); library(odin); library(dplyr); 

# Generate model instance and source relevant functiosn
source("Models/single_dose_model_ODIN.R")
source("Functions/single_dose_MCMC_Functions.R")
source("Functions/Processing_and_Plotting_Functions.R")

# Importing and Processing the Data
multi_data <- read.csv("Data/Albendazole_Pharmacokinetic_Data/multiple_dosing_Alb_PK_data.csv", stringsAsFactors = FALSE)

# Defining the MCMC Inputs
number_iterations <- 10000
sd_proposals <- c(0.002, 0.002, 0.002, 0.002, 0.002)
reparam <- c(10, 1000, 0.1, 100, 10)
start_covariance_adaptation <- 2500
informative_prior <- FALSE
raw_MCMC_output <- vector(mode = "list", length = length(unique(multi_data$Temporal_ID)))

# Looping Over the Dataset and Fitting Each of the Datasets
par(mfrow = c(5, 5))
par(mar = c(2, 2, 0.5, 1))
par(oma = c(0, 2, 0, 0))
plot_limit <- FALSE
# ts <- c(7, 8, 9, 10, 11, 12, 16, 35, 36, 42, 43, 44, 51, 53, 54)
for (k in 1:max(All_PK_Data$Temporal_ID)) {
  
  initial_parameters <- c(k_abs = 2, bioavailability = 0.001, sigma = 15, k_alb_so = 0.1300, k_alb = 0.2)
  time_series_number <- k
  Single_PK_Dataset <- filter(All_PK_Data, Temporal_ID == k) # filter by time_series_number
  dose <- Single_PK_Dataset$Dose_mg[1] 
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
  
  bloop <- MCMC_running(number_iterations, initial_parameters, sd_proposals, start_covariance_adaptation, Alb_data, Alb_SO_data, metabolite_availability, Albendazole_PK_Model, end_time, dose, FALSE, informative_prior, reparam)
  # par(mfrow = c(3, 2))
  # for (j in 1:ncol(bloop$MCMC_Output)) {
  #   cl <- rainbow(ncol(bloop$MCMC_Output))
  #   plot(bloop$MCMC_Output[1:number_iterations, j], type = "l", lwd = 2, col = cl[j])
  # }

  burn_in <- 0.5
  iterations <- number_iterations
  chain <- bloop$MCMC_Output[(burn_in*number_iterations):number_iterations, ]
  raw_MCMC_output[[k]] <- chain
  
  median_k_abs <- median(chain[, 1])/reparam[1]
  median_bioavailability <- median(chain[, 2])/reparam[2]
  median_sigma <- median(chain[, 3])/reparam[3]
  median_k_alb_so <- median(chain[, 4])/reparam[4]
  median_k_alb <- median(chain[, 5])/reparam[5]
  model_runner <- Albendazole_PK_Model(k_abs = median_k_abs, bioavailability = median_bioavailability, sigma = median_sigma, k_alb_so = median_k_alb_so, k_alb = median_k_alb, dose = dose) 
  times <- seq(0, 100, length.out = 1000)
  median_out <- model_runner$run(times)
  max <- max(c(as.vector(median_out[, 5:6]), Alb_SO_data$Alb_SO))
  
  # Generating the 95% Credible Intervals
  size_of_burned_chain <- number_iterations - (burn_in * number_iterations)
  thinned_chain <- chain # [seq(1, size_of_burned_chain, 10), ]
  Alb_SO_Storage_matrix <- matrix(nrow = length(thinned_chain[, 1]), ncol = length(median_out[, 1]))
  Alb_Storage_matrix <- matrix(nrow = length(thinned_chain[, 1]), ncol = length(median_out[, 1]))
  for (i in 1:length(thinned_chain[, 1])) {
    k_abs <- thinned_chain[i, 1]/reparam[1]
    bioavailability <- thinned_chain[i, 2]/reparam[2]
    sigma <- thinned_chain[i, 3]/reparam[3]
    k_alb_so <- thinned_chain[i, 4]/reparam[4]
    k_alb <- thinned_chain[i, 5]/reparam[5]
    model_runner <- Albendazole_PK_Model(k_abs = k_abs, bioavailability = bioavailability, sigma = sigma, k_alb_so = k_alb_so, dose = dose, k_alb = k_alb)
    out <- model_runner$run(times)
    Alb_Storage_matrix[i, ] <- out[, 5]
    Alb_SO_Storage_matrix[i, ] <- out[, 6]
  }
  Alb_Credible_Upper <- apply(Alb_Storage_matrix, 2, quantile, 0.975)
  Alb_Credible_Lower <- apply(Alb_Storage_matrix, 2, quantile, 0.025)
  Alb_SO_Credible_Upper <- apply(Alb_SO_Storage_matrix, 2, quantile, 0.975)
  Alb_SO_Credible_Lower <- apply(Alb_SO_Storage_matrix, 2, quantile, 0.025)

  # Plotting the Median Output Along With the 95% Credible Intervals 
  # par(mfrow = c(1, 1))
  if (plot_limit == TRUE) {
    plot(median_out[, 1], median_out[, 5], type = "l", col = "#E5005F",  ylim = c(0, max), xlim = c(0, 48), ylab = "Concentration (ng/ml)", xlab = "Time (Hours)", las = 1, lwd = 2)
    lines(median_out[, 1], median_out[, 6], type = "l", col = "#7609BA", lwd = 2)
    points(Alb_data$Time, Alb_data$Alb, pch = 20, col = "#E5005F", cex = 2)
    points(Alb_SO_data$Time, Alb_SO_data$Alb_SO, pch = 20, col = "#7609BA", cex = 2)
    polygon(c(median_out[, 1], rev(median_out[, 1])), c(Alb_Credible_Lower, rev(Alb_Credible_Upper)), col = adjustcolor("#E5005F", alpha.f = 0.2), border = NA)
    polygon(c(median_out[, 1], rev(median_out[, 1])), c(Alb_SO_Credible_Lower, rev(Alb_SO_Credible_Upper)), col = adjustcolor("#7609BA", alpha.f = 0.2), border = NA)
  } else {
    plot(median_out[, 1], median_out[, 5], type = "l", col = "#E5005F",  ylim = c(0, max), xlim = c(0, end_time), ylab = "Concentration (ng/ml)", xlab = "Time (Hours)", las = 1, lwd = 2)
    lines(median_out[, 1], median_out[, 6], type = "l", col = "#7609BA", lwd = 2)
    points(Alb_data$Time, Alb_data$Alb, pch = 20, col = "#E5005F", cex = 2)
    points(Alb_SO_data$Time, Alb_SO_data$Alb_SO, pch = 20, col = "#7609BA", cex = 2)
    polygon(c(median_out[, 1], rev(median_out[, 1])), c(Alb_Credible_Lower, rev(Alb_Credible_Upper)), col = adjustcolor("#E5005F", alpha.f = 0.2), border = NA)
    polygon(c(median_out[, 1], rev(median_out[, 1])), c(Alb_SO_Credible_Lower, rev(Alb_SO_Credible_Upper)), col = adjustcolor("#7609BA", alpha.f = 0.2), border = NA)
  }
  print(k)
}

saveRDS(raw_MCMC_output, file = "C:/Users/cw1716/Documents/Q_Drive_Copy/Loa_Loa/Drug Efficacy Modelling/Albendazole Modelling/MCMC_Output_26th_December.rds")
raw_MCMC_output <- readRDS("C:/Users/cw1716/Documents/Q_Drive_Copy/Loa_Loa/Drug Efficacy Modelling/Albendazole Modelling/MCMC_Output_26th_December.rds")

par(mfrow = c(5, 11))
for (k in 1:max(All_PK_Data$Temporal_ID)) {
  MCMC_output_plotting_function(raw_MCMC_output, k, plot_limit = TRUE, plot_separately = FALSE, plot_both = TRUE, credible_intervals = TRUE)  
}


# Plotting All the Alb-SO Lines On the Same Graph
par(mfrow = c(2, 4))
par(mar = c(4, 5.2, 3, 3))
sex <- time_series_information$Sex_Binary
state <- time_series_information$Feeding_State
dose_info <- time_series_information$Dose_Binary
infected <- time_series_information$Disease_Binary
drugs <- time_series_information$Drug_Binary
age_group <- ifelse(is.na(time_series_information$Age_Cont), "Unsure", ifelse(time_series_information$Age_Cont <= 18, "Young", "Old"))
weight_group <- ifelse(is.na(time_series_information$Weight), "Unsure", ifelse(time_series_information$Weight <= 50, "Light", "Heavy"))
plotting_function(All_PK_Data, raw_MCMC_output, "Sex", colours = c("#4E90CE", "#CE4E61"))
plotting_function(All_PK_Data, raw_MCMC_output, "State", colours = c("#F28123", "#4BBC00"))
plotting_function(All_PK_Data, raw_MCMC_output, "Dose", colours = c("#3ED6BA", "#725FAD"))
plotting_function(All_PK_Data, raw_MCMC_output, "Infection", colours = c("#F7C911", "#CE5A90"))
plotting_function(All_PK_Data, raw_MCMC_output, "Drug", colours = c("#ADA9A9", "#D62C2C"))
plotting_function(All_PK_Data, raw_MCMC_output, "Age", colours = c("#ED7BEB", "#E0A02A"))
plotting_function(All_PK_Data, raw_MCMC_output, "Weight", colours = c("#32720C", "#CCB5FF"))


parameters <- matrix(nrow = 55, ncol = 5)
colnames(parameters) <- c("k_abs", "bioavailability", "sigma", "k_alb_so", "k_alb")
bioavailability_vector <- c() 
k_alb_so_vector <- c()
AUC_vector <- c()
C_max_vector <- c()
times <- seq(0, 100, length.out = 400)
model_output <- matrix(nrow = 55, ncol = 400)

for (i in 1:55) {
  temp <- raw_MCMC_output[[i]]
  temp_k_abs <- median(temp[20000:40000, 1])
  temp_bioavailability <- median(temp[20000:40000, 2])
  temp_sigma <- median(temp[20000:40000, 3])
  temp_k_alb_so <- median(temp[20000:40000, 4])
  temp_k_alb <- median(temp[20000:40000, 5])
  
  bioavailability_vector[i] <- temp_bioavailability
  k_alb_so_vector[i] <- temp_k_alb_so
  parameters[i, ] <- c(temp_k_abs, temp_bioavailability, temp_sigma, temp_k_alb_so, temp_k_alb)
  dose <- time_series_information$Dose_mg[i]
  
  model_runner <- Albendazole_PK_Model(k_abs = temp_k_abs, bioavailability = temp_bioavailability, sigma = temp_sigma, k_alb_so = temp_k_alb_so, dose = dose, k_alb = temp_k_alb) 
  median_out <- model_runner$run(times)
  model_output[i, ] <- median_out[, 6]
  
  AUC_vector[i] <- DescTools::AUC(times, model_output[i, ])
  C_max_vector[i] <- max(model_output[i, ])
  print(i)
}

sex_ratio <- as.numeric(time_series_information$Sex_Ratio)
state <- as.character(time_series_information$Feeding_State)
state <- replace(state, state == "Mixture", NA)
state <- replace(state, state == "Unclear", NA)
dose_info <- time_series_information$Dose_mg
infected <- as.character(time_series_information$Disease_Binary)
infected <- replace(infected, infected == "Mixture", NA)
drugs <- as.character(time_series_information$Drug_Binary)
age <- time_series_information$Age_Cont
weight <- time_series_information$Weight
lm_compatible <- data.frame(bioavailability = parameters[, "bioavailability"],
                            half_life = 1/parameters[, "k_alb_so"],
                            AUC = AUC_vector,
                            Cmax = C_max_vector,
                            sex = sex_ratio,
                            state = state, 
                            dose = dose_info, 
                            drugs = drugs,
                            age = age,
                            weight = weight)

cor(parameters[, "bioavailability"], parameters[, "k_alb_so"])
cor(parameters[, "bioavailability"], C_max_vector)
cor(parameters[, "bioavailability"], AUC_vector)
cor(C_max_vector, parameters[, "k_alb_so"])
cor(AUC_vector, parameters[, "k_alb_so"])
cor(AUC_vector, C_max_vector)

table(time_series_information$Sex_Ratio, useNA = "ifany")
table(time_series_information$Feeding_State, useNA = "ifany")
table(time_series_information$Dose_mg, useNA = "ifany")
table(time_series_information$Drug_Binary, useNA = "ifany")
table(time_series_information$Age_Cont, useNA = "ifany")
table(time_series_information$Weight, useNA = "ifany")
table(time_series_information$Disease, useNA = "ifany")

table(time_series_information$Disease_Binary, time_series_information$Sex_Ratio, useNA = "ifany")
table(time_series_information$Disease_Binary, time_series_information$Feeding_State, useNA = "ifany")
table(time_series_information$Disease_Binary, time_series_information$Dose_mg, useNA = "ifany")
table(time_series_information$Disease_Binary, time_series_information$Drug_Binary, useNA = "ifany")
table(time_series_information$Disease_Binary, time_series_information$Age_Cont, useNA = "ifany")
table(time_series_information$Disease_Binary, time_series_information$Weight, useNA = "ifany")

table(time_series_information$Sex_Ratio, time_series_information$Feeding_State, useNA = "ifany")
table(time_series_information$Sex_Ratio, time_series_information$Dose_mg, useNA = "ifany")
table(time_series_information$Sex_Ratio, time_series_information$Drug_Binary, useNA = "ifany")
table(time_series_information$Sex_Ratio, time_series_information$Age_Cont, useNA = "ifany")
table(time_series_information$Sex_Ratio, time_series_information$Weight, useNA = "ifany")

table(time_series_information$Feeding_State, time_series_information$Dose_mg, useNA = "ifany")
table(time_series_information$Feeding_State, time_series_information$Drug_Binary, useNA = "ifany")
table(time_series_information$Feeding_State, time_series_information$Age_Cont, useNA = "ifany")
table(time_series_information$Feeding_State, time_series_information$Weight, useNA = "ifany")

table(time_series_information$Dose_mg, time_series_information$Drug_Binary, useNA = "ifany")
table(time_series_information$Dose_mg, time_series_information$Age_Cont, useNA = "ifany")
table(time_series_information$Dose_mg, time_series_information$Weight, useNA = "ifany")

table(time_series_information$Age_Cont, time_series_information$Weight, useNA = "ifany")

complete_studies <- time_series_information[!is.na(time_series_information$Sex_Ratio) &
                                            (time_series_information$Feeding_State != "Mixture" & time_series_information$Feeding_State != "Unclear") &
                                            !is.na(time_series_information$Age_Cont) &
                                            !is.na(time_series_information$Weight), ]

sex_lost_studies <- time_series_information[!is.na(time_series_information$Sex_Ratio) & 
                                            ((time_series_information$Feeding_State == "Mixture" | time_series_information$Feeding_State == "Unclear") |
                                               is.na(time_series_information$Age_Cont) | 
                                               is.na(time_series_information$Weight)), ]



# Regressions - Bioavailability
bioavailability_sex <- lm(bioavailability ~ sex, data = lm_compatible)
summary(bioavailability_sex)

bioavailability_state <- lm(bioavailability ~ state, data = lm_compatible)
summary(bioavailability_state)

bioavailability_dose <- lm(bioavailability ~ dose, data = lm_compatible)
summary(bioavailability_dose)

bioavailability_drugs <- lm(bioavailability ~ drugs, data = lm_compatible)
summary(bioavailability_drugs)

bioavailability_age <- lm(bioavailability ~ age, data = lm_compatible)
summary(bioavailability_age)

bioavailability_weight <- lm(bioavailability ~ weight, data = lm_compatible)
summary(bioavailability_weight)

bioavailability_coinfection <- lm(bioavailability ~ infected, data = lm_compatible)
summary(bioavailability_coinfection)

bioavailability_all <- lm(bioavailability ~ sex + state + dose + drugs + age + weight + infected, data = lm_compatible)
summary(bioavailability_all)

bloop <- summary(bioavailability_all)
bloop$coefficients
  
# Regressions - Half Life
half_life_sex <- lm(half_life ~ sex, data = lm_compatible)
summary(half_life_sex)

half_life_state <- lm(half_life ~ state, data = lm_compatible)
summary(half_life_state)

half_life_dose <- lm(half_life ~ dose, data = lm_compatible)
summary(half_life_dose)

half_life_drugs <- lm(half_life ~ drugs, data = lm_compatible)
summary(half_life_drugs)

half_life_age <- lm(half_life ~ age, data = lm_compatible)
summary(half_life_age)

half_life_weight <- lm(half_life ~ weight, data = lm_compatible)
summary(half_life_weight)

half_life_coinfection <- lm(half_life ~ infected, data = lm_compatible)
summary(half_life_coinfection)

half_life_all <- lm(half_life ~ sex + state + dose + drugs + age + weight + infected, data = lm_compatible)
summary(half_life_all)

# Regressions - AUC
AUC_sex <- lm(AUC ~ sex, data = lm_compatible)
summary(AUC_sex)

AUC_state <- lm(AUC ~ state, data = lm_compatible)
summary(AUC_state)

AUC_dose <- lm(AUC ~ dose, data = lm_compatible)
summary(AUC_dose)

AUC_drugs <- lm(AUC ~ drugs, data = lm_compatible)
summary(AUC_drugs)

AUC_age <- lm(AUC ~ age, data = lm_compatible)
summary(AUC_age)

AUC_weight <- lm(AUC ~ weight, data = lm_compatible)
summary(AUC_weight)

AUC_coinfection <- lm(AUC ~ infected, data = lm_compatible)
summary(AUC_coinfection)

AUC_all <- lm(AUC ~ sex + state + dose + drugs + age + weight + infected, data = lm_compatible)
summary(AUC_all)

# Univariate Regressions - Cmax
Cmax_sex <- lm(Cmax ~ sex, data = lm_compatible)
summary(Cmax_sex)

Cmax_state <- lm(Cmax ~ state, data = lm_compatible)
summary(Cmax_state)

Cmax_dose <- lm(Cmax ~ dose, data = lm_compatible)
summary(Cmax_dose)

Cmax_drugs <- lm(Cmax ~ drugs, data = lm_compatible)
summary(Cmax_drugs)

Cmax_age <- lm(Cmax ~ age, data = lm_compatible)
summary(Cmax_age)

Cmax_weight <- lm(Cmax ~ weight, data = lm_compatible)
summary(Cmax_weight)

Cmax_coinfection <- lm(Cmax ~ infected, data = lm_compatible)
summary(Cmax_coinfection)

Cmax_all <- lm(Cmax ~ sex + state + dose + drugs + age + weight + infected, data = lm_compatible)
summary(Cmax_all)


mean(time_series_information$Sex_Ratio, na.rm = TRUE)

###____________

All_PK_Data <- read.csv("C:/Users/cw1716/Documents/Q_Drive_Copy/Loa_Loa/Drug Efficacy Modelling/Albendazole Modelling/Albendazole_Pharmacokinetic_Data/Alb_PK_Data.csv")

unique(All_PK_Data$Temporal_ID)

time_series_information <- All_PK_Data[!duplicated(All_PK_Data$Temporal_ID), ]

table(time_series_information$Dose_Binary)
time_series_information$Dose_mg[order(time_series_information$Dose_mg)]
median(time_series_information$Dose_mg[order(time_series_information$Dose_mg)])

table(time_series_information$Sex_Ratio)
time_series_information$Sex_Ratio[order(time_series_information$Sex_Ratio)]
median(time_series_information$Sex_Ratio[order(time_series_information$Sex_Ratio)], na.rm = TRUE)

table(time_series_information$Age_Cont)
time_series_information$Age_Cont[order(time_series_information$Age_Cont)]
median(time_series_information$Age_Cont[order(time_series_information$Age_Cont)], na.rm = TRUE)

table(time_series_information$Weight)
time_series_information$Weight[order(time_series_information$Weight)]
median(time_series_information$Weight[order(time_series_information$Weight)], na.rm = TRUE)

table(time_series_information$Sex_Binary)
table(time_series_information$Feeding_State)
table(time_series_information$Age_Binary)
table(time_series_information$Metabolite)
sum(time_series_information$Number)
table(time_series_information$Disease_Binary)
table(time_series_information$Disease)
table(time_series_information$Drug_Binary)
table(time_series_information$Drug, time_series_information$Drug_Binary)

bloop <- time_series_information[time_series_information$Drug == "None" & time_series_information$Drug_Binary == "Yes", ]


