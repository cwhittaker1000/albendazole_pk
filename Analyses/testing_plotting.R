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

prior_string <- "Uninf"

# Plotting all the raw model fitting results
for (i in 61:max(overall$Temporal_ID)) {
  
  # Loading in results
  Single_PK_Dataset <- filter(overall, Temporal_ID == i) # filter by time_series_number
  dosing <- unique(Single_PK_Dataset$Dosing)
  if (i <= 60) {
    temp <- readRDS(paste0("Outputs/SingleDose_AlbPK_ModelFitting_Outputs/TS", i, "_priorUninf_singleDose_simpleModel_.rds"))
  } else {
    temp <- readRDS(paste0("Outputs/SingleDose_AlbPK_ModelFitting_Outputs/TS", i, "_priorUninf_multipleDose_simpleModel_.rds"))
  }
  alb_so <- temp$mcmc_output$Alb_SO[max(temp$number_iterations/2, start_covariance_adaptation):temp$number_iterations, ]
  mean_alb_so <- apply(alb_so, 2, mean)
  lower_alb_so <- apply(alb_so, 2, quantile, 0.025)
  upper_alb_so <- apply(alb_so, 2, quantile, 0.975)
  
  alb <- temp$mcmc_output$Alb[max(temp$number_iterations/2, start_covariance_adaptation):temp$number_iterations, ]
  mean_alb <- apply(alb, 2, mean)
  lower_alb <- apply(alb, 2, quantile, 0.025)
  upper_alb <- apply(alb, 2, quantile, 0.975)
  
  # Plotting model outputs
  if (i <= 60) {
    xmin <- 0
  } else {
    xmin <- max(0, min(Single_PK_Dataset$Time) - 15) 
  }
  plot(temp$mcmc_output$Times, mean_alb_so, type = "l", col = "#7609BA", ylab = "Concentration (ng/ml)", xlab = "Time (Hours)", las = 1, lwd = 2, 
       main = paste0("Time Series ", i, " - Dose = ", dosing, "; Data = ", unique(Single_PK_Dataset$Data)), 
       ylim = c(0, max(upper_alb_so, temp$info$Converted_Concentration[temp$info$Metabolite == "AlbSO"])),
       xlim = c(xmin, max(Single_PK_Dataset$Time)))
  points(temp$info$Time[temp$info$Metabolite == "AlbSO"], temp$info$Converted_Concentration[temp$info$Metabolite == "AlbSO"], pch = 20, col = "#7609BA", cex = 2)
  polygon(c(temp$mcmc_output$Times, rev(temp$mcmc_output$Times)), c(lower_alb_so, rev(upper_alb_so)), col = adjustcolor("#7609BA", alpha.f = 0.2), border = NA)

  browser()
  
}
