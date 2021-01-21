plotting_function <- function(All_PK_Data, raw_MCMC_output, colouring_variable, mean_plotting, colours) {
  text_size <- 1.5
  if (colouring_variable == "Sex") {
    
    xlab <- "Time Since Treatment (Hours)"
    ylab <- "Concentration (ng/ml)"
    Male_model_output_storage <- matrix(nrow = sum(sex == "Male"), ncol = length(seq(0, 100, length.out = 1000)))
    Other_model_output_storage <- matrix(nrow = sum(sex != "Male" & sex != "Unclear"), ncol = length(seq(0, 100, length.out = 1000)))
    male_counter <- 1
    other_counter <- 1
    
    for (k in 1:max(All_PK_Data$Temporal_ID)) { 
      
      Single_PK_Dataset <- filter(All_PK_Data, Temporal_ID == k) # filter by time_series_number
      dose <- Single_PK_Dataset$Dose_mg[1]  # need to alter the model so that it takes dose as an input, probably have it accessed as part of Single_PK_Dataset within the MCMC
      end_time <- max(Single_PK_Dataset$Time)
      
      Albendazole_Time <- Single_PK_Dataset$Time[Single_PK_Dataset$Metabolite == "Alb"]  
      Albendazole_Conc <- Single_PK_Dataset$Converted_Concentration[Single_PK_Dataset$Metabolite == "Alb"]  
      Alb_data <- data.frame(Time = Albendazole_Time, Alb = Albendazole_Conc)
      
      Albendazole_Sulfoxide_Time <- Single_PK_Dataset$Time[Single_PK_Dataset$Metabolite == "AlbSO"]  
      Albendazole_Sulfoxide_Conc <- Single_PK_Dataset$Converted_Concentration[Single_PK_Dataset$Metabolite == "AlbSO"]  
      Alb_SO_data <- data.frame(Time = Albendazole_Sulfoxide_Time, Alb_SO = Albendazole_Sulfoxide_Conc)
      
      individual_dataset_MCMC_output <- raw_MCMC_output[[k]]
      median_k_abs <- median(individual_dataset_MCMC_output[, 1])
      median_bioavailability <- median(individual_dataset_MCMC_output[, 2])
      median_sigma <- median(individual_dataset_MCMC_output[, 3])
      median_k_alb_so <- median(individual_dataset_MCMC_output[, 4])
      median_k_alb <- median(individual_dataset_MCMC_output[, 5])
      model_runner <- Albendazole_PK_Model(k_abs = median_k_abs, bioavailability = median_bioavailability, sigma = median_sigma, k_alb_so = median_k_alb_so, dose = dose, k_alb = median_k_alb) 
      times <- seq(0, 100, length.out = 1000)
      median_out <- model_runner$run(times)
      max <- max(c(as.vector(median_out[, 5:6]), Alb_SO_data$Alb_SO))
      
      if (sex[k] == "Male") {
        colour <- colours[1]
        Male_model_output_storage[male_counter, ] <- median_out[, 6]
        male_counter <- male_counter + 1
        if (k == 1) {
          plot(median_out[, 1], median_out[, 6], type = "l", col = adjustcolor(colour, alpha.f = 0.2),  ylim = c(0, 1000), xlim = c(0, 48), ylab = "", xlab = xlab, las = 1, lwd = 2, cex.lab = text_size, cex.axis = text_size)
          title(ylab = ylab, line = 4, cex.lab = text_size)
        } else {
          lines(median_out[, 1], median_out[, 6], type = "l", col = adjustcolor(colour, alpha.f = 0.2),  ylim = c(0, max), xlim = c(0, 48), ylab = "", xlab = "", las = 1, lwd = 2)
        }
      } else if (sex[k] == "Mixture") {
        colour <- colours[2]
        Other_model_output_storage[other_counter, ] <- median_out[, 6]
        other_counter <- other_counter + 1
        if (k == 1) {
          plot(median_out[, 1], median_out[, 6], type = "l", col = adjustcolor(colour, alpha.f = 0.2),  ylim = c(0, 1000), xlim = c(0, 48), ylab = "", xlab = xlab, las = 1, lwd = 2, cex.lab = text_size, cex.axis = text_size)
          title(ylab = ylab, line = 4, cex.lab = text_size)
        } else {
          lines(median_out[, 1], median_out[, 6], type = "l", col = adjustcolor(colour, alpha.f = 0.2),  ylim = c(0, max), xlim = c(0, 48), ylab = "", xlab = "", las = 1, lwd = 2)
        }
      }
    }
    Male_mean_model_output <- apply(Male_model_output_storage, 2, mean)
    Other_mean_model_output <- apply(Other_model_output_storage, 2, mean)
    lines(median_out[, 1], Male_mean_model_output, col = colours[1], lwd = 5)
    lines(median_out[, 1], Other_mean_model_output, col = colours[2], lwd = 5)
    legend("topright", inset = 0.02, legend = c(paste0("Males (n = ", male_counter - 1, ")"), paste0("Mixture (n = ", other_counter - 1, ")")), title = "Sex", col = c(colours[1], colours[2]), lty = 1, lwd = 5, cex = 1.25, title.adj = 0.5)
  }
  else if (colouring_variable == "State") {
    
    xlab <- "Time Since Treatment (Hours)"
    ylab <- "Concentration (ng/ml)"
    Fasted_model_output_storage <- matrix(nrow = sum(state == "Fasted"), ncol = length(seq(0, 100, length.out = 1000)))
    Fatty_meal_model_output_storage <- matrix(nrow = sum(state == "Fatty_meal"), ncol = length(seq(0, 100, length.out = 1000)))
    Unclear_model_output_storage <- matrix(nrow = sum(state != "Fasted" & state != "Fatty_meal"), ncol = length(seq(0, 100, length.out = 1000)))
    Fasted_counter <- 1
    Fatty_meal_counter <- 1
    Unclear_counter <- 1
    
    for (k in 1:max(All_PK_Data$Temporal_ID)) { 
      
      Single_PK_Dataset <- filter(All_PK_Data, Temporal_ID == k) # filter by time_series_number
      dose <- Single_PK_Dataset$Dose_mg[1] # need to alter the model so that it takes dose as an input, probably have it accessed as part of Single_PK_Dataset within the MCMC
      end_time <- max(Single_PK_Dataset$Time)
      
      Albendazole_Time <- Single_PK_Dataset$Time[Single_PK_Dataset$Metabolite == "Alb"]  
      Albendazole_Conc <- Single_PK_Dataset$Converted_Concentration[Single_PK_Dataset$Metabolite == "Alb"]  
      Alb_data <- data.frame(Time = Albendazole_Time, Alb = Albendazole_Conc)
      
      Albendazole_Sulfoxide_Time <- Single_PK_Dataset$Time[Single_PK_Dataset$Metabolite == "AlbSO"]  
      Albendazole_Sulfoxide_Conc <- Single_PK_Dataset$Converted_Concentration[Single_PK_Dataset$Metabolite == "AlbSO"]  
      Alb_SO_data <- data.frame(Time = Albendazole_Sulfoxide_Time, Alb_SO = Albendazole_Sulfoxide_Conc)
      
      individual_dataset_MCMC_output <- raw_MCMC_output[[k]]
      median_k_abs <- median(individual_dataset_MCMC_output[, 1])
      median_bioavailability <- median(individual_dataset_MCMC_output[, 2])
      median_sigma <- median(individual_dataset_MCMC_output[, 3])
      median_k_alb_so <- median(individual_dataset_MCMC_output[, 4])
      median_k_alb <- median(individual_dataset_MCMC_output[, 5])
      model_runner <- Albendazole_PK_Model(k_abs = median_k_abs, bioavailability = median_bioavailability, sigma = median_sigma, k_alb_so = median_k_alb_so, dose = dose, k_alb = median_k_alb) 
      times <- seq(0, 100, length.out = 1000)
      median_out <- model_runner$run(times)
      max <- max(c(as.vector(median_out[, 5:6]), Alb_SO_data$Alb_SO))
      
      if (state[k] == "Fasted") {
        colour <- colours[1]
        Fasted_model_output_storage[Fasted_counter, ] <- median_out[, 6]
        Fasted_counter <- Fasted_counter + 1
        if (k == 1) {
          plot(median_out[, 1], median_out[, 6], type = "l", col = adjustcolor(colour, alpha.f = 0.2),  ylim = c(0, 1000), xlim = c(0, 48), ylab = "", xlab = xlab, las = 1, lwd = 2, cex.lab = text_size, cex.axis = text_size)
          title(ylab = ylab, line = 4, cex.lab = text_size)
        }
        else {
          lines(median_out[, 1], median_out[, 6], type = "l", col = adjustcolor(colour, alpha.f = 0.2),  ylim = c(0, max), xlim = c(0, 48), ylab = "", xlab = xlab, las = 1, lwd = 2)
        }
      } else if (state[k] == "Fatty_meal") {
        colour <- colours[2]
        Fatty_meal_model_output_storage[Fatty_meal_counter, ] <- median_out[, 6]
        Fatty_meal_counter <- Fatty_meal_counter + 1
        if (k == 1) {
          plot(median_out[, 1], median_out[, 6], type = "l", col = adjustcolor(colour, alpha.f = 0.2),  ylim = c(0, 1000), xlim = c(0, 48), ylab = "", xlab = xlab, las = 1, lwd = 2, cex.lab = text_size, cex.axis = text_size)
          title(ylab = ylab, line = 4, cex.lab = text_size)
        }
        else {
          lines(median_out[, 1], median_out[, 6], type = "l", col = adjustcolor(colour, alpha.f = 0.2),  ylim = c(0, max), xlim = c(0, 48), ylab = "", xlab = xlab, las = 1, lwd = 2)
        }
      } 
    }
    Fasted_mean_model_output <- apply(Fasted_model_output_storage, 2, mean)
    Fatty_meal_mean_model_output <- apply(Fatty_meal_model_output_storage, 2, mean)
    #Unclear_mean_model_output <- apply(Unclear_model_output_storage, 2, mean)
    lines(median_out[, 1], Fasted_mean_model_output, col = colours[1], lwd = 5)
    lines(median_out[, 1], Fatty_meal_mean_model_output, col = colours[2], lwd = 5)
    #lines(median_out[, 1], Unclear_mean_model_output, col = "grey", lwd = 5)
    legend("topright", inset = 0.02, legend = c(paste0("Fasted (n = ", Fasted_counter - 1, ")"), paste0("Fatty Meal (n = ", Fatty_meal_counter - 1, ")")), title = "State", col = c(colours[1], colours[2]), lty = 1, lwd = 5, cex = 1.25)
  }
  else if (colouring_variable == "Dose") {
    xlab <- "Time Since Treatment (Hours)"
    ylab <- "Concentration (ng/ml)"
    High_model_output_storage <- matrix(nrow = sum(dose_info == "High"), ncol = length(seq(0, 100, length.out = 1000)))
    Low_model_output_storage <- matrix(nrow = sum(dose_info == "Low"), ncol = length(seq(0, 100, length.out = 1000)))
    high_counter <- 1
    low_counter <- 1
    
    for (k in 1:max(All_PK_Data$Temporal_ID)) { 
      
      Single_PK_Dataset <- filter(All_PK_Data, Temporal_ID == k) # filter by time_series_number
      dose <- Single_PK_Dataset$Dose_mg[1] # need to alter the model so that it takes dose as an input, probably have it accessed as part of Single_PK_Dataset within the MCMC
      end_time <- max(Single_PK_Dataset$Time)
      
      Albendazole_Time <- Single_PK_Dataset$Time[Single_PK_Dataset$Metabolite == "Alb"]  
      Albendazole_Conc <- Single_PK_Dataset$Converted_Concentration[Single_PK_Dataset$Metabolite == "Alb"]  
      Alb_data <- data.frame(Time = Albendazole_Time, Alb = Albendazole_Conc)
      
      Albendazole_Sulfoxide_Time <- Single_PK_Dataset$Time[Single_PK_Dataset$Metabolite == "AlbSO"]  
      Albendazole_Sulfoxide_Conc <- Single_PK_Dataset$Converted_Concentration[Single_PK_Dataset$Metabolite == "AlbSO"]  
      Alb_SO_data <- data.frame(Time = Albendazole_Sulfoxide_Time, Alb_SO = Albendazole_Sulfoxide_Conc)
      
      individual_dataset_MCMC_output <- raw_MCMC_output[[k]]
      median_k_abs <- median(individual_dataset_MCMC_output[, 1])
      median_bioavailability <- median(individual_dataset_MCMC_output[, 2])
      median_sigma <- median(individual_dataset_MCMC_output[, 3])
      median_k_alb_so <- median(individual_dataset_MCMC_output[, 4])
      median_k_alb <- median(individual_dataset_MCMC_output[, 5])
      model_runner <- Albendazole_PK_Model(k_abs = median_k_abs, bioavailability = median_bioavailability, sigma = median_sigma, k_alb_so = median_k_alb_so, dose = dose, k_alb = median_k_alb) 
      times <- seq(0, 100, length.out = 1000)
      median_out <- model_runner$run(times)
      max <- max(c(as.vector(median_out[, 5:6]), Alb_SO_data$Alb_SO))
      
      if (dose_info[k] == "High") {
        colour <- colours[1]
        High_model_output_storage[high_counter, ] <- median_out[, 6]
        high_counter <- high_counter + 1
        if (k == 1) {
          plot(median_out[, 1], median_out[, 6], type = "l", col = adjustcolor(colour, alpha.f = 0.2),  ylim = c(0, 1000), xlim = c(0, 48), ylab = "", xlab = xlab, las = 1, lwd = 2, cex.lab = text_size, cex.axis = text_size)
          title(ylab = ylab, line = 4, cex.lab = text_size)
        } else {
          lines(median_out[, 1], median_out[, 6], type = "l", col = adjustcolor(colour, alpha.f = 0.2),  ylim = c(0, max), xlim = c(0, 48), ylab = "", xlab = xlab, las = 1, lwd = 2)
        }
      } else if (dose_info[k] == "Low") {
        colour <- colours[2]
        Low_model_output_storage[low_counter, ] <- median_out[, 6]
        low_counter <- low_counter + 1
        if (k == 1) {
          plot(median_out[, 1], median_out[, 6], type = "l", col = adjustcolor(colour, alpha.f = 0.2),  ylim = c(0, 1000), xlim = c(0, 48), ylab = "", xlab = xlab, las = 1, lwd = 2, cex.lab = text_size, cex.axis = text_size)
          title(ylab = ylab, line = 4, cex.lab = text_size)
        } else {
          lines(median_out[, 1], median_out[, 6], type = "l", col = adjustcolor(colour, alpha.f = 0.2),  ylim = c(0, max), xlim = c(0, 48), ylab = "", xlab = xlab, las = 1, lwd = 2)
        }
      }
    }
    High_mean_model_output <- apply(High_model_output_storage, 2, mean)
    Low_mean_model_output <- apply(Low_model_output_storage, 2, mean)
    lines(median_out[, 1], High_mean_model_output, col = colours[1], lwd = 5)
    lines(median_out[, 1], Low_mean_model_output, col = colours[2], lwd = 5)
    legend("topright", inset = 0.02, legend = c(paste0("High (n = ", high_counter - 1, ")"), paste0("Low (n = ", low_counter - 1, ")")), col = c(colours[1], colours[2]), title = "Dose",  lty = 1, lwd = 5, cex = 1.25)
  }
  else if (colouring_variable == "Infection") {
    xlab <- "Time Since Treatment (Hours)"
    ylab <- "Concentration (ng/ml)"
    Healthy_model_output_storage <- matrix(nrow = sum(infected == "Healthy"), ncol = length(seq(0, 100, length.out = 1000)))
    Infected_model_output_storage <- matrix(nrow = sum(infected == "Infected"), ncol = length(seq(0, 100, length.out = 1000)))
    healthy_counter <- 1
    infected_counter <- 1
    
    for (k in 1:max(All_PK_Data$Temporal_ID)) { 
      
      Single_PK_Dataset <- filter(All_PK_Data, Temporal_ID == k) # filter by time_series_number
      dose <- Single_PK_Dataset$Dose_mg[1] # need to alter the model so that it takes dose as an input, probably have it accessed as part of Single_PK_Dataset within the MCMC
      end_time <- max(Single_PK_Dataset$Time)
      
      Albendazole_Time <- Single_PK_Dataset$Time[Single_PK_Dataset$Metabolite == "Alb"]  
      Albendazole_Conc <- Single_PK_Dataset$Converted_Concentration[Single_PK_Dataset$Metabolite == "Alb"]  
      Alb_data <- data.frame(Time = Albendazole_Time, Alb = Albendazole_Conc)
      
      Albendazole_Sulfoxide_Time <- Single_PK_Dataset$Time[Single_PK_Dataset$Metabolite == "AlbSO"]  
      Albendazole_Sulfoxide_Conc <- Single_PK_Dataset$Converted_Concentration[Single_PK_Dataset$Metabolite == "AlbSO"]  
      Alb_SO_data <- data.frame(Time = Albendazole_Sulfoxide_Time, Alb_SO = Albendazole_Sulfoxide_Conc)
      
      individual_dataset_MCMC_output <- raw_MCMC_output[[k]]
      median_k_abs <- median(individual_dataset_MCMC_output[, 1])
      median_bioavailability <- median(individual_dataset_MCMC_output[, 2])
      median_sigma <- median(individual_dataset_MCMC_output[, 3])
      median_k_alb_so <- median(individual_dataset_MCMC_output[, 4])
      median_k_alb <- median(individual_dataset_MCMC_output[, 5])
      model_runner <- Albendazole_PK_Model(k_abs = median_k_abs, bioavailability = median_bioavailability, sigma = median_sigma, k_alb_so = median_k_alb_so, dose = dose, k_alb = median_k_alb) 
      times <- seq(0, 100, length.out = 1000)
      median_out <- model_runner$run(times)
      max <- max(c(as.vector(median_out[, 5:6]), Alb_SO_data$Alb_SO))
      
      if (infected[k] == "Healthy") {
        colour <- colours[1]
        Healthy_model_output_storage[healthy_counter, ] <- median_out[, 6]
        healthy_counter <- healthy_counter + 1
        if (k == 1) {
          plot(median_out[, 1], median_out[, 6], type = "l", col = adjustcolor(colour, alpha.f = 0.2),  ylim = c(0, 1000), xlim = c(0, 48), ylab = "", xlab = xlab, las = 1, lwd = 2, cex.lab = text_size, cex.axis = text_size)
          title(ylab = ylab, line = 4, cex.lab = text_size)
        } else {
          lines(median_out[, 1], median_out[, 6], type = "l", col = adjustcolor(colour, alpha.f = 0.2),  ylim = c(0, max), xlim = c(0, 48), ylab = "", xlab = xlab, las = 1, lwd = 2)
        }
      } else if (infected[k] == "Infected") {
        colour <- colours[2]
        Infected_model_output_storage[infected_counter, ] <- median_out[, 6]
        infected_counter <- infected_counter + 1
        if (k == 1) {
          plot(median_out[, 1], median_out[, 6], type = "l", col = adjustcolor(colour, alpha.f = 0.2),  ylim = c(0, 1000), xlim = c(0, 48), ylab = "", xlab = xlab, las = 1, lwd = 2, cex.lab = text_size, cex.axis = text_size)
          title(ylab = ylab, line = 4, cex.lab = text_size)
        } else {
          lines(median_out[, 1], median_out[, 6], type = "l", col = adjustcolor(colour, alpha.f = 0.2),  ylim = c(0, max), xlim = c(0, 48), ylab = "", xlab = xlab, las = 1, lwd = 2)
        }
      }
    }
    Healthy_mean_model_output <- apply(Healthy_model_output_storage, 2, mean)
    Infected_mean_model_output <- apply(Infected_model_output_storage, 2, mean)
    lines(median_out[, 1], Healthy_mean_model_output, col = colours[1], lwd = 5)
    lines(median_out[, 1], Infected_mean_model_output, col = colours[2], lwd = 5)
    legend("topright", inset = 0.02, legend = c(paste0("Healthy (n = ", healthy_counter - 1, ")"), paste0("Infected (n = ", infected_counter - 1, ")")), col = c(colours[1], colours[2]), title = "Infection",  lty = 1, lwd = 5, cex = 1.25)
  }
  else if (colouring_variable == "Drug") {
    
    xlab <- "Time Since Treatment (Hours)"
    ylab <- "Concentration (ng/ml)"
    Yes_model_output_storage <- matrix(nrow = sum(drugs == "Yes"), ncol = length(seq(0, 100, length.out = 1000)))
    None_model_output_storage <- matrix(nrow = sum(drugs == "None"), ncol = length(seq(0, 100, length.out = 1000)))
    yes_counter <- 1
    none_counter <- 1
    
    for (k in 1:max(All_PK_Data$Temporal_ID)) { 
      
      Single_PK_Dataset <- filter(All_PK_Data, Temporal_ID == k) # filter by time_series_number
      dose <- Single_PK_Dataset$Dose_mg[1] # need to alter the model so that it takes dose as an input, probably have it accessed as part of Single_PK_Dataset within the MCMC
      end_time <- max(Single_PK_Dataset$Time)
      
      Albendazole_Time <- Single_PK_Dataset$Time[Single_PK_Dataset$Metabolite == "Alb"]  
      Albendazole_Conc <- Single_PK_Dataset$Converted_Concentration[Single_PK_Dataset$Metabolite == "Alb"]  
      Alb_data <- data.frame(Time = Albendazole_Time, Alb = Albendazole_Conc)
      
      Albendazole_Sulfoxide_Time <- Single_PK_Dataset$Time[Single_PK_Dataset$Metabolite == "AlbSO"]  
      Albendazole_Sulfoxide_Conc <- Single_PK_Dataset$Converted_Concentration[Single_PK_Dataset$Metabolite == "AlbSO"]  
      Alb_SO_data <- data.frame(Time = Albendazole_Sulfoxide_Time, Alb_SO = Albendazole_Sulfoxide_Conc)
      
      individual_dataset_MCMC_output <- raw_MCMC_output[[k]]
      median_k_abs <- median(individual_dataset_MCMC_output[, 1])
      median_bioavailability <- median(individual_dataset_MCMC_output[, 2])
      median_sigma <- median(individual_dataset_MCMC_output[, 3])
      median_k_alb_so <- median(individual_dataset_MCMC_output[, 4])
      median_k_alb <- median(individual_dataset_MCMC_output[, 5])
      model_runner <- Albendazole_PK_Model(k_abs = median_k_abs, bioavailability = median_bioavailability, sigma = median_sigma, k_alb_so = median_k_alb_so, dose = dose, k_alb = median_k_alb) 
      times <- seq(0, 100, length.out = 1000)
      median_out <- model_runner$run(times)
      max <- max(c(as.vector(median_out[, 5:6]), Alb_SO_data$Alb_SO))
      
      if (drugs[k] == "Yes") {
        colour <- colours[1]
        Yes_model_output_storage[yes_counter, ] <- median_out[, 6]
        yes_counter <- yes_counter + 1
      } else if (drugs[k] == "None") {
        colour <- colours[2]
        None_model_output_storage[none_counter, ] <- median_out[, 6]
        none_counter <- none_counter + 1
      }
      if (k == 1) {
        plot(median_out[, 1], median_out[, 6], type = "l", col = adjustcolor(colour, alpha.f = 0.2),  ylim = c(0, 1000), xlim = c(0, 48), ylab = "", xlab = xlab, las = 1, lwd = 2, cex.lab = text_size, cex.axis = text_size)
        title(ylab = ylab, line = 4, cex.lab = text_size)
      } else {
        lines(median_out[, 1], median_out[, 6], type = "l", col = adjustcolor(colour, alpha.f = 0.2),  ylim = c(0, max), xlim = c(0, 48), ylab = "", xlab = xlab, las = 1, lwd = 2)
      }
    }
    Yes_mean_model_output <- apply(Yes_model_output_storage, 2, mean)
    None_mean_model_output <- apply(None_model_output_storage, 2, mean)
    lines(median_out[, 1], Yes_mean_model_output, col = colours[1], lwd = 5)
    lines(median_out[, 1], None_mean_model_output, col = colours[2], lwd = 5)
    legend("topright", inset = 0.02, legend = c(paste0("Yes (n = ", yes_counter - 1, ")"), paste0("None (n = ", none_counter - 1, ")")), col = c(colours[1], colours[2]), title = "Drug",  lty = 1, lwd = 5, cex = 1.25)
  }
  else if (colouring_variable == "Age") {
    
    xlab <- "Time Since Treatment (Hours)"
    ylab <- "Concentration (ng/ml)"
    Young_model_output_storage <- matrix(nrow = sum(age_group == "Young", na.rm = TRUE), ncol = length(seq(0, 100, length.out = 1000)))
    Old_model_output_storage <- matrix(nrow = sum(age_group == "Old", na.rm = TRUE), ncol = length(seq(0, 100, length.out = 1000)))
    Young_counter <- 1
    Old_counter <- 1

    for (k in 1:max(All_PK_Data$Temporal_ID)) { 
      
      Single_PK_Dataset <- filter(All_PK_Data, Temporal_ID == k) # filter by time_series_number
      dose <- Single_PK_Dataset$Dose_mg[1] # need to alter the model so that it takes dose as an input, probably have it accessed as part of Single_PK_Dataset within the MCMC
      end_time <- max(Single_PK_Dataset$Time)
      
      Albendazole_Time <- Single_PK_Dataset$Time[Single_PK_Dataset$Metabolite == "Alb"]  
      Albendazole_Conc <- Single_PK_Dataset$Converted_Concentration[Single_PK_Dataset$Metabolite == "Alb"]  
      Alb_data <- data.frame(Time = Albendazole_Time, Alb = Albendazole_Conc)
      
      Albendazole_Sulfoxide_Time <- Single_PK_Dataset$Time[Single_PK_Dataset$Metabolite == "AlbSO"]  
      Albendazole_Sulfoxide_Conc <- Single_PK_Dataset$Converted_Concentration[Single_PK_Dataset$Metabolite == "AlbSO"]  
      Alb_SO_data <- data.frame(Time = Albendazole_Sulfoxide_Time, Alb_SO = Albendazole_Sulfoxide_Conc)
      
      individual_dataset_MCMC_output <- raw_MCMC_output[[k]]
      median_k_abs <- median(individual_dataset_MCMC_output[, 1])
      median_bioavailability <- median(individual_dataset_MCMC_output[, 2])
      median_sigma <- median(individual_dataset_MCMC_output[, 3])
      median_k_alb_so <- median(individual_dataset_MCMC_output[, 4])
      median_k_alb <- median(individual_dataset_MCMC_output[, 5])
      model_runner <- Albendazole_PK_Model(k_abs = median_k_abs, bioavailability = median_bioavailability, sigma = median_sigma, k_alb_so = median_k_alb_so, dose = dose, k_alb = median_k_alb) 
      times <- seq(0, 100, length.out = 1000)
      median_out <- model_runner$run(times)
      max <- max(c(as.vector(median_out[, 5:6]), Alb_SO_data$Alb_SO))
      
      if (age_group[k] == "Young") {
        colour <- colours[1]
        Young_model_output_storage[Young_counter, ] <- median_out[, 6]
        Young_counter <- Young_counter + 1
        if (k == 1) {
          plot(median_out[, 1], median_out[, 6], type = "l", col = adjustcolor(colour, alpha.f = 0.2),  ylim = c(0, 1000), xlim = c(0, 48), ylab = "", xlab = xlab, las = 1, lwd = 2, cex.lab = text_size, cex.axis = text_size)
          title(ylab = ylab, line = 4, cex.lab = text_size)
        }
        else {
          lines(median_out[, 1], median_out[, 6], type = "l", col = adjustcolor(colour, alpha.f = 0.2),  ylim = c(0, max), xlim = c(0, 48), ylab = "", xlab = xlab, las = 1, lwd = 2)
        }
      } else if (age_group[k] == "Old") {
        colour <- colours[2]
        Old_model_output_storage[Old_counter, ] <- median_out[, 6]
        Old_counter <- Old_counter + 1
        if (k == 1) {
          plot(median_out[, 1], median_out[, 6], type = "l", col = adjustcolor(colour, alpha.f = 0.2),  ylim = c(0, 1000), xlim = c(0, 48), ylab = "", xlab = xlab, las = 1, lwd = 2, cex.lab = text_size, cex.axis = text_size)
          title(ylab = ylab, line = 4, cex.lab = text_size)
        }
        else {
          lines(median_out[, 1], median_out[, 6], type = "l", col = adjustcolor(colour, alpha.f = 0.2),  ylim = c(0, max), xlim = c(0, 48), ylab = "", xlab = xlab, las = 1, lwd = 2)
        }
      } 
    }
    Young_model_output_storage <- apply(Young_model_output_storage, 2, mean)
    Old_model_output_storage <- apply(Old_model_output_storage, 2, mean)
    #Unclear_mean_model_output <- apply(Unclear_model_output_storage, 2, mean)
    lines(median_out[, 1], Young_model_output_storage, col = colours[1], lwd = 5)
    lines(median_out[, 1], Old_model_output_storage, col = colours[2], lwd = 5)
    #lines(median_out[, 1], Unclear_mean_model_output, col = "grey", lwd = 5)
    legend("topright", inset = 0.02, legend = c(paste0("Young (n = ", Young_counter - 1,")"), paste0("Old (n = ", Old_counter - 1, ")")), title = "Age Group", col = c(colours[1], colours[2]), lty = 1, lwd = 5, cex = 1.25)
  }
  else if (colouring_variable == "Weight") {
    
    xlab <- "Time Since Treatment (Hours)"
    ylab <- "Concentration (ng/ml)"
    Light_model_output_storage <- matrix(nrow = sum(weight_group == "Light", na.rm = TRUE), ncol = length(seq(0, 100, length.out = 1000)))
    Heavy_model_output_storage <- matrix(nrow = sum(weight_group == "Heavy", na.rm = TRUE), ncol = length(seq(0, 100, length.out = 1000)))
    Light_counter <- 1
    Heavy_counter <- 1
    
    for (k in 1:max(All_PK_Data$Temporal_ID)) { 
      
      Single_PK_Dataset <- filter(All_PK_Data, Temporal_ID == k) # filter by time_series_number
      dose <- Single_PK_Dataset$Dose_mg[1] # need to alter the model so that it takes dose as an input, probably have it accessed as part of Single_PK_Dataset within the MCMC
      end_time <- max(Single_PK_Dataset$Time)
      
      Albendazole_Time <- Single_PK_Dataset$Time[Single_PK_Dataset$Metabolite == "Alb"]  
      Albendazole_Conc <- Single_PK_Dataset$Converted_Concentration[Single_PK_Dataset$Metabolite == "Alb"]  
      Alb_data <- data.frame(Time = Albendazole_Time, Alb = Albendazole_Conc)
      
      Albendazole_Sulfoxide_Time <- Single_PK_Dataset$Time[Single_PK_Dataset$Metabolite == "AlbSO"]  
      Albendazole_Sulfoxide_Conc <- Single_PK_Dataset$Converted_Concentration[Single_PK_Dataset$Metabolite == "AlbSO"]  
      Alb_SO_data <- data.frame(Time = Albendazole_Sulfoxide_Time, Alb_SO = Albendazole_Sulfoxide_Conc)
      
      individual_dataset_MCMC_output <- raw_MCMC_output[[k]]
      median_k_abs <- median(individual_dataset_MCMC_output[, 1])
      median_bioavailability <- median(individual_dataset_MCMC_output[, 2])
      median_sigma <- median(individual_dataset_MCMC_output[, 3])
      median_k_alb_so <- median(individual_dataset_MCMC_output[, 4])
      median_k_alb <- median(individual_dataset_MCMC_output[, 5])
      model_runner <- Albendazole_PK_Model(k_abs = median_k_abs, bioavailability = median_bioavailability, sigma = median_sigma, k_alb_so = median_k_alb_so, dose = dose, k_alb = median_k_alb) 
      times <- seq(0, 100, length.out = 1000)
      median_out <- model_runner$run(times)
      max <- max(c(as.vector(median_out[, 5:6]), Alb_SO_data$Alb_SO))
      
      if (weight_group[k] == "Light") {
        colour <- colours[1]
        Light_model_output_storage[Light_counter, ] <- median_out[, 6]
        Light_counter <- Light_counter + 1
        if (k == 1) {
          plot(median_out[, 1], median_out[, 6], type = "l", col = adjustcolor(colour, alpha.f = 0.2),  ylim = c(0, 1000), xlim = c(0, 48), ylab = "", xlab = xlab, las = 1, lwd = 2, cex.lab = text_size, cex.axis = text_size)
          title(ylab = ylab, line = 4, cex.lab = text_size)
        }
        else {
          lines(median_out[, 1], median_out[, 6], type = "l", col = adjustcolor(colour, alpha.f = 0.2),  ylim = c(0, max), xlim = c(0, 48), ylab = "", xlab = xlab, las = 1, lwd = 2)
        }
      } else if (weight_group[k] == "Heavy") {
        colour <- colours[2]
        Heavy_model_output_storage[Heavy_counter, ] <- median_out[, 6]
        Heavy_counter <- Heavy_counter + 1
        if (k == 1) {
          plot(median_out[, 1], median_out[, 6], type = "l", col = adjustcolor(colour, alpha.f = 0.2),  ylim = c(0, 1000), xlim = c(0, 48), ylab = "", xlab = xlab, las = 1, lwd = 2, cex.lab = text_size, cex.axis = text_size)
          title(ylab = ylab, line = 4, cex.lab = text_size)
        }
        else {
          lines(median_out[, 1], median_out[, 6], type = "l", col = adjustcolor(colour, alpha.f = 0.2),  ylim = c(0, max), xlim = c(0, 48), ylab = "Concentration (ng/ml)", xlab = "Time (Hours)", las = 1, lwd = 2)
        }
      } 
    }
    Light_model_output_storage <- apply(Light_model_output_storage, 2, mean)
    Heavy_model_output_storage <- apply(Heavy_model_output_storage, 2, mean)
    lines(median_out[, 1], Light_model_output_storage, col = colours[1], lwd = 5)
    lines(median_out[, 1], Heavy_model_output_storage, col = colours[2], lwd = 5)
    legend("topright", inset = 0.02, legend = c(paste0("Light (n = ", Light_counter - 1,")"), paste0("Heavy (n = ", Heavy_counter - 1, ")")), title = "Weight Group", col = c(colours[1], colours[2]), lty = 1, lwd = 5, cex = 1.25)
  }
}

MCMC_output_plotting_function <- function(raw_MCMC_output, time_series_index, plot_limit, plot_separately, plot_both, credible_intervals) {
  k <- time_series_index
  Single_PK_Dataset <- filter(All_PK_Data, Temporal_ID == k) # filter by time_series_number
  end_time <- max(Single_PK_Dataset$Time)
  dose <- Single_PK_Dataset$Dose_mg[1] # need to alter the model so that it takes dose as an input, probably have it accessed as part of Single_PK_Dataset within the MCMC
  
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
  
  chain <- raw_MCMC_output[[k]][20000:40000, ]
  median_k_abs <- median(chain[, 1])
  median_bioavailability <- median(chain[, 2])
  median_sigma <- median(chain[, 3])
  median_k_alb_so <- median(chain[, 4])
  median_k_alb <- median(chain[, 5])
  model_runner <- Albendazole_PK_Model(k_abs = median_k_abs, bioavailability = median_bioavailability, sigma = median_sigma, k_alb_so = median_k_alb_so, k_alb = median_k_alb, dose = dose) 
  times <- seq(0, 100, length.out = 1000)
  median_out <- model_runner$run(times)
  max <- max(c(as.vector(median_out[, 5:6]), Alb_SO_data$Alb_SO))
  
  # Generating the 95% Credible Intervals
  if (credible_intervals) {
    #size_of_burned_chain <- number_iterations - (burn_in * number_iterations)
    thinned_chain <- chain # [seq(1, size_of_burned_chain, 10), ]
    Alb_SO_Storage_matrix <- matrix(nrow = length(thinned_chain[, 1]), ncol = length(median_out[, 1]))
    Alb_Storage_matrix <- matrix(nrow = length(thinned_chain[, 1]), ncol = length(median_out[, 1]))
    for (i in 1:length(thinned_chain[, 1])) {
      k_abs <- thinned_chain[i, 1]
      bioavailability <- thinned_chain[i, 2]
      sigma <- thinned_chain[i, 3]
      k_alb_so <- thinned_chain[i, 4]
      k_alb <- thinned_chain[i, 5]
      model_runner <- Albendazole_PK_Model(k_abs = k_abs, bioavailability = bioavailability, sigma = sigma, k_alb_so = k_alb_so, dose = dose, k_alb = k_alb)
      out <- model_runner$run(times)
      Alb_Storage_matrix[i, ] <- out[, 5]
      Alb_SO_Storage_matrix[i, ] <- out[, 6]
    }
    Alb_Credible_Upper <- apply(Alb_Storage_matrix, 2, quantile, 0.975)
    Alb_Credible_Lower <- apply(Alb_Storage_matrix, 2, quantile, 0.025)
    Alb_SO_Credible_Upper <- apply(Alb_SO_Storage_matrix, 2, quantile, 0.975)
    Alb_SO_Credible_Lower <- apply(Alb_SO_Storage_matrix, 2, quantile, 0.025)
  }
  
  # Plotting the Median Output Along With the 95% Credible Intervals 
  if (plot_separately == FALSE) {
    if (plot_limit == TRUE) {
      if (plot_both == TRUE) {
        plot(median_out[, 1], median_out[, 5], type = "l", col = "#E5005F",  ylim = c(0, max), xlim = c(0, 48), ylab = "", xlab = "", las = 1, lwd = 2)
        lines(median_out[, 1], median_out[, 6], type = "l", col = "#7609BA", lwd = 2)
        points(Alb_data$Time, Alb_data$Alb, pch = 20, col = "#E5005F", cex = 2)
        points(Alb_SO_data$Time, Alb_SO_data$Alb_SO, pch = 20, col = "#7609BA", cex = 2)
        if(credible_intervals) {
          polygon(c(median_out[, 1], rev(median_out[, 1])), c(Alb_Credible_Lower, rev(Alb_Credible_Upper)), col = adjustcolor("#E5005F", alpha.f = 0.2), border = NA)
          polygon(c(median_out[, 1], rev(median_out[, 1])), c(Alb_SO_Credible_Lower, rev(Alb_SO_Credible_Upper)), col = adjustcolor("#7609BA", alpha.f = 0.2), border = NA)
        }
      } else {
        plot(median_out[, 1], median_out[, 6], type = "l", col = "#7609BA", ylim = c(0, max), xlim = c(0, 48), ylab = "", xlab = "", las = 1, lwd = 2)
        points(Alb_SO_data$Time, Alb_SO_data$Alb_SO, pch = 20, col = "#7609BA", cex = 2)
        if(credible_intervals) {
          polygon(c(median_out[, 1], rev(median_out[, 1])), c(Alb_SO_Credible_Lower, rev(Alb_SO_Credible_Upper)), col = adjustcolor("#7609BA", alpha.f = 0.2), border = NA)
        }
     }
    } else {
      if (plot_both == TRUE) {
        plot(median_out[, 1], median_out[, 5], type = "l", col = "#E5005F",  ylim = c(0, max), xlim = c(0, end_time), ylab = "", xlab = "", las = 1, lwd = 2)
        lines(median_out[, 1], median_out[, 6], type = "l", col = "#7609BA", lwd = 2)
        points(Alb_data$Time, Alb_data$Alb, pch = 20, col = "#E5005F", cex = 2)
        points(Alb_SO_data$Time, Alb_SO_data$Alb_SO, pch = 20, col = "#7609BA", cex = 2)
        if (credible_intervals) {
          polygon(c(median_out[, 1], rev(median_out[, 1])), c(Alb_Credible_Lower, rev(Alb_Credible_Upper)), col = adjustcolor("#E5005F", alpha.f = 0.2), border = NA)
          polygon(c(median_out[, 1], rev(median_out[, 1])), c(Alb_SO_Credible_Lower, rev(Alb_SO_Credible_Upper)), col = adjustcolor("#7609BA", alpha.f = 0.2), border = NA)
        }
      } else {
        plot(median_out[, 1], median_out[, 6], type = "l", col = "#7609BA",  ylim = c(0, max), xlim = c(0, end_time), ylab = "", xlab = "", las = 1, lwd = 2)
        points(Alb_SO_data$Time, Alb_SO_data$Alb_SO, pch = 20, col = "#7609BA", cex = 2)
        if (credible_intervals) {
          polygon(c(median_out[, 1], rev(median_out[, 1])), c(Alb_SO_Credible_Lower, rev(Alb_SO_Credible_Upper)), col = adjustcolor("#7609BA", alpha.f = 0.2), border = NA)
        }
      }
    }
  } 
  else {
    if (plot_limit == TRUE) {
      if (plot_both == TRUE) {
        plot(median_out[, 1], median_out[, 5], type = "l", col = "#E5005F",  ylim = c(0, max(c(median_out[, 4], Alb_data$Alb))), xlim = c(0, 48), ylab = "", xlab = "", las = 1, lwd = 2)
        points(Alb_data$Time, Alb_data$Alb, pch = 20, col = "#E5005F", cex = 2)
        if (credible_intervals) {
          polygon(c(median_out[, 1], rev(median_out[, 1])), c(Alb_Credible_Lower, rev(Alb_Credible_Upper)), col = adjustcolor("#E5005F", alpha.f = 0.2), border = NA)
        }
        plot(median_out[, 1], median_out[, 6], type = "l", col = "#7609BA", ylim = c(0, max), xlim = c(0, 48), ylab = "", xlab = "", las = 1, lwd = 2)
        points(Alb_SO_data$Time, Alb_SO_data$Alb_SO, pch = 20, col = "#7609BA", cex = 2)
        if (credible_intervals) {
          polygon(c(median_out[, 1], rev(median_out[, 1])), c(Alb_SO_Credible_Lower, rev(Alb_SO_Credible_Upper)), col = adjustcolor("#7609BA", alpha.f = 0.2), border = NA)
        }
      } else {
        plot(median_out[, 1], median_out[, 6], type = "l", col = "#7609BA", ylim = c(0, max), xlim = c(0, 48), ylab = "", xlab = "", las = 1, lwd = 2)
        points(Alb_SO_data$Time, Alb_SO_data$Alb_SO, pch = 20, col = "#7609BA", cex = 2)
        if (credible_intervals) {
          polygon(c(median_out[, 1], rev(median_out[, 1])), c(Alb_SO_Credible_Lower, rev(Alb_SO_Credible_Upper)), col = adjustcolor("#7609BA", alpha.f = 0.2), border = NA)
        }
      }
    } else {
      if (plot_both == TRUE) {
        plot(median_out[, 1], median_out[, 5], type = "l", col = "#E5005F",  ylim = c(0, max(c(median_out[, 4], Alb_data$Alb))), xlim = c(0, end_time), ylab = "Concentration (ng/ml)", xlab = "Time (Hours)", las = 1, lwd = 2)
        points(Alb_data$Time, Alb_data$Alb, pch = 20, col = "#E5005F", cex = 2)
        if (credible_intervals) {
          polygon(c(median_out[, 1], rev(median_out[, 1])), c(Alb_Credible_Lower, rev(Alb_Credible_Upper)), col = adjustcolor("#E5005F", alpha.f = 0.2), border = NA)
        }
        plot(median_out[, 1], median_out[, 6], type = "l", col = "#7609BA",   ylim = c(0, max), xlim = c(0, end_time), ylab = "", xlab = "", las = 1, lwd = 2)
        points(Alb_SO_data$Time, Alb_SO_data$Alb_SO, pch = 20, col = "#7609BA", cex = 2)
        if (credible_intervals) {
          polygon(c(median_out[, 1], rev(median_out[, 1])), c(Alb_SO_Credible_Lower, rev(Alb_SO_Credible_Upper)), col = adjustcolor("#7609BA", alpha.f = 0.2), border = NA)
        }
      } else {
        plot(median_out[, 1], median_out[, 6], type = "l", col = "#7609BA",  ylim = c(0, max), xlim = c(0, end_time), ylab = "", xlab = "", las = 1, lwd = 2)
        points(Alb_SO_data$Time, Alb_SO_data$Alb_SO, pch = 20, col = "#7609BA", cex = 2)
        if (credible_intervals) {
          polygon(c(median_out[, 1], rev(median_out[, 1])), c(Alb_SO_Credible_Lower, rev(Alb_SO_Credible_Upper)), col = adjustcolor("#7609BA", alpha.f = 0.2), border = NA)
        }
      }
    }
  }
  print(k)
}
