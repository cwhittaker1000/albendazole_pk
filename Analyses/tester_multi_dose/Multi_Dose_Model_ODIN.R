library(DescTools)

Albendazole_PK_Model <- odin::odin({
  
  ## Derivatives
  deriv(Gut_1) <- - k_abs * Gut_1
  deriv(Gut_2) <- (k_abs * Gut_1) - (k_abs * Gut_2)
  deriv(Liver) <- (bioavailability * k_abs * Gut_2) - (sigma * Liver) - (15 * Liver) + (15 * Blood_Alb)
  deriv(Blood_Alb) <- (15 * Liver) - (15 * Blood_Alb) - (k_alb * Blood_Alb)
  deriv(Blood_Alb_SO) <- (sigma * Liver) - (k_alb_so * Blood_Alb_SO) - (regime_k * k_alb_so * Blood_Alb_SO)
  
  ## Initial conditions
  initial(Gut_1) <- dose
  initial(Gut_2) <- gut_2
  initial(Liver) <- liver
  initial(Blood_Alb) <- blood_alb
  initial(Blood_Alb_SO) <- blood_alb_so
  
  ## Initial condition parameters
  dose <- user()
  gut_2 <- user()
  liver <- user()
  blood_alb <- user()
  blood_alb_so <- user()
  regime_k <- user()

  ## Parameters
  k_abs <- user()
  bioavailability <- user()
  sigma <- user()
  k_alb_so <- user()
  k_alb <- user()
  
})


extract_final_values <- function(model_output) {
  end <- dim(model_output)[1]
  dose <- model_output[end, 2]
  gut_2 <- model_output[end, 3]
  liver <- model_output[end, 4]
  blood_alb <- model_output[end, 5]
  blood_alb_so <- model_output[end, 6]
  return(c(dose, gut_2, liver, blood_alb, blood_alb_so))
}

# 0, 24, 48

run_multi_dose <- function(parameters, dose_amounts, dose_times) {
  
  doses <- length(dose_amounts)
  k_abs <- parameters["k_abs"]
  bioavailability <- parameters["bioavailability"]
  sigma <- parameters["sigma"]
  k_alb_so <- parameters["k_alb_so"]
  k_alb <- parameters["k_alb"]
  
  for (i in 1:doses) {
    if (i == 1) {
      times <- seq(0, dose_times[2], length.out = 1 + (dose_times[2]/0.1))
      times_proper <- times
      model_runner <- Albendazole_PK_Model(k_abs = k_abs, bioavailability = bioavailability, sigma = sigma, k_alb_so = k_alb_so, k_alb = k_alb, 
                                           dose = dose_amounts[1], gut_2 = 0, liver = 0, blood_alb = 0, blood_alb_so = 0, 
                                           regime_k = 0)
      y_one <- model_runner$run(times)
      end <- dim(y_one)[1]
      vals <- extract_final_values(y_one)
      dose <- vals[1]
      gut_2 <- vals[2]
      liver <- vals[3]
      blood_alb <- vals[4]
      blood_alb_so <- vals[5]
      print(i)
    }
    else if ((i != 1) & (i < doses)) {
      times <- seq(0, dose_times[i + 1] - dose_times[i], length.out = 1 + (dose_times[i + 1] - dose_times[i])/0.1)
      times_proper <- c(times_proper[-length(times_proper)], times + max(times_proper))
      model_runner <- Albendazole_PK_Model(k_abs = k_abs, bioavailability = bioavailability, sigma = sigma, k_alb_so = k_alb_so, k_alb = k_alb, 
                                           dose = dose + dose_amounts[i], gut_2 = gut_2, liver = liver, blood_alb = blood_alb, blood_alb_so = blood_alb_so, 
                                           regime_k = 0)
      y <- model_runner$run(times)
      end <- dim(y)[1]
      vals <- extract_final_values(y)
      dose <- vals[1]
      gut_2 <- vals[2]
      liver <- vals[3]
      blood_alb <- vals[4]
      blood_alb_so <- vals[5]
      y_one <- rbind(y_one, y[-1, ])
      print(i)
    } else {
      times <- seq(0, 24, length.out = 241)
      times_proper <- c(times_proper[-length(times_proper)], times + max(times_proper))
      model_runner <- Albendazole_PK_Model(k_abs = k_abs, bioavailability = bioavailability, sigma = sigma, k_alb_so = k_alb_so, k_alb = k_alb, 
                                           dose = dose + dose_amounts[i], gut_2 = gut_2, liver = liver, blood_alb = blood_alb, blood_alb_so = blood_alb_so, 
                                           regime_k = 0)
      y <- model_runner$run(times)
      y_one <- rbind(y_one, y[-1, ])
      print(i)
      
    }
  }
  return(cbind(times_proper, y_one))
}

kloop <- tester_multi_dose[tester_multi_dose$ts == 3, ]
parameters <- c(k_abs = 3,
                bioavailability = 0.02,
                sigma = 70,
                k_alb_so = 0.25,
                k_alb = 0.12)
blopp <- run_multi_dose(parameters, c(80000, 80000, 80000), c(0, 24, 48))
plot(blopp[, 1], blopp[, 7], type = "l") #, ylim = c(0, max(c(blopp$conc, blopp[, 7]))))
points(kloop$time, kloop$conc, pch = 20)

times <- seq(0, 24, length.out = 100)
model_runner <- Albendazole_PK_Model(k_abs = 3, bioavailability = 0.03, sigma = 50, k_alb_so = 0.2, k_alb = 0.12, 
                                     dose = 45000, gut_2 = 0, liver = 0, blood_alb = 0, blood_alb_so = 0, 
                                     regime_k = 0)
y_one <- model_runner$run(times)
plot(y_one[, 6], type = "l")
extract_final_values(y_one)
end <- dim(y_one)[1]
dose <- y_one[end, 2]
gut_2 <- y_one[end, 3]
liver <- y_one[end, 4]
blood_alb <- y_one[end, 5]
blood_alb_so <- y_one[end, 6]
model_runner <- Albendazole_PK_Model(k_abs = 3, bioavailability = 0.01, sigma = 50, k_alb_so = 0.2, k_alb = 0.12, 
                                     dose = (dose + 45000), gut_2 = gut_2, liver = liver, blood_alb = blood_alb, blood_alb_so = blood_alb_so, 
                                     regime_k = 1)
y_two <- model_runner$run(times[-1])
y_overall <- rbind(y_one, y_two)
plot(y_overall[, 6], type = "l")




bloop <- seq(1:200)
n <- 0.5
half <- 25
hill <- 1 / (1 + (half/bloop)^n)
plot(hill, ylim = c(0, 1))

