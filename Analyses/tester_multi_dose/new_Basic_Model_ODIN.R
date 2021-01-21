Albendazole_PK_Model <- odin::odin({
 
  ## Derivatives
  deriv(Gut_1) <- - k_abs * Gut_1
  deriv(Gut_2) <- (k_abs * Gut_1) - (k_abs * Gut_2)
  deriv(Liver) <- (k_abs * Gut_2) - (sigma * Liver) - (15 * Liver) + (15 * Blood_Alb)
  deriv(Blood_Alb) <- (15 * Liver) - (15 * Blood_Alb)
  deriv(Blood_Alb_SO) <- (sigma * Liver) - (k_alb_so * Blood_Alb_SO) - (k_periph * Blood_Alb_SO) + (k_central * Peripheral)
  deriv(Peripheral) <- (k_periph * Blood_Alb_SO) - (k_central * Peripheral)
  
  ## Initial conditions
  initial(Gut_1) <- (bioavailability/1e+5) * (dose * 1e+6)
  initial(Gut_2) <- 0.0
  initial(Liver) <- 0.0
  initial(Blood_Alb) <- 0.0
  initial(Blood_Alb_SO) <- 0.0
  initial(Peripheral) <- 0.0
  
  ## Parameters
  dose <- user()
  k_abs <- user()
  bioavailability <- user()
  sigma <- user()
  k_alb_so <- user()
  k_periph <- user()
  k_central <- user()
})


model_runner <- Albendazole_PK_Model(k_abs = 0.5, bioavailability = 0.03, sigma = 2, 
                                     k_alb_so = 0.1, k_periph = 0, k_central = 0, dose = 400)
times <- seq(0, 48, length.out = 1000)
y <- model_runner$run(times)
plot(y)

par(mfrow = c(1, 1))
plot(y[, 1], y[, 6], type = "l")
lines(y[, 1], y[, 5], type = "l", col = "red")
