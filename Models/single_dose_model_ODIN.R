Albendazole_PK_Model <- odin::odin({
 
  ## Derivatives
  deriv(Gut_1) <- - k_abs * Gut_1
  deriv(Gut_2) <- (k_abs * Gut_1) - (k_abs * Gut_2)
  deriv(Liver) <- (k_abs * Gut_2) - (sigma * Liver) - (15 * Liver) + (15 * Blood_Alb)
  deriv(Blood_Alb) <- (15 * Liver) - (15 * Blood_Alb) - (k_alb * Blood_Alb)
  deriv(Blood_Alb_SO) <- (sigma * Liver) - (k_alb_so * Blood_Alb_SO)
  
  ## Initial conditions
  initial(Gut_1) <- (bioavailability/1e+5) * (dose * 1e+6)
  initial(Gut_2) <- 0.0
  initial(Liver) <- 0.0
  initial(Blood_Alb) <- 0.0
  initial(Blood_Alb_SO) <- 0.0
  
  ## Parameters
  dose <- user()
  k_abs <- user()
  bioavailability <- user()
  sigma <- user()
  k_alb_so <- user()
  k_alb <- user()
})
