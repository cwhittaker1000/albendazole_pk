Albendazole_PK_Model <- odin::odin({
 
  ## Derivatives
  deriv(Gut_1) <- - k_abs * Gut_1
  deriv(Gut_2) <- (k_abs * Gut_1) - (k_abs * Gut_2)
  deriv(Liver) <- (k_abs * Gut_2) - (sigma * Liver) - (15 * Liver) + (15 * Blood_Alb)
  deriv(Blood_Alb) <- (15 * Liver) - (15 * Blood_Alb) - (k_alb * Blood_Alb)
  deriv(Blood_Alb_SO) <- (sigma * Liver) - (k_alb_so * Blood_Alb_SO)
  
  ## Initial conditions
  initial(Gut_1) <- gut_1
  initial(Gut_2) <- gut_2
  initial(Liver) <- liver
  initial(Blood_Alb) <- blood_alb
  initial(Blood_Alb_SO) <- blood_alb_so
  
  ## Initial condition parameters
  gut_1 <- user()
  gut_2 <- user()
  liver <- user()
  blood_alb <- user()
  blood_alb_so <- user()
  
  ## Parameters
  k_abs <- user()
  sigma <- user()
  k_alb <- user()
  k_alb_so <- user()

})
