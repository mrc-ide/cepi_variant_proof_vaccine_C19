### ------------------------------------------------------------------------ ###
###     Code to generate vaccine parameters: 
###     1. Coverage matricesper dos
###     2. Next dose priority matrix
###     3. Daily vaccines vector
### ------------------------------------------------------------------------ ###


get_vaccine_strategy <- function(days_to_vacc_vec,  #Vector with time point where vaccination pace changes 
                                 doses_per_day,     #vector with the doses per day 
                                 time_period, 
                                 age_groups_covered_vec,
                                 vaccine_doses, 
                                 vacc_per_week,
                                 matrix_d3,
                                 matrix_booster,
                                 matrix_primary){
  
  # 1. Priority matrices for each dose ####
  
  # Reading coverage matrices for primary series, first boost and onwards boosters 
  strategy_matrixd1_d2 <- as.matrix(read.csv(matrix_primary))
  strategy_matrix_d3 <- as.matrix(read.csv(matrix_d3)) 
  strategy_matrix_booster <- as.matrix(read.csv(matrix_booster)) 
  
  
  
  vaccine_coverage_strategy <- vector("list" , length = vaccine_doses)
  
  if (vaccine_doses == 2) {
    vaccine_coverage_strategy[[1]] <- strategy_matrixd1_d2[1:age_groups_covered_vec[1],]
    vaccine_coverage_strategy[[2]] <- strategy_matrixd1_d2[1:age_groups_covered_vec[1],]
    
  } else if (vaccine_doses == 3 ) {
    
    vaccine_coverage_strategy[[1]] <- strategy_matrixd1_d2[1:age_groups_covered_vec[1],]
    vaccine_coverage_strategy[[2]] <-strategy_matrixd1_d2[1:age_groups_covered_vec[1],]
    vaccine_coverage_strategy[[3]] <-strategy_matrix_d3[1:age_groups_covered_vec[2],]
    
  }  else if (vaccine_doses > 3){
    
    vaccine_coverage_strategy[[1]] <- strategy_matrixd1_d2[1:age_groups_covered_vec[1],]
    vaccine_coverage_strategy[[2]] <-strategy_matrixd1_d2[1:age_groups_covered_vec[1],]
    vaccine_coverage_strategy[[3]] <-strategy_matrix_d3[1:age_groups_covered_vec[2],]
    
    for (i in 4:length(vaccine_coverage_strategy)) {
      vaccine_coverage_strategy[[i]] <- strategy_matrix_booster[1:age_groups_covered_vec[(i-1)],]
    }
    
  }
  
  
  #2. Next dose priority matrix ####
  
  next_dose_priority <- matrix(data = 1, nrow = vaccine_doses - 1, ncol = ncol(vaccine_coverage_strategy[[1]]))
  
  
  if (vaccine_doses == 2) {
    next_dose_priority <- matrix(data = 1, nrow = vaccine_doses - 1, ncol = ncol(vaccine_coverage_strategy[[1]]))
  } else if (vaccine_doses ==3 ) {
    
    next_dose_priority <- matrix(data = 0, nrow = vaccine_doses - 1, ncol = ncol(vaccine_coverage_strategy[[1]]))
    next_dose_priority[1,(17 - age_groups_covered_vec[1] + 1):17] <- 1
    
  } else if (vaccine_doses >3) {
    
    next_dose_priority <- matrix(data = 0, nrow = vaccine_doses - 1, ncol = ncol(vaccine_coverage_strategy[[1]]))
    next_dose_priority[1,(17 - age_groups_covered_vec[1] + 1):17] <- 1
    
    for (i in 3:length(next_dose_priority[,1])){
      next_dose_priority[i,(17 - age_groups_covered_vec[(i)] + 1):17] <- 1
    }
  }
  
  # 3. Creating vector of doses per day.  ####
  
  vaccine_set <- c(rep(0, days_to_vacc_vec[1])) # No vaccination until vaccination start day 
  
  #If there are more than one vaccination paces 
  if (length(doses_per_day) >1){
    
    for (i in 1:(length(doses_per_day)-1)){
      vaccine_set <- c(vaccine_set, rep(doses_per_day[i], (days_to_vacc_vec[i+1] - days_to_vacc_vec[i])))
    }
    
    vaccine_set <- c(vaccine_set, rep(doses_per_day[length(doses_per_day)], time_period- days_to_vacc_vec[length(doses_per_day)]))
    
  } else if (length(doses_per_day) ==1)  {
    vaccine_set <- c(vaccine_set, rep(doses_per_day[1], time_period- days_to_vacc_vec[1])) #If there is only one vaccination pace
  }
  
  
  # Output #### 
  
  return(
    list(
      vaccine_set = vaccine_set,
      vaccine_coverage_strategy = vaccine_coverage_strategy,
      next_dose_priority = next_dose_priority))
  
  
  
  
}