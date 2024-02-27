### ------------------------------------------------------------------------ ###
###          Functions needed for the run_function_main.R code
### ------------------------------------------------------------------------ ###


add_cols <- function(x, df1){
  if (!x %in% names(df1)){
    n <- nrow(df1)
    df1[x] <- rep(0,n)
  } else {
    df1 <- df1
  }
  return(df1)
}

# Get representative country for each income group 
get_representative_country <- function(country_type){
  case_when(country_type == "LIC" ~ "Chad",
            country_type == "LMIC" ~ "Pakistan",
            country_type == "UMIC" ~ "Argentina",
            country_type == "HIC" ~ "Canada" )
}


# Get representative mixing matrix for each income group
get_representative_contacts <- function(country_type){
  case_when(country_type == "LIC" ~ "Chad",
            country_type == "LMIC" ~ "Pakistan",
            country_type == "UMIC" ~ "Argentina",
            country_type == "HIC" ~ "Canada" )
}

# Get healthcare capacity (ICU and Hospital beds) based on income group and country specific data 
get_capacity <- function(country, income_group, pop, hs_constraints){
 
   hc <- squire::get_healthcare_capacity(country = country)
  
  # Unconstrained healthcare
  if(hs_constraints == "Absent"){
    hc$hosp_beds <- 1
    hc$ICU_beds <- 1
  }
  
  
  hc$hosp_beds <- round(hc$hosp_beds * sum(pop) / 1000)
  hc$ICU_beds <- round(hc$ICU_beds * sum(pop) / 1000)
  
  return(hc)
}



get_prob_non_severe_death_treatment <- function(income_group, hs_constraints){
  psdt <- squire:::probs$prob_non_severe_death_treatment
  
  if(income_group  == "LIC" & hs_constraints == "Present"){
    psdt <- c(rep(0.25, 16), 0.5804312)
  }
  return(psdt)
}
