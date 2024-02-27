library(purrr)
library(tidyverse)
library(dplyr)

name <- "UMIC"

post_process <- function(name){
  
  # get a few parameters we need - specify start of time period for when you want to calculate totals
  R0_t0 <- as.Date(x = "1/2/2020", format = "%d/%m/%Y")
  start_sum <- as.Date(x = "2/2/2020", format = "%d/%m/%Y")
  days_to_start_sum <- as.integer(difftime(start_sum, R0_t0))
  
  # read in scenarios
  scenarios <- read_csv(paste0("scenarios_", name, ".csv"), show_col_types = FALSE)
  
  # read in raw outputs
  df_all <- list.files(path = paste0("./raw_outputs/",name,"/"), pattern = ".rds")
  df_all <- map(paste0("./raw_outputs/",name,"/", df_all), readRDS)
  scenario_num <- data.frame(scenario_num = list.files(path = paste0("./raw_outputs/",name,"/"), pattern = ".rds")) %>%
    separate(scenario_num, c("A", "scenario_num"), sep = "o") %>%
    separate(scenario_num, c("scenario", "B"), sep = ".rds") %>%
    select(scenario) %>%
    mutate(scenario = as.double(scenario))
  

  
  # 1. Summarise each run ####
  for (i in 1:length(df_all)){
   
     saf_reps_summarise <- df_all[[i]] %>%
      mutate(IMild_count = IMild_count + IAsymp_count) %>%
      dplyr::select(-IAsymp_count) %>%
      pivot_longer(cols = contains(c("count", "Rt","incidence","nat","sp", "vax","hosp","ICU")), names_to = "compartment") %>%
      filter(compartment %in% c("D_count", "X1_count", "X2_count", "X3_count", "X4_count", "X5_count", "X6_count", "X7_count", "X8_count", "X9_count","R_count","E_count", "IMild_count", "ICase_count", "Rt", "incidence","vax_ab_median", "vax_ab_lower", "vax_ab_upper", "nat_ab_median", "nat_ab_lower", "nat_ab_upper", "nat_median", "nat_lower","nat_upper","sp_median","sp_lower", "sp_upper", "hosp", "hosp_all", "ICU", "ICU_all" )) %>%
      group_by(compartment) %>%
      mutate(value = if_else(compartment == "D_count", value - lag(value), value),
             value = if_else(is.na(value), 0, value)) %>%
      ungroup() %>%
      pivot_wider(id_cols = timestep, names_from = "compartment", values_from = "value")  %>%
      mutate(deaths = sum(D_count[timestep >= days_to_start_sum]),
             cum_hosp = sum(hosp[timestep >= days_to_start_sum]),
             cum_hosp_all = sum(hosp_all[timestep >= days_to_start_sum]),
             cum_ICU = sum(ICU[timestep >= days_to_start_sum]),
             cum_ICU_all = sum(ICU_all[timestep >= days_to_start_sum]),
             inc = sum(incidence[timestep >= days_to_start_sum])) %>%
      ungroup()
    
    cols <- c(X1_count = 0, X2_count = 0, X3_count = 0,
              X4_count = 0, X5_count = 0, X6_count = 0,
              X7_count = 0, X8_count = 0, X9_count = 0)
    
    saf_reps_summarise<- add_column(saf_reps_summarise, !!!cols[setdiff(names(cols), names(saf_reps_summarise))])
    
    saf_reps_summarise <- saf_reps_summarise %>%
      mutate(doses = X1_count + X2_count * 2 + X3_count * 3 + X4_count * 4 + X5_count * 5 + X6_count * 6 + X7_count * 7 + X8_count * 8 + X9_count * 9,
             total_doses = max(doses)) %>%
      nest(cols = c(timestep, D_count, R_count, E_count, IMild_count, ICase_count, Rt, X1_count, X2_count, X3_count, X4_count, X5_count, X6_count, X7_count, X8_count, X9_count, incidence, vax_ab_median, vax_ab_lower, vax_ab_upper, nat_ab_median, nat_ab_lower, nat_ab_upper, nat_median, nat_lower, nat_upper, sp_median, sp_lower, sp_upper, hosp, hosp_all, ICU, ICU_all, doses)) 
    
    saf_reps_summarise$scenario <- scenario_num$scenario[i]
    df_all[[i]] <- saf_reps_summarise
  }
  
  # join the runs and link to parameters
  df_all <- do.call(rbind,df_all)
  
  df <- left_join(df_all, scenarios, by = "scenario") %>% mutate(age_groups_covered_d4=5)
  
  
  # name the options
  df <- df %>%
    mutate(epi_scenario = "undefined",
           vaccine_type = "undefined",
    ) %>%
    mutate(
      epi_scenario = case_when(hosp_scal_vfr2 == 0.3 ~ "baseline",
                               hosp_scal_vfr2 == 0.65 ~ "moderate",
                               hosp_scal_vfr2 == 1 ~ "worst_case",
                               TRUE ~ epi_scenario),
      vaccine_type = case_when(variant_specific == 1 ~ "Variant Specific",
                               variant_specific == 0 ~ "Variant Proof",
                               TRUE ~ vaccine_type)) %>%
    mutate(strategy_name = paste0(epi_scenario, ",",vaccine_type))
  
  # 2. Summarise temporal dynamics over repetitions ####
  df_summarise <- df %>%
    unnest(cols) %>%
    select(-c(deaths, cum_hosp, cum_hosp_all, cum_ICU, cum_ICU_all, total_doses, inc)) %>%
    group_by(timestep,
             income_group,
             target_pop,
             age_groups_covered,
             age_groups_covered_d3,
             age_groups_covered_d4,
             age_groups_covered_d5, 
             age_groups_covered_d6,
             age_groups_covered_d7,
             age_groups_covered_d8,
             age_groups_covered_d9,
             vaccine_doses,
             vacc_start_vec,
             strategy_name,
             epi_scenario,
             vaccine_type,
             vacc_per_week,
             t_d3, t_d4, t_d5, t_d6, t_d7, t_d8, t_d9,
             vfr,
             vfr2,
             vfr_time1,
             vfr_time2,
             vfr2_time1,
             vfr2_time2,
             mu_ab_infection,
             hosp_scal_vfr,
             ICU_scal_vfr,
             hosp_scal_vfr2,
             ICU_scal_vfr2,
             variant_specific,
             vfr_drift_factor,
             rt_drift_factor
    ) %>%
    summarise(deaths_t = median(D_count),
              deaths_tmin = quantile(D_count, 0.025),
              deaths_tmax = quantile(D_count, 0.975),
              hosp_t = median(hosp),
              hosp_tmin = quantile(hosp, 0.025),
              hosp_tmax = quantile(hosp, 0.975),
              ICU_t = median(ICU),
              ICU_tmin = quantile(ICU, 0.025),
              ICU_tmax = quantile(ICU, 0.975),
              E_t = median(E_count),
              IMild_t = median(IMild_count),
              ICase_t = median(ICase_count),
              inc_t = median(incidence),
              inc_tmin = quantile(incidence, 0.025),
              inc_tmax = quantile(incidence, 0.975),
              vax_ab_med = median(vax_ab_median),
              vax_ab_lower = median(vax_ab_lower),
              vax_ab_upper = median(vax_ab_upper),
              nat_ab_med = median(nat_ab_median),
              nat_ab_lower = median(nat_ab_lower),
              nat_ab_upper = median(nat_ab_upper),
              nat_med = median(nat_median),
              nat_lower = median(nat_lower),
              nat_upper = median(nat_upper),
              sp_med = median(sp_median),
              sp_lower = median(sp_lower),
              sp_upper = median(sp_upper),
              vaccines_t = median(X1_count + X2_count * 2 + X3_count * 3 + X4_count * 4 
                                  + X5_count * 5 + X6_count * 6 + X7_count * 7 +
                                    X8_count * 8 + X9_count * 9),
              dose1_t = median(X1_count),
              dose2_t = median(X2_count),
              dose3_t = median(X3_count),
              dose4_t = median(X4_count),
              dose5_t = median(X5_count),
              dose6_t = median(X6_count),
              dose7_t = median(X7_count),
              dose8_t = median(X8_count),
              dose9_t = median(X9_count),
              Rt = median(Rt),
              .groups = 'drop') %>%
    unique() %>%
    mutate(date = timestep + as.Date("2020-02-01"))
  
  
  saveRDS(df_summarise, file = paste0("./processed_outputs/df_summarise_dynamics_", name, ".rds"))
  
  
  
  # 3. Estimate averted per epi_scenatio  ####
  
  # Separating counterfactual (Variant proof) from Variant specific vaccine
  v_proof <- df %>% select(c("deaths", "cum_hosp", "cum_hosp_all","cum_ICU","cum_ICU_all","inc","total_doses", "scenario",
                             "income_group","repetition","vfr2","hosp_scal_vfr2","ICU_scal_vfr2", "variant_specific","epi_scenario",
                             "vaccine_type","strategy_name" )) %>% filter(vaccine_type == "Variant Proof") %>% group_by(epi_scenario)
  
  counter <- df %>% select(c("deaths", "cum_hosp", "cum_hosp_all","cum_ICU","cum_ICU_all","inc","total_doses", "scenario",
                                "income_group","repetition","vfr2","hosp_scal_vfr2","ICU_scal_vfr2", "variant_specific","epi_scenario",
                                "vaccine_type","strategy_name" )) %>% filter(vaccine_type == "Variant Specific") %>% group_by(epi_scenario)
  
  
  summary_averted <- counter %>% left_join(v_proof,by = c("income_group","vfr2","hosp_scal_vfr2","ICU_scal_vfr2","epi_scenario")) %>%
    group_by(epi_scenario) %>% sample_n(size = 200) %>% 
    mutate(deaths_averted = deaths.x- deaths.y,
           hops_averted= cum_hosp_all.x-cum_hosp_all.y,
           cases_averted = inc.x-inc.y,
           icu_averted= cum_ICU_all.x - cum_ICU_all.y) %>% 
    mutate(deaths_med = median(deaths_averted),
           deaths_lower = quantile(deaths_averted, 0.025),
           deaths_upper = quantile(deaths_averted, 0.975),
           hosp_med = median(hops_averted),
           hosp_lower = quantile(hops_averted, 0.025),
           hosp_upper = quantile(hops_averted, 0.975),
           icu_med = median(icu_averted),
           icu_lower = quantile(icu_averted, 0.025),
           icu_upper = quantile(icu_averted, 0.975),
           inc_med = median(cases_averted),
           inc_lower = quantile(cases_averted, 0.025),
           inc_upper = quantile(cases_averted, 0.975))%>%
    select(c(income_group,vfr2,hosp_scal_vfr2,ICU_scal_vfr2,epi_scenario,deaths_med,deaths_lower,deaths_upper,
             hosp_med,hosp_lower,hosp_upper,icu_med,icu_lower,icu_upper,inc_med,inc_lower,inc_upper))%>%
    unique()
  
  
  # 4. YLL analysis ####
  life_expectancy <- read.csv("./data/life_expectancy.csv") %>%  select(age_group,name) %>% rename (expectancy = name, age_band=age_group) %>% mutate(age_group = row_number())
  
  # read in scenarios
  scenarios2 <- read_csv(paste0("scenarios_", name, ".csv"), show_col_types = FALSE) %>% mutate(age_groups_covered_d4=5) %>%
    select(c("age_groups_covered_d4", "vfr2","max_Rt_var2_scal", "hosp_scal_vfr2", "variant_specific", "vfr_drift_factor", "t_d5","name", "repetition", "vaccine_doses","scenario","matrix_booster"))
  
  # read in raw outputs
  df_all_yll <- list.files(path = paste0("./raw_outputs/",name,"_age/"), pattern = ".rds")
  df_all_yll <- map(paste0("./raw_outputs/",name,"_age/", df_all_yll), readRDS)
  scenario_num <- data.frame(scenario_num = list.files(path = paste0("./raw_outputs/",name,"/"), pattern = ".rds")) %>%
    separate(scenario_num, c("A", "scenario_num"), sep = "o") %>%
    separate(scenario_num, c("scenario", "B"), sep = ".rds") %>%
    select(scenario) %>%
    mutate(scenario = as.double(scenario))
  
  total_deaths <- data.frame(matrix(0,nrow=0,ncol=19))
  
  # Get the last row of the output as it has deaths on the last time step
  for (i in 1:length(df_all_yll)){
    
    temp <- df_all_yll[[i]]
    tmp <- temp  %>%select(contains("_D_age"), timestep)  %>% drop_na() %>%
      filter(row_number()==n() | row_number() == days_to_start_sum)
    
    deaths_time_selected <- tmp[2,] - tmp[1,]
    deaths_time_selected$scenario <- scenario_num[i,1]
    
    total_deaths <- rbind(total_deaths,deaths_time_selected)
  }
  
  
  # Joing the deaths per age group to the scenario 
  df <- left_join(total_deaths, scenarios2, by = "scenario") 
  
  
  # name the options
  df <- df %>%
    mutate(epi_scenario = "undefined",
           vaccine_type = "undefined",
    ) %>%
    mutate(
      epi_scenario = case_when(hosp_scal_vfr2 == 0.3 ~ "baseline",
                               hosp_scal_vfr2 == 0.65 ~ "moderate",
                               hosp_scal_vfr2 == 1 ~ "worst_case",
                               TRUE ~ epi_scenario),
      vaccine_type = case_when(variant_specific == 1 ~ "Variant Specific",
                               variant_specific == 0 ~ "Variant Proof",
                               TRUE ~ vaccine_type)) %>%
    mutate(strategy_name = paste0(epi_scenario, ",",vaccine_type))
  
  
  counter_yll <- df %>% filter(vaccine_type == "Variant Specific") %>% select(!c("timestep", "age_groups_covered_d4", "max_Rt_var2_scal", "vfr_drift_factor","t_d5",
                                                                              "vaccine_doses", "matrix_booster")) %>%
    pivot_longer(cols = contains("compartment_D_age"),names_to = c("age_group","stat"), names_pattern = "compartment_D_age_?(.*)_(.*)") %>%
    group_by(vfr2,hosp_scal_vfr2,epi_scenario,name,age_group)
  
  v_yll <- df %>% filter(vaccine_type == "Variant Proof")  %>% select(!c("timestep", "age_groups_covered_d4", "max_Rt_var2_scal", "vfr_drift_factor","t_d5",
                                                                            "vaccine_doses", "matrix_booster")) %>%
    pivot_longer(cols = contains("compartment_D_age"),names_to = c("age_group","stat"), names_pattern = "compartment_D_age_?(.*)_(.*)") %>%
    group_by(vfr2,hosp_scal_vfr2,epi_scenario,name,age_group)
  
  
  df_yll <- v_yll %>% left_join(counter_yll,suffix = c("_vp","_counter"),by= c("vfr2","hosp_scal_vfr2","epi_scenario","name","age_group"))%>%
    group_by (epi_scenario,age_group) %>% sample_n(size = 200) %>%
    mutate(deaths_averted = value_counter-value_vp) %>%
    mutate(deaths_med = median(deaths_averted),
           deaths_lower = quantile(deaths_averted, 0.025),
           deaths_upper = quantile(deaths_averted, 0.975),
           age_group = as.numeric(age_group)) %>%
    select (name,vfr2,hosp_scal_vfr2,epi_scenario,age_group,deaths_med,deaths_lower,deaths_upper) %>% unique() %>% ungroup() %>%
    group_by(epi_scenario) %>% left_join(life_expectancy, by = "age_group") %>%
    mutate( YLL_median = expectancy*deaths_med,
            YLL_lower= expectancy*deaths_lower,
            YLL_higher= expectancy*deaths_upper)
  
  
  
  
  # 5. Economic analysis ####
  
  median_hospital_days <-  7.5
  prod_ages <- c("15-19","20-24","25-29","30-34","35-39","40-44","45-49","50-54","55-59","60-64")
  # Loading economic parameters 
  econ <- read.csv("./data/vsly.csv") %>% filter (income_group==name)
  
  
  prod_YLL <- df_yll %>%  filter(age_band %in% prod_ages )%>% group_by(epi_scenario) %>%
    summarise(YLL_prod_lower = sum(YLL_lower),
              YLL_prod_med = sum(YLL_median),
              YLL_prod_upper =sum(YLL_higher))
  
  all_YLL <- df_yll %>% group_by(epi_scenario) %>% summarise(YLL_lower = sum(YLL_lower),
                                                             YLL_med=sum(YLL_median),
                                                             YLL_upper=sum(YLL_higher)) %>% 
    left_join(prod_YLL)
  
  
  summary_averted <- summary_averted %>% left_join(all_YLL, by= "epi_scenario") %>%
    mutate( vsl = deaths_med*econ$vsl,
            vsl_low= deaths_lower*econ$vsl,
            vsl_upper =deaths_upper*econ$vsl,
            vsly = YLL_med * econ$vsly,
            vsly_low = YLL_lower * econ$vsly,
            vsly_upper = YLL_upper * econ$vsly,
            prod_loss =YLL_prod_med *econ$gni,
            prod_loss_low =YLL_prod_lower *econ$gni,
            prod_loss_upper =YLL_prod_upper *econ$gni,
            hosp_cost = hosp_med*median_hospital_days*econ$coi_1,
            hosp_cost_low = hosp_lower*median_hospital_days*econ$coi_1,
            hosp_cost_upper = hosp_upper*median_hospital_days*econ$coi_1)
  
  
  return(summary_averted)
  
  
}




averted_LIC <- post_process("LIC")
averted_LMIC <- post_process("LMIC")
averted_UMIC <- post_process("UMIC")
averted_HIC <- post_process("HIC")

#averted_UMIC <- post_process("UMIC_AG_new")
#df_all <- averted_UMIC
df_all <- rbind(averted_LIC,averted_LMIC,averted_UMIC,averted_HIC)
mill <- 1e6


table <- df_all %>%
  mutate(deaths_averted_str = paste0(round(deaths_med,digits = 1), " (",round(deaths_lower,digits = 1),"-",round(deaths_upper,digits = 1),")"),
         hosp_averted_str =paste0(round(hosp_med,digits = 1), " (",round(hosp_lower,digits = 1),"-",round(hosp_upper,digits = 1),")"),
         cases_averted = paste0(round(inc_med,digits = 1), " (",round(inc_lower,digits = 1),"-",round(inc_upper,digits = 1),")"),
         vsl_str = paste0(round(vsl/mill,digits = 1), " (",round(vsl_low/mill,digits = 1),"-",round(vsl_upper/mill,1),")"),
         vsly_str = paste0(round(vsly/mill,digits = 1), " (",round(vsly_low/mill,digits = 1),"-",round(vsly_upper/mill,digits = 1),")"),
         prod_loss_str = paste0(round(prod_loss/mill,digits = 1), " (",round(prod_loss_low/mill,digits=1),"-",round(prod_loss_upper/mill,digits = 1),")"),
         coi = paste0(round(hosp_cost,digits = 1), " (",round(hosp_cost_low,digits = 1),"-",round(hosp_cost_upper,digits = 1),")")) 

write.csv(table,file = "./processed_outputs/table_output.csv")

