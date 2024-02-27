library (tidyverse)
library(safir)
library(future)
library(furrr)
library(purrr)
library(dplyr)
library(didehpc)

name <- "LMIC"
fit <- "UKHSA_v6pn_65+_20220702_AZPD2=FALSE_SB=FALSE_NewDecay=TRUE_AddBst=FALSE_AltSev=FALSE_mcmc_chain"
#### Get vaccine parameters  ##############################################
vaccine <- "Moderna"
vacc_names <- data.frame(vaccine = c("Moderna", "Oxford-AstraZeneca"), vacc = c("MD", "AZ"))
vaccine_s <- vaccine
vacc_params <- readRDS(paste0("data/param_list_",fit,".rds")) %>%
  rename(vacc = vaccine) %>%
  left_join(vacc_names, by = "vacc") %>%
  filter(vaccine == vaccine_s) %>%
  mutate(std10 = 0.44) %>%
  filter(vfr > 1) %>%
  select(-c(vacc))

#### Set up other simulation parameters  ##############################################
target_pop <- 1e6
income_group <- "LMIC"
hs_constraints <- "Present"
dt <- 0.25
repetition <-  1:100
vacc_start_vec <- "06/23/2021" 
vacc_per_week <- "0.01"
vaccine_doses <- 6
age_groups_covered <- 16
age_groups_covered_d3 <- 16
age_groups_covered_d4 <- 5
seeding_cases <- 10
matrix_primary<- "./data/LMIC_coverage_primary.csv"
matrix_d3 <-  "./data/LMIC_coverage_primary.csv"
matrix_booster<- "./data/boosters_coverage.csv"
t_d2 <- 28
t_d3 <- 365
t_d4 <- 365
vfr_time1 <- "11/27/2021"
vfr_time2 <- "12/31/2021"
vfr2_time1 <- "10/1/2024" # wont have any effect if vfr2 <- vfr, hosp_scale_vfr <- hosp_scale_vfr2 and ICU_scal_vfr <- ICU_scal_vfr2
vfr2_time2 <- "10/31/2024"
vfr <- sort(unique(vacc_params$vfr))[2]
vfr2 <- c(vfr,7.5,10)
max_Rt_var2_scal <- c(1,1.1,1.3)
ICU_scal_vfr <-  0.3
hosp_scal_vfr <- 0.3
ICU_scal_vfr2 <- c(0.3,0.65,1)
hosp_scal_vfr2 <- c(0.3,0.65,1)
mu_ab_infection <- 1
mu_ab_inf_scal_vfr <- 0.5
max_ab <- 5
variant_specific <- c(0,1)
vp_date<- c(NA, "01/01/2024")
dose_3_fold_increase <- 1.227400
vfr_drift_factor <- 1.025
rt_drift_factor <- 1 
end_date <- "12/31/2026"


#### Create scenarios ##########################################################

scenarios <- expand_grid(fit=fit,
                         matrix_primary =matrix_primary,
                         matrix_d3 =matrix_d3,
                         matrix_booster =matrix_booster,
                         income_group = income_group,
                         target_pop = target_pop,
                         hs_constraints = hs_constraints,
                         vaccine_doses = vaccine_doses,
                         vaccine = vaccine,
                         age_groups_covered = age_groups_covered,
                         age_groups_covered_d3 = age_groups_covered_d3,
                         vacc_start_vec = vacc_start_vec,
                         dt = dt,
                         repetition = repetition,
                         seeding_cases = seeding_cases,
                         t_d3 = t_d3,
                         t_d4= t_d4,
                         vfr = vfr,
                         vfr2 = vfr2,
                         max_Rt_var2_scal = max_Rt_var2_scal,
                         vfr_time1 = vfr_time1,
                         vfr_time2 = vfr_time2,
                         vfr2_time1 = vfr2_time1,
                         vfr2_time2 = vfr2_time2,
                         mu_ab_infection = mu_ab_infection,
                         mu_ab_inf_scal_vfr = mu_ab_inf_scal_vfr,
                         max_ab = max_ab,
                         hosp_scal_vfr = hosp_scal_vfr,
                         ICU_scal_vfr = ICU_scal_vfr,
                         hosp_scal_vfr2 = hosp_scal_vfr2,
                         ICU_scal_vfr2 = ICU_scal_vfr2,
                         variant_specific = variant_specific,
                         vp_date=vp_date,
                         dose_3_fold_increase = dose_3_fold_increase,
                         vfr_drift_factor = vfr_drift_factor,
                         rt_drift_factor = rt_drift_factor,
                         end_date= end_date,
                         vacc_per_week = vacc_per_week
) %>%
  mutate(age_groups_covered_d5 = age_groups_covered_d4,
         age_groups_covered_d6 = age_groups_covered_d4,
         age_groups_covered_d7 = age_groups_covered_d4,
         age_groups_covered_d8 = age_groups_covered_d4,
         age_groups_covered_d9 = age_groups_covered_d4,
         age_groups_covered_d10 = age_groups_covered_d4 ) %>%
  mutate(t_d4 = 365,  t_d5 = 365,t_d6 = 365,  t_d7 = 365, t_d8 = 365,  t_d9 = 365, t_d10=365) %>%
  filter((ICU_scal_vfr2==0.65 & hosp_scal_vfr2==0.65 & max_Rt_var2_scal==1.1 & vfr2==7.5) | ( ICU_scal_vfr2 == 1 & hosp_scal_vfr2==1 & max_Rt_var2_scal==1.3 & vfr2==10)| (ICU_scal_vfr2 == 0.3 & hosp_scal_vfr2==0.3& max_Rt_var2_scal==1 & vfr2==vfr)) %>%
  filter ((is.na(vp_date) & variant_specific ==1) | (!is.na(vp_date) & variant_specific==0))  #Filter new variants scenarios 


scenarios$scenario <- 1:nrow(scenarios)
scenarios$name <- name
scenarios <- left_join(scenarios, vacc_params %>% select(-dose_3_fold_increase), by = c("vaccine", "vfr")) # Remove dose3_fold_increase


nrow(scenarios)

write.csv(scenarios, paste0("./scenarios_", name, ".csv"), row.names = FALSE)
dir.create(paste0("./raw_outputs/", name,"_age/"))
dir.create(paste0("./raw_outputs/", name))

# Checking which scenarios are missing ####
scenarios <- read_csv(paste0("scenarios_", name, ".csv"), show_col_types = FALSE)
#
# # read in raw outputs
df_all <- list.files(path = paste0("./raw_outputs/",name,"/"), pattern = ".rds")
df_all <- map(paste0("./raw_outputs/",name,"/", df_all), readRDS)
scenario_num <- data.frame(scenario_num = list.files(path = paste0("./raw_outputs/",name,"/"), pattern = ".rds")) %>%
  separate(scenario_num, c("A", "scenario_num"), sep = "o") %>%
  separate(scenario_num, c("scenario", "B"), sep = ".rds") %>%
  select(scenario) %>%
  mutate(scenario = as.double(scenario)) # ID of scenarios that have output


scenarios_torun <- data.frame(seq(1:nrow(scenarios)))
colnames(scenarios_torun) <- c("seq")

missing <-  scenarios_torun$seq[!scenarios_torun$seq %in% scenario_num$scenario]

scenarios <- scenarios %>% filter(scenario %in% missing)



# # # ## test on PC
#  source("./R/load_packages.R")
#  source("./R/run_function_main.R")
# source("./R/utils.R")
#  source("./R/vaccine_strategy.R")
#  source("./R/generate_rt.R")
#  source("./R/generate_external_foi.R")
#  plan(multicore, workers = 2)
#  system.time({out <- future_pmap(scenarios[1,], run_scenario, .progress = TRUE)})
#  out <-system.timepmap(scenarios[1,], run_scenario, .progress = TRUE)

#### Run the model on cluster ###############################################
sources <- c("R/run_function_main.R", "R/utils.R", "R/vaccine_strategy.R", "R/generate_rt.R", "R/generate_external_foi.R")
src <- conan::conan_sources(c("github::mrc-ide/safir@vp_vaccine", "mrc-ide/squire", "mrc-ide/nimue"))
ctx <- context::context_save("context",
                             sources = sources,
                             packages = c("tibble", "dplyr", "tidyr", "countrycode", "safir", "nimue", "squire", "data.table"),
                             package_sources = src)


config <- didehpc::didehpc_config(use_rrq = FALSE, use_workers = FALSE, cluster="fi--didemrchnb")

# Create the queue
run <- didehpc::queue_didehpc(ctx, config = config)



t <- run$enqueue(sessionInfo()) # will have packages from context attached
t$status()

t$result() # if I run a function, the output will be returned into the results. only return what you need!
# see cluster load: 
run$cluster_load(TRUE)

# Run
runs2 <- run$enqueue_bulk(scenarios, run_scenario, do_call = TRUE, progress = TRUE)
runs2$status()




