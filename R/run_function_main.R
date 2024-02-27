# main running function
run_scenario <-
  function(fit="UKHSA_v6pn_65+_20220702_AZPD2=FALSE_SB=FALSE_NewDecay=TRUE_AddBst=FALSE_AltSev=FALSE_mcmc_chain",
           scenario = 1,
           target_pop = 1e6,
           income_group = "LMIC",
           hs_constraints = "Present",
           vacc_start_vec,
           vaccine_doses,
           vaccine,
           age_groups_covered = 14,
           age_groups_covered_d3 = 14,
           age_groups_covered_d4 = 14,
           age_groups_covered_d5 = 14,
           age_groups_covered_d6 = 14,
           age_groups_covered_d7 = 14,
           age_groups_covered_d8 = 14,
           age_groups_covered_d9 = 14,
           age_groups_covered_d10= 14,
           matrix_d3= NA,
           matrix_booster=NA,
           matrix_primary = NA,
           dt = 0.25,
           repetition = 1,
           seeding_cases = 10,
           vfr,
           vfr2 = 1,
           vfr_time1,
           vfr_time2,
           vfr2_time1 = NA,
           vfr2_time2 = NA,
           max_Rt_var2_scal = NA,
           dose_3_fold_increase = 1,
           dose_4_fold_increase = 1,
           vacc_per_week,
           name = "scenario1",
           std10 = 0.44,
           std10_infection = 0.44,
           t_d2 = 28,
           t_d3=NA,
           t_d4=NA,
           t_d5=NA,
           t_d6=NA,
           t_d7=NA,
           t_d8=NA,
           t_d9=NA,
           t_d10=NA,
           mu_ab_d1,
           mu_ab_d2,
           k,
           hl_s,
           hl_l,
           period_s,
           period_l,
           ab_50,
           ab_50_severe,
           mu_ab_infection,
           mu_ab_inf_scal_vfr = 1,
           max_ab,
           hosp_scal_vfr = 1,
           ICU_scal_vfr = 1,
           hosp_scal_vfr2 = 1,
           ICU_scal_vfr2 = 1,
           variant_specific= 0,
           vp_date = NA,
           vaccine_vfr = 1,
           end_date = "12/31/2025",
           vfr_drift_factor = 1,
           rt_drift_factor = 1,
           ...){
    
    # set up transmission
    R0_t0 <- as.Date(x = "2/1/2020", format = "%m/%d/%Y")
    tmax_date <- as.Date(x = end_date, format = "%m/%d/%Y")
    time_period <- as.integer(difftime(tmax_date, R0_t0 - 1))
    
    rt_out <- generate_Rt_lmic(max_Rt_var2_scal = max_Rt_var2_scal,
                               tmax_date = end_date,
                               vfr_time1 = vfr_time1,
                               vfr_time2 = vfr_time2,
                               vfr2_time1 = vfr2_time1,
                               vfr2_time2 = vfr2_time2,
                               name = name)
    
    # set up daily per-capita prob of external infection with external pulse
    lambda_external <- get_external_foi_lmic(time_period = time_period,
                                             R0_t0 = R0_t0)
    
    
    # implement gradual variant replacement in form of small rt increase every 6 months
    drift_start <- as.Date("2022-04-01", format = "%Y-%m-%d")
    drift_start_day <- as.integer(difftime(drift_start, R0_t0 - 1))
    
    x <- rep(rt_drift_factor, 122)
    for (i in 1:20){
      x <- c(x, rep(tail(x,1)*rt_drift_factor, 122))
    }
    
    rt_drift_multiplier <- c(rep(1, drift_start_day - 1), x)[1:time_period]
    rt_out$Rt <- rt_out$Rt * rt_drift_multiplier
    
    # BASE PARAMETERS
    
    # Population and mixing
    rep_country <- get_representative_country(country_type = income_group)
    iso3c <- countrycode(rep_country, origin = "country.name", destination = "iso3c")
    pop <- squire::get_population(country = rep_country)
    pop_standardise <- target_pop / sum(pop$n)
    pop$n <- as.integer(pop$n * pop_standardise)
    
    contact_country <- get_representative_contacts(country_type = income_group)
    contact_mat <- squire::get_mixing_matrix(country = contact_country)
    
    # Hospital capacity
    hc <- get_capacity(country = rep_country, income_group = income_group, pop = pop$n, hs_constraints = hs_constraints)
    
    # Poorer health outcomes for LMICs and LICs
    pnsdt = get_prob_non_severe_death_treatment(income_group, hs_constraints)
    
    # base parameters
    parameters <- safir::get_parameters(
      population = pop$n,
      contact_matrix_set = contact_mat,
      iso3c = iso3c,
      R0 = rt_out$Rt[1:time_period],
      tt_R0 = rt_out$Rt_tt[1:time_period],
      time_period = time_period,
      dt = dt,
      hosp_bed_capacity = hc$hosp_beds,
      ICU_bed_capacity = hc$ICU_beds,
      prob_non_severe_death_treatment = pnsdt,
      seeding_cases = seeding_cases,
      lambda_external = lambda_external[1:time_period]
    )
    
    # --------------------------------------------------------------------------------
    # Get vaccine and immune parameters
    # --------------------------------------------------------------------------------
    
    # doses available each day
    vacc_per_week <- as.numeric(strsplit(vacc_per_week,",")[[1]])
    
    doses_per_day <- floor(sum(pop$n) * vacc_per_week / 7)
    
    # Day to change vaccination pace
    vacc_start_vec <- strsplit(vacc_start_vec,",")[[1]]
    vacc_start_vec <- as.Date(x = vacc_start_vec, format = "%m/%d/%Y")
    days_to_vacc_vec <- as.integer(difftime(vacc_start_vec, R0_t0))
    
    age_groups_covered_vec <- c(age_groups_covered, age_groups_covered_d3,age_groups_covered_d4,age_groups_covered_d5, 
                                age_groups_covered_d6,  age_groups_covered_d7, age_groups_covered_d8, age_groups_covered_d9,
                                age_groups_covered_d10)
    vaccine_out <-
      get_vaccine_strategy(
        days_to_vacc_vec = days_to_vacc_vec,
        age_groups_covered_vec=age_groups_covered_vec,
        doses_per_day = doses_per_day,
        time_period = time_period,
        vaccine_doses = vaccine_doses,
        vacc_per_week = vacc_per_week,
        matrix_d3= matrix_d3,
        matrix_booster=matrix_booster,
        matrix_primary=matrix_primary)
    
    
    vaccine_set <- vaccine_out$vaccine_set
    vaccine_coverage_strategy <- vaccine_out$vaccine_coverage_strategy
    next_dose_priority <- vaccine_out$next_dose_priority
    
    
    
    # Dosing
    dose_period <- c(NaN,28, na.omit(c(t_d3, t_d4, t_d5, t_d6, t_d7, t_d8, t_d9,t_d10)))
    dose_period <- dose_period[1:vaccine_doses]
    
    
    # make VFR reduction vector and attach
    vfr_time1 <- as.Date(x = vfr_time1, format = "%m/%d/%Y")
    vfr_time2 <- as.Date(x = vfr_time2, format = "%m/%d/%Y")
    vfr2_time1 <- as.Date(x = vfr2_time1, format = "%m/%d/%Y")
    vfr2_time2 <- as.Date(x = vfr2_time2, format = "%m/%d/%Y")
    
    stopifnot(vfr_time1 < tmax_date)
    stopifnot(vfr_time2 < tmax_date)
    stopifnot(vfr2_time1 < tmax_date)
    stopifnot(vfr2_time2 < tmax_date)
    
    
    vfr_time1_day <- as.integer(difftime(vfr_time1, R0_t0 - 1))
    vfr_time2_day <- as.integer(difftime(vfr_time2, R0_t0 - 1))
    vfr2_time1_day <- as.integer(difftime(vfr2_time1, R0_t0 - 1))
    vfr2_time2_day <- as.integer(difftime(vfr2_time2, R0_t0 - 1))
    
    vfr_vector <- c(rep(1, (vfr_time1_day - 1)),
                    seq(1, vfr, length = vfr_time2_day - vfr_time1_day + 1),
                    rep(vfr, (vfr2_time1_day - vfr_time2_day)),
                    seq(vfr, vfr2, length = vfr2_time2_day - vfr2_time1_day + 1),
                    rep(vfr2, time_period - vfr2_time2_day-1))
    
    
    
    # include additional gradual antigenic drift in the form of small change in VFR every 4 months to represent immune escape
    x <- rep(vfr_drift_factor, 122)
    for (i in 1:60){
      x <- c(x, rep(tail(x,1)*vfr_drift_factor, 122))
    }
    vfr_drift_multiplier <- c(rep(1, drift_start_day - 1),x)
    vfr_drift_multiplier <- vfr_drift_multiplier[1:time_period]
    
    
    # set up mu_ab_infection vector
    mu_ab_infection_vector <- c(rep(mu_ab_infection, (vfr_time1_day - 1)),
                                seq(mu_ab_infection, mu_ab_infection*vfr, length = vfr_time2_day - vfr_time1_day + 1),
                                rep(mu_ab_infection*vfr, (vfr2_time1_day - vfr_time2_day)),
                                seq(mu_ab_infection*vfr, mu_ab_infection*vfr2, length = vfr2_time2_day - vfr2_time1_day + 1),
                                rep(mu_ab_infection*vfr2, time_period - vfr2_time2_day-1))
    
    # first way of doing this applied the 50% in protection against reinfection again each time a new variant emerged.
    
    mu_ab_inf_scal_vfr_vector <- c(rep(1, (vfr_time1_day - 1)),
                                   seq(1, mu_ab_inf_scal_vfr, length = vfr_time2_day - vfr_time1_day + 1),
                                   rep(mu_ab_inf_scal_vfr, time_period - vfr_time2_day ))
    
    mu_ab_infection_vector_in <- t(as.matrix(mu_ab_infection_vector * mu_ab_inf_scal_vfr_vector * vfr_drift_multiplier))
    
    
    
    vfr_final <-  vfr_vector * vfr_drift_multiplier
    
    
    ## Correcting vaccine induced immunity based on vaccine type for boosters: 
    #  If variant proof vaccine, mu_ab_d4 and onward s to create vaccine efficacy ~80% against hospitalisarion for all variants. 
    
    time_to_dose <- numeric(length(dose_period))
    time_to_dose[1] <- days_to_vacc_vec[1]
    for (i in 2:length(dose_period)) {
      time_to_dose[i] <- time_to_dose[i - 1] + dose_period[i]
    }
    
    
    mu_ab_d3 <-  mu_ab_d2 * dose_3_fold_increase
    
    if(variant_specific ==0 ){                # If vaccine is not variant specific 
      
      if(!is.na(vp_date)){                    # If vaccine is variant proof 
        mu_ab_d4 = mu_ab_d2 *0.4448183
      }else {mu_ab_d4=mu_ab_d3}
      
      mu_ab_list <- data.frame(name = vaccine,
                               mu_ab_d1 = mu_ab_d1,
                               mu_ab_d2 = mu_ab_d2) %>%
        mutate(mu_ab_d3 = mu_ab_d3) %>%
        mutate(mu_ab_d4 = mu_ab_d4 ) %>%
        mutate(mu_ab_d5 = mu_ab_d4)%>%
        mutate(mu_ab_d6 = mu_ab_d4) %>%
        mutate(mu_ab_d7 = mu_ab_d4) %>%
        mutate(mu_ab_d8 = mu_ab_d4) %>%
        mutate(mu_ab_d9 = mu_ab_d4) %>%
        mutate(mu_ab_d10 = mu_ab_d4)
      
      vax_pars <- safir::get_vaccine_ab_titre_parameters(
        vaccine = vaccine, max_dose = vaccine_doses, correlated = TRUE,
        hl_s = hl_s, hl_l = hl_l, period_s = period_s, t_period_l = period_l,
        ab_50 = ab_50, ab_50_severe = ab_50_severe, std10 = std10, k = k,
        mu_ab_list = mu_ab_list
      )
      
      vax_pars$max_ab <- max_ab # max titre on natural log scale
    }else if (variant_specific ==1 ){
      
      
      mu_ab_list <- data.frame(name = vaccine,
                               mu_ab_d1 = mu_ab_d1,
                               mu_ab_d2 = mu_ab_d2) %>%
        mutate(mu_ab_d3 = mu_ab_d2 * dose_3_fold_increase) %>%
        mutate(mu_ab_d4 = ifelse(is.na(time_to_dose[4]), mu_ab_d3, mu_ab_d3 * vfr_final[(time_to_dose[4] - 365)])) %>%
        mutate(mu_ab_d5 =  ifelse(is.na(time_to_dose[5]), mu_ab_d3, mu_ab_d3 * vfr_final[(time_to_dose[5] - 365)])) %>%
        mutate(mu_ab_d6 =  ifelse(is.na(time_to_dose[6]), mu_ab_d3, mu_ab_d3 * vfr_final[(time_to_dose[6] - 365)])) %>%
        mutate(mu_ab_d7 =  ifelse(is.na(time_to_dose[7]), mu_ab_d3, mu_ab_d3 * vfr_final[(time_to_dose[7] - 365)])) %>%
        mutate(mu_ab_d8 =  ifelse(is.na(time_to_dose[8]), mu_ab_d3, mu_ab_d3 * vfr_final[(time_to_dose[8] - 365)])) %>%
        mutate(mu_ab_d9 =  ifelse(is.na(time_to_dose[9]), mu_ab_d3, mu_ab_d3 * vfr_final[(time_to_dose[9] - 365)])) %>%
        mutate(mu_ab_d10 =  ifelse(is.na(time_to_dose[10]), mu_ab_d3, mu_ab_d3 * vfr_final[(time_to_dose[10] - 365)])) 
      
      
      vax_pars <- safir::get_vaccine_ab_titre_parameters(
        vaccine = vaccine, max_dose = vaccine_doses, correlated = TRUE,
        hl_s = hl_s, hl_l = hl_l, period_s = period_s, t_period_l = period_l,
        ab_50 = ab_50, ab_50_severe = ab_50_severe, std10 = std10, k = k,
        mu_ab_list = mu_ab_list
      )
      
      vax_pars$max_ab <- max_ab # max titre on natural log scale
      
      
    }
    
    
    # Variant proof vaccine 
    if (!is.na(vp_date)){
      vp_date <- as.Date(x = vp_date, format = "%m/%d/%Y")
      vp_time =  as.integer(difftime(vp_date, R0_t0))
    } else {vp_time =-1}
    
    
    
    # combine parameters and verify
    parameters <- safir:: make_vaccine_parameters(
      safir_parameters = parameters,
      vaccine_ab_parameters = vax_pars,
      vaccine_set = vaccine_set,
      dose_period = dose_period,
      strategy_matrix = vaccine_coverage_strategy,
      next_dose_priority_matrix = next_dose_priority,
      vp_time = vp_time
    )
    
    parameters$mu_ab_infection <- mu_ab_infection
    
    # Include immune parameters
    parameters <-
      make_immune_parameters(
        parameters = parameters,
        vfr =vfr_final,
        mu_ab_infection = mu_ab_infection_vector_in,
        std10_infection = std10_infection
      )
    
    # --------------------------------------------------------------------------------
    #  Apply modification to rate of hospitalisation and ICU
    # --------------------------------------------------------------------------------
    
    v <- parameters$prob_hosp
    parameters$prob_hosp <-
      matrix(v,
             nrow = time_period,
             ncol = length(v),
             byrow = TRUE)
    mult <-
      c(
        rep(1, vfr_time1_day - 1),
        seq(1, hosp_scal_vfr, length.out = vfr_time2_day - vfr_time1_day + 1),
        rep(hosp_scal_vfr, (vfr2_time1_day - vfr_time2_day)),
        seq(hosp_scal_vfr, hosp_scal_vfr2, length = vfr2_time2_day - vfr2_time1_day + 1),
        rep(hosp_scal_vfr2, time_period - vfr2_time2_day-1)
      )
    new <- parameters$prob_hosp * mult
    parameters$prob_hosp <- t(new)
    
    v <- parameters$prob_severe
    parameters$prob_severe <-
      matrix(v,
             nrow = time_period,
             ncol = length(v),
             byrow = TRUE)
    mult <-
      c(
        rep(1, vfr_time1_day - 1),
        seq(1, ICU_scal_vfr, length.out = vfr_time2_day - vfr_time1_day + 1),
        rep(ICU_scal_vfr, (vfr2_time1_day - vfr_time2_day)),
        seq(ICU_scal_vfr, ICU_scal_vfr2, length = vfr2_time2_day - vfr2_time1_day + 1),
        rep(ICU_scal_vfr2, time_period - vfr2_time2_day-1)
      )
    new <- parameters$prob_severe * mult
    parameters$prob_severe <- t(new)
    
    # read in decay rate vector
    dr_vec <- readRDS(paste0("data/dr_vec_",fit,"_SD1.rds"))
    dr_vec_vaccine <- dr_vec %>%
      select(-t,-N)
    
    # assume that the decay rate for natural infection is the same as for the vaccine
    dr_vec_inf <- dr_vec %>%
      select(N)
    
    dr_vec_doses_m <- data.matrix(dr_vec_vaccine)
    dr_vec_inf_m <- data.matrix(dr_vec_inf)
    
    parameters <-
      make_independent_vaccine_infection_nat_parameters(
        parameters = parameters,
        dr_vec_doses = dr_vec_doses_m,
        dr_vec_inf = dr_vec_inf_m,
        max_ab_inf = max_ab
      )
    
    ######################################################
    # run the simulation
    ######################################################
    
    # create variables
    timesteps <- parameters$time_period/dt
    
    # creates the categorical states and ages for the simulated population
    variables <- create_variables(pop = pop, parameters = parameters)
    variables <- create_vaccine_variables(variables = variables, parameters = parameters)
    variables <- create_natural_immunity_variables(variables = variables, parameters = parameters)
    variables <- create_independent_nat_variables(variables = variables, parameters = parameters)
    
    # creates the list of events and attaches listeners which handle state changes and queue future events
    events <- create_events(parameters = parameters)
    events <- create_events_vaccination(events = events, parameters = parameters)
    attach_event_listeners(variables = variables, events = events, parameters = parameters, dt = dt)
    attach_event_listeners_vaccination(variables = variables,events = events,parameters = parameters,dt = dt)
    # attach_event_listeners_natural_immunity(variables = variables, events = events, parameters = parameters, dt = dt, additive = TRUE)
    attach_event_listeners_independent_nat(variables = variables, events = events, parameters = parameters, dt = dt)
    
    # renderer object is made
    renderer <- individual::Render$new(parameters$time_period)
    
    age_incidence_renderer <- individual::Render$new(timesteps)
    #attach_tracking_listener_incidence(events=events, renderer = incidence_renderer)
    attach_tracking_listener_age_incidence(
      events = events,
      renderer = age_incidence_renderer,
      age = variables$discrete_age,
      parameters = parameters
    )
    
    
    vaxx_renderer <- individual::Render$new(parameters$time_period)
    inf_renderer <- individual::Render$new(parameters$time_period)
    incidence_renderer <- individual::Render$new(timesteps)
    attach_tracking_listener_incidence(events=events, renderer = incidence_renderer)
    nat_renderer <- individual::Render$new(parameters$time_period)
    nat_inf_renderer <- individual::Render$new(parameters$time_period)
    sp_renderer <- individual::Render$new(parameters$time_period)
    hosp_render <- create_hosp_renderers(parameters = parameters)
    
    # track incidence
    attach_hosp_listeners(renderers = hosp_render, events = events)
    
    double_count_render_process_daily <- function(renderer, variable, dt) {
      stopifnot(inherits(variable, "DoubleVariable"))
      stopifnot(inherits(renderer, "Render"))
      function(t) {
        if ((t * dt) %% 1 == 0) {
          day <- as.integer(t * dt)
          nat <- exp(variable$get_values())
          quantiles_nat <- quantile(x = nat, probs = c(0.025, 0.5, 0.975))
          renderer$render(name = "lower", value = quantiles_nat[[1]], timestep = day)
          renderer$render(name = "median", value = quantiles_nat[[2]], timestep = day)
          renderer$render(name = "upper", value = quantiles_nat[[3]], timestep = day)
          renderer$render(name = "mean", value = mean(x = nat), timestep = day)
        }
      }
    }
    
    sp_render_process_daily <- function(renderer, variable1, variable2, dt) {
      stopifnot(inherits(variable1, "DoubleVariable"))
      stopifnot(inherits(variable2, "DoubleVariable"))
      stopifnot(inherits(renderer, "Render"))
      function(t) {
        if ((t * dt) %% 1 == 0) {
          day <- as.integer(t * dt)
          nat <- exp(variable1$get_values()) + exp(variable2$get_values())
          quantiles <- quantile(x = nat, probs = c(0.025, 0.5, 0.975))
          renderer$render(name = "nat_lower", value = quantiles[[1]], timestep = day)
          renderer$render(name = "nat_median", value = quantiles[[2]], timestep = day)
          renderer$render(name = "nat_upper", value = quantiles[[3]], timestep = day)
          renderer$render(name = "nat_mean", value = mean(x = nat), timestep = day)
          sp <-ifelse(nat>=0.2,1,0)
          quantiles <- quantile(x = sp, probs = c(0.025, 0.5, 0.975))
          renderer$render(name = "sp_lower", value = quantiles[[1]], timestep = day)
          renderer$render(name = "sp_median", value = quantiles[[2]], timestep = day)
          renderer$render(name = "sp_upper", value = quantiles[[3]], timestep = day)
          renderer$render(name = "sp_mean", value = mean(x = sp), timestep = day)
          
        }
      }
    }
    
    # processes
    processes <- list(
      independent_ab_titre_process(parameters=parameters, variables = variables, dt = dt),
      vaccination_process(parameters = parameters,variables = variables,events = events, dt = dt),
      infection_process_vaccine_cpp(parameters = parameters,variables = variables,events = events, dt = dt),
      categorical_count_renderer_process_daily(renderer = renderer, variable = variables$states, categories = variables$states$get_categories(),dt = dt),
      double_count_render_process_daily(renderer = nat_renderer, variable = variables$ab_titre, dt = dt ),
      double_count_render_process_daily(renderer = nat_inf_renderer, variable = variables$ab_titre_inf, dt = dt ),
      sp_render_process_daily(renderer = sp_renderer, variable1 = variables$ab_titre, variable2 = variables$ab_titre_inf, dt = dt),
      integer_count_render_process_daily(renderer = vaxx_renderer, variable = variables$dose_num, margin = 1:vaccine_doses, dt = dt),
      compartments_age_render_process_daily(
        renderer = age_incidence_renderer,
        age = variables$discrete_age,
        compartments = variables$states,
        parameters = parameters,
        dt = dt))
    
    setup_events(parameters = parameters, events=events, variables = variables, dt=dt)
    
    simulation_loop_safir(
      variables = variables,
      events = events,
      processes = processes,
      timesteps = timesteps,
      variables_dont_update = c("discrete_age", "phase"),
      progress = TRUE
    )
    
    
    df <- renderer$to_dataframe()
    df_vacc <- vaxx_renderer$to_dataframe()
    
    
    
    
    df_inc <- incidence_renderer$to_dataframe()
    df_nat <- nat_renderer$to_dataframe() %>%
      rename(vax_ab_mean = mean,
             vax_ab_median = median,
             vax_ab_lower = lower,
             vax_ab_upper = upper)
    
    
    df_nat_inf <- nat_inf_renderer$to_dataframe() %>%
      rename(nat_ab_mean = mean,
             nat_ab_median = median,
             nat_ab_lower = lower,
             nat_ab_upper = upper)
    df_sp <- sp_renderer$to_dataframe()
    
    df_inc <- df_inc %>%
      mutate(timestep = floor((timestep-1)*dt)) %>%
      mutate(incidence = coalesce(incidence,0))%>%
      group_by(timestep) %>%
      summarise(incidence = sum(incidence))
    
    df_hosp <- process_hosp_renderers(renderers = hosp_render, parameters = parameters)
    
    df_hosp <- df_hosp %>%
      mutate(hosp = hosp_get_live + hosp_get_die,
             hosp_all = hosp + hosp_not_get_live + hosp_not_get_die,
             ICU = ICU_get_live + ICU_get_die,
             ICU_all = ICU + ICU_not_get_live + ICU_not_get_die,
             timestep = day-1) %>%
      select(timestep,hosp,hosp_all,ICU,ICU_all)
    
    df <- left_join(df, df_vacc, by = c("timestep"))
    df <- left_join(df, df_inc, by = c("timestep"))
    df <- left_join(df, df_nat, by = c("timestep"))
    df <- left_join(df, df_nat_inf, by = c("timestep"))
    df <- left_join(df, df_sp, by = c("timestep"))
    
    df_rt <- as.data.frame(rt_out) %>%
      rename("timestep" = "Rt_tt")
    df <- left_join(df, df_rt, by = c("timestep"))
    df <- left_join(df, df_hosp, by = c("timestep"))
    
    tmp <- as.data.frame(age_incidence_renderer$to_dataframe())
    
    tmp <- tmp %>%
      mutate(timestep = floor((timestep-1)*dt)) %>%
      select(contains(c("incidence", "_D_age","_ICase_age","_IOx","_IMV")), timestep)
    
    
    # # Plots to tesT:
    # #
    # 
    # df_vacc$timestep <- df_vacc$timestep + R0_t0
    # dose_out <- melt(as.data.table(df_vacc),id.vars="timestep")
    # setnames(dose_out, "variable", "dose")
    # 
    # ggplot(data = dose_out) +
    #   geom_line(aes(x=timestep,y=value,color=dose)) +
    #   theme_bw()
    
    # 
    # df_inc$timestep <- df_inc$timestep + R0_t0
    # 
    # p1 <-  ggplot(data = df_novacc) +
    #   geom_line(aes(x=timestep,y=incidence)) +
    #   theme_bw()
    # p2 <-  ggplot(data = df_novacc_drift) +
    #   geom_line(aes(x=timestep,y=incidence)) +
    #   theme_bw()
    # 
    # p3 <-  ggplot(data = df_inc) +
    #   geom_line(aes(x=timestep,y=incidence)) +
    #   theme_bw()
    
    
    # Save output
    saveRDS(tmp, file=paste0("./raw_outputs/", name,"_age/", "scenario" ,scenario, ".rds"))
    saveRDS(df,file=paste0("./raw_outputs/", name,"/", "scenario" ,scenario, ".rds"))
    
  }
