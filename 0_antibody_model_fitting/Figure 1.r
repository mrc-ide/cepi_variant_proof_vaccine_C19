library(rstudioapi)
library(tidyverse)
library(patchwork)
library(dplyr)
library(drjacoby)

setwd(dirname(getActiveDocumentContext()$path)) 

source("R/create_profile.R")

##################################
##### LOAD THE PARAMETERS AND DATA
##################################

## main fits
load("Model Fits/UKHSA_v6pn_65+_20220702_AZPD2=FALSE_SB=FALSE_NewDecay=TRUE_AddBst=FALSE_AltSev=FALSE_mcmc_chain.Rdata")  

chain <- sample_chains(mcmc, 10000)

posterior_median <- chain %>%
  summarise( 
  across(where(is.numeric), median)
  )

name <- c("Original","Variant-Specific Drift","Variant-Specific Moderate","Variant-Specific Worst Case", "Variant-Proof") 
vaccine_number <-5


log10_d2_PF <- posterior_median$d2_PF
log10_d2_AZ <- posterior_median$d2_AZ
log10_d2_MD <- posterior_median$d2_MD


d1_AZ       <- log10_d2_AZ + posterior_median$d1_AZ
d1_PF       <- log10_d2_PF + posterior_median$d1_PF
d1_MD       <- log10_d2_MD + posterior_median$d1_MD

d2_AZ    <- posterior_median$d2_AZ
d2_PF    <- posterior_median$d2_PF
d2_MD    <- posterior_median$d2_MD

d3_AZ    <- posterior_median$bst_AZ
d3_PF    <- posterior_median$bst_PF
d3_MD    <- posterior_median$bst_MD

ab_50       <- posterior_median$ni50 
ab_50_severe <- posterior_median$ns50
ab_50_death  <- posterior_median$nd50

k           <- posterior_median$k
hl_s        <- posterior_median$hl_s
hl_l        <- posterior_median$hl_l
period_s    <- posterior_median$period_s

# fixed parameter
std10 <- 0.44 # Pooled standard deviation of antibody level on log10 scale

om_red <- posterior_median$om_red
ratio <-1.85
bivalent_scaling <- (1/ratio)
om_red <- om_red*bivalent_scaling

vs_reduction_drift <- log10(5.815209) - log10(6.262344) 
vs_reduction_moderate <- log10(5.815209) - log10(9.138022)
vs_reduction_worst <- log10(5.815209) - log10(12.184029)
vp_reduction <- log10(0.4448183)

mu_ab_d1 <- c(d1_MD, d1_MD, d1_MD, d1_MD, d1_MD)
mu_ab_d2 <- c(d2_MD, d2_MD, d2_MD, d2_MD, d2_MD)
mu_ab_d3 <- c(d3_MD, d3_MD, d3_MD, d3_MD, d3_MD)
mu_ab_d4 <- c(d3_MD, d3_MD + vs_reduction_drift, d3_MD + vs_reduction_moderate, d3_MD + vs_reduction_worst, d2_MD + vp_reduction)   


# transforms
dr_s <- -log(2)/hl_s  # Corresponding decay rate in days for half life above
dr_l <- -log(2)/hl_l

# Timing of doses
max_t <- 365*4 # number of days to model
t <- 0:max_t
t_d2 <- 84 # timing of second dose relative to first
t_d3 <- 180 # timing of third dose relative to second dose
t_d4 <- 365 # timing of fourth dose relative to third dose

param_list <-
  data.frame(
    name,
    mu_ab_d1,
    mu_ab_d2,
    mu_ab_d3,
    mu_ab_d4,
    t_d2,
    t_d3,
    t_d4,
    dr_s,
    dr_l,
    period_s,
    ab_50,
    ab_50_severe,
    ab_50_death
  ) 
  

# initialise other parameters


r1_summary <- NULL
summary_stats <- NULL

##########################
###### OUTPUT OPTIONS
##########################

#####################################################################################
### loop through the vaccines to calculate the profiles of NAT and efficacy over time
#####################################################################################

for (j in 1:vaccine_number){
    
    mu_ab_d1 <- param_list$mu_ab_d1[j] 
    mu_ab_d2 <- param_list$mu_ab_d2[j] 
    mu_ab_d3 <- param_list$mu_ab_d3[j]
    mu_ab_d4 <- param_list$mu_ab_d4[j]
    dr_s <- param_list$dr_s[j]
    dr_l <- param_list$dr_l[j]
    period_s <- param_list$period_s[j]
    t_d2 <- param_list$t_d2[j]
    t_d3 <- param_list$t_d3[j]
    t_d4 <- param_list$t_d4[j]
    ab_50 <- param_list$ab_50[j]
    ab_50_severe <- param_list$ab_50_severe[j]
    ab_50_death <- param_list$ab_50_death[j]
    
    out0 <-
      draw(
        mu_ab_d1 = mu_ab_d1,
        mu_ab_d2 = mu_ab_d2,
        mu_ab_d3 = mu_ab_d3,
        mu_ab_d4 = mu_ab_d4,
        std10 = 0,
        ab_50 = ab_50,
        ab_50_severe = ab_50_severe,
        ab_50_death = ab_50_death,
        dr_s = dr_s,
        dr_l = dr_l,
        period_s = period_s,
        t = t,
        t_d2 = t_d2,
        t_d3 = t_d3,
        t_d4 = t_d4,
        k = k
      )
    
    sub0 <-
      data.frame(
        t = t,
        run = 0,
        Titre = out0$titre,
        Efficacy = out0$VE,
        Efficacy_Severe = out0$VE_severe,
        Death = out0$VE_death
      )
    
    sub0 <- sub0 %>%
      mutate(Efficacy = Efficacy * 100,
             Efficacy_Severe = Efficacy_Severe * 100,
             Death = Death * 100) %>%
      pivot_longer(cols = c("Titre", "Efficacy", "Efficacy_Severe", "Death"), names_to = "type") %>%
      mutate(type = factor(type, levels = c("Titre", "Efficacy", "Efficacy_Severe", "Death"))) %>%
      mutate(j=j)
    
    summary_stats <-rbind(summary_stats,sub0)
}

df <- summary_stats %>%
  mutate(type = factor(type, levels = c("Titre", "Efficacy", "Efficacy_Severe", "Death"), labels = c("IL", "mild disease", "hospitalisation", "death"))) 

data_to_output <- df %>%
  mutate(vaccines = case_when(j==1 ~ "Original",
                              j==2 ~ "Variant-Specific Drift",
                              j==3 ~ "Variant-Specific Moderate",
                              j==4 ~ "Variant-Specific Worst Case",
                              j==5 ~ "Variant-Proof")
  ) %>%
  filter(t<=2*365)
#saveRDS(data_to_output, "../Figures/vaccine_profiles.rds")


initial_efficacy <- data_to_output %>%
  filter((t==30)|(t==60)|(t==90)|(t==180)|(t==365))


plots <- NULL

plot_efficacy_function <- function(d, i){
  plot_out <- ggplot(data = filter(d, type != "IL"), aes(x = t, y = value, col = factor(vaccines))) +
    geom_line() +
#    geom_point(aes(x = t, y=VE), size = 1.5) +
 #   geom_errorbar(aes(x = t, ymin=L95, ymax=U95), width=15) +
    facet_wrap(~ type, ncol = 3) +
    lims(y = c(0,100)) +
    theme_bw() +
    theme(strip.background = element_rect(fill = NA),
          panel.border = element_blank(),
          axis.line = element_line(),
          legend.text.align = 0) +
    scale_x_continuous(breaks = c(0, 365, 365*2, 365*3, 365*4)) +
    scale_y_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100))
    labs(x = "time (days)", y = "efficacy (%)", col = "variant")
  return(plot_out)
}

plot_nat_function <- function(d, i){
  plot_out <- ggplot(data = filter(d, type == "IL"), aes(x = t, y = value, col = factor(variant))) +
    geom_line() +
    geom_point(aes(x = t, y=VE), size = 1.5) +
    geom_errorbar(aes(x = t, ymin=L95, ymax=U95), width=15) +
    facet_wrap(~ type, ncol = 3) +
    theme_bw() +
    theme(strip.background = element_rect(fill = NA),
          panel.border = element_blank(),
          axis.line = element_line(),
          legend.text.align = 0) +
    scale_x_continuous(breaks = c(0, 365, 365*2, 365*3, 365*4)) +
    scale_y_log10(limits = c(1e-3,1e1)) +
    labs(x = "time (days)", y = "IL", col = "variant", title = name[i])
  return(plot_out)
}

plot_efficacy_function2 <- function(d, i){
  plot_out <- ggplot(data = filter(d, type != "IL"), aes(x = t, y = value, col = factor(variant))) +
    geom_line() +
    geom_point(aes(x = t, y=VE), size = 1.5) +
    geom_errorbar(aes(x = t, ymin=L95, ymax=U95), width=15) +
    facet_wrap(~ type, ncol = 3) +
    lims(y = c(0,100)) +
    theme_bw() +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(),
          axis.title.x=element_blank(),
          legend.text.align = 0) +
    scale_x_continuous(breaks = c(0, 365, 365*2, 365*3, 365*4)) +
    labs(x = "time (days)", y = "effectiveness (%)", col = "variant")
  return(plot_out)
}

plot_nat_function2 <- function(d, i){
  plot_out <- ggplot(data = filter(d, type == "IL"), aes(x = t, y = value, col = factor(variant))) +
    geom_line() +
    geom_point(aes(x = t, y=VE), size = 1.5) +
    geom_errorbar(aes(x = t, ymin=L95, ymax=U95), width=15) +
    facet_wrap(~ type, ncol = 3) +
    theme_bw() +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(),
          axis.title.x=element_blank(),
          legend.text.align = 0) +
    scale_x_continuous(breaks = c(0, 365, 365*2, 365*3, 365*4)) +
    scale_y_log10(limits = c(1e-3,1e1)) +
    labs(x = "time (days)", y = "IL", col = "variant")
  return(plot_out)
}


test_plot <- plot_efficacy_function(data_to_output,1)
test_plot

for (i in 1:3){
  d1 <- filter(df, j==i )
  plots[[i*2]] <- plot_efficacy_function2(d1, i)
  plots[[i*2-1]] <- plot_nat_function2(d1, i)
}
