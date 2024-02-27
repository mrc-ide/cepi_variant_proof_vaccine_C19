library(ggplot2)
library(tidyverse)
library(dplyr)
library(scales)
library(RColorBrewer)
library(patchwork)


# Color palettes for Figures 
color_ab <- c('mediumpurple4','maroon4', 'orangered3')
col_set_fig2 <- c("#433E85FF","#1F9A8AFF", "#4DC36BFF","#BBDF27FF","#440154FF","#2E6F8EFF" )

col_country <- c ("#5A5353", "#A07178", "#E6CCBE","#776274", "#C8CC92")

########### 0. Load postprocessed results ##### 

df_LIC <- readRDS(paste0("processed_outputs/df_summarise_dynamics_LIC.rds"))
df_LMIC <- readRDS(paste0("processed_outputs/df_summarise_dynamics_LMIC.rds"))
df_UMIC <- readRDS(paste0("processed_outputs/df_summarise_dynamics_UMIC.rds"))
df_HIC <- readRDS(paste0("processed_outputs/df_summarise_dynamics_HIC.rds"))


df_vac <- rbind(df_LIC, df_LMIC,df_UMIC,df_HIC ) %>%  mutate(atemp = date - as.Date("2023-07-01")) %>%
  mutate(income_group=fct_relevel(income_group,c("LIC","LMIC","UMIC","HIC")))


d4 <- df_vac%>% filter(strategy_name=="baseline,Variant Proof") %>% 
  group_by(income_group) %>% filter(dose4_t >0) %>% 
  slice(1) %>% select(income_group, date) %>% rename(New_vaccine = date) %>% mutate(atemp =New_vaccine - as.Date("2023-07-01") )
d4$atemp <- as.numeric(d4$atemp/365)


# Create a data frame for immune escape annotations
immune_escape <- data.frame(
  label = c("VFR1 (Omicron)", "Drift begins", "New Variant"),  # Labels for immune escape dates
  date = c("2021-11-27", "2022-04-01", "2024-10-1"))  # Immune escape dates

new_variant <- as.numeric(as.Date("2024-10-01")-as.Date("2023-07-01"))/365

########### 1.Dynamics plots #####

to_plot <-  "moderate"  #"baseline" #"worst_case"#

deaths_df <- df_vac %>% 
  select ("income_group", "epi_scenario", "vaccine_type", "deaths_t" , "deaths_tmin","deaths_tmax","date") %>%
  separate(date, into = c("year", "month", "day"), sep="-") %>% drop_na() %>% 
  group_by(epi_scenario,vaccine_type,year,month,income_group) %>%
  summarise(across(deaths_t:deaths_tmax, ~sum(.x, na.rm = TRUE))) %>%
  mutate(date= as.Date(paste0(year,"-",month,"-01")))


plot_deaths <- ggplot(data = deaths_df %>% filter(epi_scenario==to_plot & as.Date(date)>as.Date("2022-08-20"))) + 
  geom_line(aes(x=(as.Date(date)-as.Date("2023-07-01"))/365, y=deaths_t, color=vaccine_type),size=0.8) +
  #geom_ribbon(aes(x=(as.Date(date)-as.Date("2023-07-01"))/365,ymin = deaths_tmin, ymax = deaths_tmax, fill = vaccine_type), alpha = 0.5, col = NA) +
  geom_vline(data=d4, aes(xintercept = atemp),linetype = "dashed", color= "gray30") +
  geom_vline(xintercept= new_variant, color= "#2D4654") +
  facet_grid(~income_group,scales = "free_y") + 
  theme(panel.grid.major = element_blank(),
        strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0)+  lims(x = c(-1,3)) +
  scale_color_manual(values = col_set_fig2) + scale_fill_manual(values = col_set_fig2) +
  #ylim(0,5) +
  theme_bw(base_size = 13) + labs (x="Time (years)", y="Monthly deaths per million",color="Vaccine type")

plot_deaths

hosp_df <- df_vac %>% 
  select ("income_group", "epi_scenario", "vaccine_type", "hosp_t" , "hosp_tmin","hosp_tmax","date") %>%
  separate(date, into = c("year", "month", "day"), sep="-") %>% drop_na() %>% 
  group_by(epi_scenario,vaccine_type,year,month,income_group) %>%
  summarise(across(hosp_t:hosp_tmax, ~sum(.x, na.rm = TRUE))) %>%
  mutate(date= as.Date(paste0(year,"-",month,"-01")))

plot_hosp <- ggplot(data = hosp_df %>% filter(epi_scenario==to_plot& as.Date(date)>as.Date("2022-08-20") )) + 
  geom_line(aes(x=(as.Date(date)-as.Date("2023-07-01"))/365, y=hosp_t, color=vaccine_type),size=0.8) +
  geom_vline(xintercept= new_variant, color= "#2D4654") +
 # geom_line(aes(x=(as.Date(date)-as.Date("2023-07-01"))/365, y=hosp_tmin, color=vaccine_type),size=0.8, linetype=2) +
 # geom_line(aes(x=(as.Date(date)-as.Date("2023-07-01"))/365, y=hosp_tmax, color=vaccine_type),size=0.8, linetype=2) +
  geom_vline(data=d4, aes(xintercept = atemp),linetype = "dashed", color= "gray30") +
  facet_grid(~income_group,scales = "free_y") + 
  theme(panel.grid.major = element_blank(),
        strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0)+#  lims(x = c(-1,3), y=c(0,800)) +
  scale_color_manual(values = col_set_fig2) + theme_bw(base_size = 13) + labs (x="Time (years)", y="Monthly hospitalisations per million",color="Vaccine type")

plot_hosp

# Infection dynamics   

inf_df <- df_vac %>% 
  select ("income_group", "epi_scenario", "vaccine_type", "inc_t" , "inc_tmin","inc_tmax","date") %>%
  separate(date, into = c("year", "month", "day"), sep="-") %>% drop_na() %>% 
  group_by(epi_scenario,vaccine_type,year,month,income_group) %>%
  summarise(across(inc_t:inc_tmax, ~sum(.x, na.rm = TRUE))) %>%
  mutate(date= as.Date(paste0(year,"-",month,"-01")))


plot_inf <- ggplot(data = inf_df ) + #%>% filter(as.Date(date)>as.Date("2022-08-20"))
  geom_line(aes(x=(as.Date(date)-as.Date("2023-07-01"))/365, y=inc_t, color=vaccine_type),size=0.8) +
  #geom_ribbon(aes(x=(as.Date(date)-as.Date("2023-07-01"))/365,ymin = deaths_tmin, ymax = deaths_tmax, fill = vaccine_type), alpha = 0.5, col = NA) +
  facet_grid(epi_scenario~income_group,scales = "free_y") + 
  theme(panel.grid.major = element_blank(),
        strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0)+  #lims(x = c(-1,3), y=c(0,50)) +
  scale_color_manual(values = col_set_fig2) + scale_fill_manual(values = col_set_fig2) +
  theme_bw(base_size = 13) + labs (x="Time (years)", y="Monthly infections per million",color="Vaccine type")  


########### 2. Summary error bars #####

df_all <- read.csv("./processed_outputs/table_output.csv") %>%  mutate(income_group=fct_relevel(income_group,c("LIC","LMIC","UMIC","HIC")))

deahts_bar <- ggplot(data = df_all %>% filter(epi_scenario ==to_plot), aes(x = income_group, y = deaths_med, fill = income_group)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single"))+
  geom_text(aes(label = paste0(round(deaths_med,2),"")), size = 4,  vjust = -2, , hjust=-0.3, position = position_dodge(1), color="black") +
  geom_errorbar( aes(x=income_group, ymin=deaths_lower, ymax=deaths_upper), colour="black",width=0, size=0.8)+
  #facet_grid(~epi_scenario,scales = "free_y") +
  theme_bw(base_size = 13) +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_blank(),
        legend.text.align = 0) +
  #scale_fill_brewer(palette = "Dark2") +
  scale_fill_manual(values = col_country) +  scale_y_continuous(labels = scales::comma)+
  labs(x = "Country category", y = "Deaths averted \n per million individuals", fill = "Country category") 

hosp_bar <- ggplot(data = df_all %>% filter(epi_scenario ==to_plot), aes(x = income_group, y = hosp_med, fill = income_group)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single"))+
  geom_text(aes(label = paste0(round(hosp_med,2),"")), size = 4,  vjust =- 1.3, , hjust=-0.3, position = position_dodge(1), color="black") +
  geom_errorbar( aes(x=income_group, ymin=hosp_lower, ymax=hosp_upper), colour="black",width=0, size=0.8)+
  #facet_grid(~epi_scenario,scales = "free_y") +
  theme_bw(base_size = 13) +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_blank(),
        legend.text.align = 0) +
  #scale_fill_brewer(palette = "Dark2") +
  scale_fill_manual(values = col_country) +  scale_y_continuous(labels = scales::comma)+
  labs(x = "Country category", y = "Hospitalisations averted \n per million individuals", fill = "Country category") 


layout <- "
AA
BB
CD
"
combined_epi <- plot_deaths + 
  plot_hosp + 
  deahts_bar + 
  hosp_bar + 
  plot_annotation(tag_levels = "A") + 
  plot_layout(guides = "collect", design = layout)

combined_epi


# Economic plots '#### 
vsl_bar <- ggplot(data = df_all %>% filter(epi_scenario =="moderate"), aes(x = income_group, y = vsl/1e6, fill = income_group)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single"))+
  geom_text(aes(label = paste0(round(vsl/1e6,2),"")), size = 4,  vjust = 1, , hjust=1, position = position_dodge(1), color="black") +
  geom_errorbar( aes(x=income_group, ymin=vsl_low/1e6, ymax=vsl_upper/1e6), colour="black",width=0, size=0.8)+
  #facet_grid(~epi_scenario,scales = "free_y") +
  theme_bw(base_size = 13) +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_blank(),
        legend.text.align = 0) +
  #scale_fill_brewer(palette = "Dark2") +
  scale_fill_manual(values = col_country) +  scale_y_continuous(labels = scales::comma)+
  labs(x = "Country category", y = "VSL averted per million individuals \n(Million GBP)", fill = "Country category") 

vsly_bar <- ggplot(data = df_all %>% filter(epi_scenario =="moderate"), aes(x = income_group, y = vsly/1e6, fill = income_group)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single"))+
  geom_text(aes(label = paste0(round(vsly/1e6,2),"")), size = 4,  vjust = 1, , hjust=1, position = position_dodge(1), color="black") +
  geom_errorbar( aes(x=income_group, ymin=vsly_low/1e6, ymax=vsly_upper/1e6), colour="black",width=0, size=0.8)+
  #facet_grid(~epi_scenario,scales = "free_y") +
  theme_bw(base_size = 13) +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_blank(),
        legend.text.align = 0) +
  #scale_fill_brewer(palette = "Dark2") +
  scale_fill_manual(values = col_country) +  scale_y_continuous(labels = scales::comma)+
  labs(x = "Country category", y = "VSLY averted per million individuals \n(Million GBP)", fill = "Country category") 


prod_loss_bar <- ggplot(data = df_all %>% filter(epi_scenario =="moderate"), aes(x = income_group, y = prod_loss/1e6, fill = income_group)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single"))+
  geom_text(aes(label = paste0(round(prod_loss/1e6,2),"")), size = 4,  vjust = 1, , hjust=1, position = position_dodge(1), color="black") +
  geom_errorbar( aes(x=income_group, ymin=prod_loss_low/1e6, ymax=prod_loss_upper/1e6), colour="black",width=0, size=0.8)+
  #facet_grid(~epi_scenario,scales = "free_y") +
  theme_bw(base_size = 13) +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_blank(),
        legend.text.align = 0) +
  #scale_fill_brewer(palette = "Dark2") +
  scale_fill_manual(values = col_country) +  scale_y_continuous(labels = scales::comma)+
  labs(x = "Country category", y = "Productivity loss per million individuals \n(Million GBP)", fill = "Country category") 


hosp_cost_bar <- ggplot(data = df_all %>% filter(epi_scenario =="moderate"), aes(x = income_group, y = hosp_cost/1e6, fill = income_group)) +
  geom_col(position = position_dodge2(width = 0.9, preserve = "single"))+
  geom_text(aes(label = paste0(round(hosp_cost/1e6,2),"")), size = 4,  vjust = 1, , hjust=1, position = position_dodge(1), color="black") +
  geom_errorbar( aes(x=income_group, ymin=hosp_cost_low/1e6, ymax=hosp_cost_upper/1e6), colour="black",width=0, size=0.8)+
  #facet_grid(~epi_scenario,scales = "free_y") +
  theme_bw(base_size = 13) +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        #axis.text.x = element_blank(),
        legend.text.align = 0) +
  #scale_fill_brewer(palette = "Dark2") +
  scale_fill_manual(values = col_country) +  scale_y_continuous(labels = scales::comma)+
  labs(x = "Country category", y = "Hospitalisations costs averted per million individuals \n(Million GBP)", fill = "Country category") 

layout2 <- "
AB
CD
"
combined_eco <- vsl_bar + 
  vsly_bar + 
  prod_loss_bar + 
  hosp_cost_bar + 
  plot_annotation(tag_levels = "A") + 
  plot_layout(guides = "collect", design = layout2)

combined_eco


########### 3. Blue-green doses barplot #####
df_vac$income_group <- factor(df_vac$income_group, levels =c("LIC","LMIC", "UMIC","HIC"))
df_vac$atemp <- as.numeric(df_vac$atemp/365)

df_vac <- df_vac  %>% filter(epi_scenario =="baseline") %>%
  rename("Dose 1" = "dose1_t", "Dose 2" = "dose2_t", "Booster 1" = "dose3_t","Booster 2" = "dose4_t" , "Booster 3" = "dose5_t","Booster 4" = "dose6_t",  "Booster 5" = "dose7_t") %>% 
  select(-"Booster 5") %>%  pivot_longer(cols = c("Dose 1", "Dose 2", "Booster 1","Booster 2","Booster 3","Booster 4"), names_to = "dose") %>%
  mutate(dose = factor(dose, levels = c("Dose 1", "Dose 2", "Booster 1","Booster 2","Booster 3","Booster 4"), ordered = TRUE)) 

df1_doses_month <- df_vac %>%
  # filter to last date of each month
  mutate(year = lubridate::year(date),
         month = lubridate::month(date),
         day = lubridate::day(date),
         date = lubridate::floor_date(date, "month")) %>%
  group_by(income_group, target_pop, vaccine_doses, dose, year, month) %>% 
  mutate(max_day = max(day)) %>%
  ungroup() %>%
  filter(day == max_day)

plot_doses <- ggplot(data = df1_doses_month %>%  filter(vaccine_type== "Variant Proof"), aes(x = (as.Date(date)-as.Date("2023-07-01"))/365, y = value/target_pop*100, fill = dose)) +
  geom_bar(stat = "identity",width=0.5) +
  geom_vline(data=d4, aes(xintercept = atemp),linetype = "dashed", color= "gray30") +
  facet_grid(~income_group, scales="free_y") +
  theme_bw(base_size = 16) +
  theme(strip.background = element_rect(fill = NA),
        panel.border = element_blank(),
        axis.line = element_line(),
        legend.text.align = 0) + #axis.text.x = element_text(angle = 335, vjust = 0.3, hjust=0.2)
  #lims(x= c(-1,3)) +
  scale_fill_brewer(palette = "PuBuGn",direction = -1) +
  labs(x = "Time (years)", y = "Vaccinated (%)", fill = "Dose number") 

plot_doses



