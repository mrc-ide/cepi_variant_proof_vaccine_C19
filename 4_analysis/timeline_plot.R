library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(lubridate)


####  0. Load postprocessed outputs ####  
df_LIC <- readRDS(paste0("processed_outputs/df_summarise_dynamics_LIC.rds")) 
df_LMIC <- readRDS(paste0("processed_outputs/df_summarise_dynamics_LMIC.rds")) 
df_UMIC <- readRDS(paste0("processed_outputs/df_summarise_dynamics_UMIC.rds"))
df_HIC <- readRDS(paste0("processed_outputs/df_summarise_dynamics_HIC.rds"))

# Create big data frame 
df_all <- rbind(df_LIC,df_LMIC,df_UMIC,df_HIC)
df_all$epi_scenario <- as.factor(df_all$epi_scenario)
levels(df_all$epi_scenario) <- c("Continous drifting","Moderate", "Wost case")



####  1. Create a data frame for vaccine dates ####
d1 <- df_all %>% filter(strategy_name=="baseline,Variant Proof")%>%
  group_by(income_group) %>% filter(dose1_t >0) %>% 
  slice(1) %>% select(income_group, date) %>% rename(Primary = date)


d3 <- df_all %>% filter(strategy_name=="baseline,Variant Proof")%>%
  group_by(income_group) %>% filter(dose3_t >0) %>% 
  slice(1) %>% select(income_group, date) %>% rename(Booster = date)

d4 <- df_all%>% filter(strategy_name=="baseline,Variant Proof") %>% 
  group_by(income_group) %>% filter(dose4_t >0) %>% 
  slice(1) %>% select(income_group, date) %>% rename(New_vaccine = date)

vax <- left_join(d1,d3, by="income_group") %>% left_join(d4, by="income_group")

#### 2.  Create a data frame for immune escape annotations ####
immune_escape <- data.frame(
  label = c("VFR1 (Omicron)", "Drift begins", "New Variant"),  # Labels for immune escape dates
  date = c("2021-11-27", "2022-04-01", "2024-10-1"))  # Immune escape dates

immune_escape$date <- as.Date(immune_escape$date)


# 3. Function that create Rt plot with annotations for one income group #### 

plot_rt_annotated <- function(cat){
  
  # Select income group
  df<- df_all %>% filter(income_group==cat)%>% select(date,Rt,income_group,epi_scenario) %>%  rename(Scenario=epi_scenario)
  
  
  vacc_period <- vax %>% filter(income_group == cat) %>% 
    pivot_longer(cols= 2:4, names_to= "Vaccine", values_to = "date_start")  %>% mutate (date_end = lead(date_start))
  vacc_period$date_end[3] <- "2027-01-01"

  # Create the plot
  p <- ggplot(df) + geom_line(data= df, aes(x = date, y = Rt, linetype=Scenario), size =1) +
    labs(x = "Time (Years) ", y = "Transmission (Rt)",title = cat) +
    # Add dose annotations
    # geom_segment(
    #   data = vacc_period,
    #   aes(x = date_start, xend = date_start, y = 10.5, yend = 9.5),
    #   arrow = arrow(length = unit(0.25, "cm")),
    #   size=1,
    #   color = "#183A37"
    # ) +
    geom_text(
      data = vacc_period,
      aes(x = date_start, y = 8, label = Vaccine),
      vjust = 1.5,
      size = 3.7,
      position = position_nudge(x = c(-100 ,150, 250)),
      color="#183A37"
    ) +
    geom_vline(data=vacc_period,aes(xintercept = date_start), linetype = "dashed", color= "gray65") +
    # Add immune escape annotations
    geom_segment(
      data = immune_escape,
      aes(x = date, xend = date, y = 1.5, yend = 0.01),
      arrow = arrow(length = unit(0.25, "cm")),
      size=1,
      color = "#231123"
    ) +
    geom_text(
      data = immune_escape,
      aes(x = date, y = 0.2, label = label),
      vjust = 0,
      color = "#231123",
      position = position_nudge(x = c(-280,250, 250)),
      size = 3.5
    ) + 
    # Add Vaccine period bars
    geom_rect(data = vacc_period,
              aes(xmin = as.Date(date_start), xmax = as.Date(date_end), ymin = 9, ymax = 9.5, fill = Vaccine),
              color = NA,
              alpha = 0.7) + scale_fill_manual(values = c("#A69F98", "#8C6057","#B07C9E")) +
    scale_color_manual(values = c("#A28F9D", "#558C8C")) +
    # Customize the appearance of the plot
    theme_bw(base_size = 15) +  # Use a white background
    theme(panel.grid.major = element_blank(),
          strip.background = element_rect(fill = NA),
          panel.border = element_blank(),
          axis.line = element_line(),
          plot.margin = margin(b = 1, t = 5, l = 10, r = 50))+ 
    scale_x_continuous(breaks = seq(ymd("2020-07-01"),ymd("2027-07-01"),"years"),
      labels = as.character(seq(-3,4))) 
 
  return(p) 
  
  
}

#### 4. Aggregating plots for all income categories ####

p1 <- plot_rt_annotated("LIC")
p2 <- plot_rt_annotated("LMIC")
p3 <- plot_rt_annotated("UMIC")
p4 <- plot_rt_annotated("HIC")



library(patchwork)
layout <- "
AB
CD
"
combined <- p1+ p2 + p3 + p4 +
  plot_annotation(tag_levels = "A") + 
  plot_layout(guides = "collect", design = layout) &  theme(legend.position = 'bottom')

combined
ggsave("./plots/Rt_annotated.png", combined, height=10, width=8)

