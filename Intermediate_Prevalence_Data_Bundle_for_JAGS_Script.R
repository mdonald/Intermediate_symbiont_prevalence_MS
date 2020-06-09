## Intermediate symbiont prevalence manuscript
## Bundling data for parameterization of the JAGS model

rm(list=ls())
## Required packages
library(tidyverse)
library(dplyr)
library(readxl)

# either load the workspace which has all of the data files already read in - both options come from Data_prep_intermediate_prev_MS.R
#load("AGHY_data_clean.RData") ## this has all of the data for bundling the data for the JAGS model

## or read in the exported datasets
# Read in merged data from the data prep script -----------------------------------------------------
AGHY.plots <- read.csv("AGHY.plots_check.csv")
AGHY.merge <- read.csv("AGHY_merge.csv") 
AGHY.demo.final <- read.csv("AGHY_demo_final.csv")
AGHY.linreg<-read_excel("AGHY_SFAEF_life_history_expt.xlsx", 
                        sheet="Linear Regression")
recruitment_raw <- read_excel("AGHY_SFAEF_life_history_expt.xlsx", sheet = "Subplot-level data")



# Bundle data for population prevalence estimate --------------------------
AGHY_subplot_dat <- AGHY.merge%>%
  select(newplot,year_t,subplot,water,transmission,target_init_freq,E_plus_liberal,total)%>% 
  mutate(water = as.factor(case_when(water == "Add" ~ "Irrigated",
                                     water == "Control" ~ "Ambient"))) %>% 
  na.omit()%>%
  arrange(year_t,newplot,subplot) %>% 
  dplyr::rename(initial_prev = target_init_freq)



AGHY_plot_dat <- AGHY.merge %>%
  select(newplot, water, target_init_freq) %>% 
  distinct()%>% 
  mutate(water = case_when(water == "Add" ~ "Irrigated",
                           water == "Control" ~ "Ambient")) %>% 
  mutate(water = as.integer(as.factor(water))) %>% 
  arrange(newplot) %>% 
  dplyr::rename(initial_prev = target_init_freq)

## levels
N.trt<-length(levels(AGHY_subplot_dat$water))
N.years<-length(unique(AGHY_subplot_dat$year_t))
N.plots<-length(unique(AGHY_subplot_dat$newplot))
N.obs<-nrow(AGHY_subplot_dat)

## data - predictor variables at the plot level
water<-AGHY_plot_dat$water
initial_prev<-AGHY_plot_dat$initial_prev

## data - response variable at subplot level
y.pos<-AGHY_subplot_dat$E_plus_liberal
N.samples<-AGHY_subplot_dat$total
plot<-AGHY_subplot_dat$newplot
year<-AGHY_subplot_dat$year_t-2013

## data for prediction
x.levels<-seq(0,1,0.01)
N.x.levels<-length(x.levels)


# Bundle data for survival vital rate estimate ----------------------------

## survival data for individuals of known endo status
survival_dat <- AGHY.demo.final %>%
  filter(spring_survival_t1 != "NA") %>%
  select(newplot, subplot, ID, year_t, year_t1, spring_survival_t1, endo.stat) %>%
  left_join(AGHY.plots %>% 
              select(newplot, water)) 

## select the plants have have BOTH survival and endo data
survival_known <- survival_dat %>%
  filter(!is.na(endo.stat))%>%
  mutate(year = year_t -(min(year_t)-1),
         plot = newplot - (min(newplot)-1),
         endo = as.integer(endo.stat)+1,
         water = as.integer(as.factor(water))) %>%
  select(year, plot, endo, water, spring_survival_t1, subplot, ID)%>%
  filter(!is.na(spring_survival_t1))

## select the plants for which there is survival info but we don't know their endo stat
## select the unknowns
survival_unk <- survival_dat %>%
  filter(is.na(endo.stat))%>%
  mutate(year = year_t - (min(year_t)-1),
         plot = newplot - (min(newplot)-1),
         water = as.integer(as.factor(water))) %>%
  select(year, plot, endo.stat, water, spring_survival_t1,subplot, ID)%>%
  filter(!is.na(spring_survival_t1))%>%
  unite(unique_ID, plot, subplot,year, sep = "_", remove = F) 


## combine to the subplot level so that we have number of individuals that survived and total alive last year
## these are just individuals for which we don't know their endo stat  
survival_unk_sbplt <- survival_unk %>% 
  group_by(unique_ID) %>% 
  summarize(alive.t1 = sum(spring_survival_t1),
            alive.t = length(unique_ID),
            water = unique(water)) %>% 
  separate(unique_ID, c("plot", "subplot", "year"), sep = "_", remove = F) %>% 
  mutate(year = (as.numeric(year))) 


## adjust the population prev

survival_known_plot <- survival_dat %>%
  filter(!is.na(endo.stat))%>%
  mutate(year = year_t -(min(year_t)-1),
         plot = newplot - (min(newplot)-1),
         endo = as.integer(endo.stat)+1,
         water = as.integer(as.factor(water))) %>%
  select(year, plot, endo, water, spring_survival_t1, subplot, ID, endo.stat)%>%
  filter(!is.na(spring_survival_t1)) %>% 
  group_by(year,plot,water) %>% 
  summarize(total_known = length(spring_survival_t1),
            total_known_ep = sum(endo.stat),
            total_known_em = total_known - total_known_ep)

## select the plants for which there is survival info but we don't know their endo stat
## select the unknowns
survival_unk_plot <- survival_dat %>%
  filter(is.na(endo.stat))%>%
    mutate(year = year_t - (min(year_t)-1),
         plot = newplot - (min(newplot)-1),
         water = as.integer(as.factor(water))) %>%
  select(year, plot, endo.stat, water, spring_survival_t1, subplot, ID)%>%
  filter(!is.na(spring_survival_t1))%>%
  unite(unique_ID, plot, subplot,year, sep = "_", remove = F) %>% 
  group_by(year,plot,water) %>% 
  summarize(total_unknown = length(spring_survival_t1))

survival_plot <- survival_known_plot %>% 
  full_join(survival_unk_plot) %>% 
  replace_na(list(total_known = 0, total_unknown = 0, total_known_ep = 0, total_known_em = 0)) %>% 
  group_by(year, plot, water) %>% 
  mutate(total_demo_plants = total_known + total_unknown) %>% 
  filter(total_unknown != 0) %>% 
  arrange(plot)



# Bundle data for fertility vital rate estimate ------------------------------------

seedmass_dat <- left_join(
  AGHY.demo.final %>% 
    filter(inf_number_t > 0 & inf_collected_t > 0) %>% 
    select(newplot, subplot, ID, year_t, inf_number_t, inf_collected_t, seed_mass_t,endo.stat),
  AGHY.plots %>% 
    select(newplot,water),
  by="newplot")  

## there are 4 instances for which we have seed masses but the number of infs collected > infs observed
## to avoid making judgement calls - we're going to drop these from the seed mass analysis
seedmass_dat_issues <-  seedmass_dat %>% 
  filter(inf_collected_t > inf_number_t)


seedmass_dat <- seedmass_dat %>% 
  anti_join(seedmass_dat_issues) ## drop these rows for which the number of infs collected is greater than observed
  
##  individuals for which we have seed mass but no endo score
seedmass_unk <- seedmass_dat %>% 
  filter(is.na(endo.stat)) %>% 
  mutate(year = year_t - (min(year_t)-1),
         plot = newplot - (min(newplot)-1),
         water = as.integer(as.factor(water))) %>% 
  select(year,plot,endo.stat,water,seed_mass_t,inf_collected_t) %>% 
  filter(!is.na(seed_mass_t)) 

seedmass_dat <- seedmass_dat%>% 
  mutate(year = year_t - (min(year_t)-1),
         plot = newplot - (min(newplot)-1),
         endo = as.integer(endo.stat+1),
         water = as.integer(as.factor(water))) %>% 
  select(year,plot,endo,water,seed_mass_t,inf_collected_t) %>% 
  na.omit()

### Infs data

infs_dat <- left_join(
  AGHY.demo.final %>% 
    #filter(inf_number_t > 0) %>% 
    filter(!is.na(inf_number_t)) %>% 
    select(newplot, subplot, ID, year_t, inf_number_t, inf_collected_t,seed_mass_t,endo.stat),
  AGHY.plots %>% 
    select(newplot,water),
  by="newplot") 


## find the "issues" - when the number of infs collected > number of infs observed -- there are 8
## In some cases, this was due to the observation occurring ~ 1 week before the infs were collected - so it's possible that another inf emerged during that week
## and we did indeed collect that many 
## For the other 6 instances, the raw data still disagree - to avoid making judgement calls and changing the data, we're going to drop these 8
issues_inf <- infs_dat %>% 
  filter(inf_collected_t > inf_number_t)## there are 8 instances when collected infs > infs observed


## add in the mean prev to the infs_dat
infs_dat_pop_prev <- infs_dat %>% 
  anti_join(issues_inf) %>% ## drop the 8 instances when collected infs > infs observed
  mutate(flowering_t = if_else(inf_number_t >= 1, 1, 0)) ## indicate whether the plant flowered (1) or not (0)

## data for adjusting population prevalence for inf production 

total_plants_infs <- infs_dat_pop_prev %>% 
  group_by(year_t, newplot) %>% 
  select(inf_number_t, inf_collected_t, year_t, newplot, endo.stat, subplot, ID, flowering_t) %>% 
  mutate(Ep_flow = case_when(endo.stat == 1 & inf_number_t >= 1 ~ 1),
         Em_flow = case_when(endo.stat == 0 & inf_number_t >= 1 ~ 0)) %>% 
  summarize(total_plants_flow = length(inf_number_t), flowering_plants = sum(flowering_t), 
            Ep_known_flow = sum(endo.stat, na.rm = T),
            Ep_flow = sum(Ep_flow, na.rm = T),
            Em_flow = sum(Em_flow, na.rm = T),
            total_flowering_plants_known = Ep_flow + Em_flow,
            total_flowering_plants_unknown = total_plants_flow - total_flowering_plants_known) %>% 
  ungroup() %>% 
  mutate(year_t = ifelse(year_t == 2014, 1 | year_t == 2015, 2)) %>% 
  filter(total_flowering_plants_unknown != 0)


## pull out the data for which we don't have endo scores
inf.unk <- infs_dat_pop_prev %>% 
  group_by(year_t, newplot) %>% 
  filter(is.na(endo.stat)) %>% 
  summarize(total_unknown_infs = length(endo.stat),
            unk_flow = sum(flowering_t))

inf.known <- infs_dat_pop_prev %>% 
  group_by(year_t, newplot) %>% 
  filter(!is.na(endo.stat)) %>% 
  summarize(total_known_infs = length(endo.stat),
            known_flow = sum(flowering_t))


## flowering individuals for which we do not know their endo stat
infs_unk <- infs_dat_pop_prev %>% 
  filter(is.na(endo.stat),
         inf_number_t >= 1) %>% ## select only plants that flowered 
  mutate(year = year_t - (min(year_t)-1),
         plot = newplot - (min(newplot)-1),
         endo = 9999,
         water = as.integer(as.factor(water))) %>% 
  select(year,plot,endo,water,inf_number_t) 

infs_dat <- infs_dat_pop_prev %>% 
  filter(inf_number_t >= 1) %>% 
  mutate(year = year_t - (min(year_t)-1),
         plot = newplot - (min(newplot)-1),
         endo = as.integer(endo.stat+1),
         water = as.integer(as.factor(water))) %>% 
  select(year,plot,endo,water,inf_number_t) %>% 
  na.omit()


## because we want the contribution of eps.infs.overdisp in the weighted mean
## we need to make the vector length start at the end of the N.infs.unk.vector
## an easy way to do this is to stack the dfs on top of one another
## then tell the for loop to start runing at the term that starts the unk df

infs_unk_stacked <- infs_dat %>% 
  bind_rows(infs_unk) 



# Bundle data for flowering vital rate estimate ---------------------------

flowering<- AGHY.merge %>% 
  select(year_t, water, newplot, subplot, spring_flowering_recruit_count_per_subplot,spring_vegetative_recruit_count_per_subplot) %>% 
  filter(!is.na(spring_vegetative_recruit_count_per_subplot), !is.na(spring_flowering_recruit_count_per_subplot)) %>%  ## drop subplots for which there is an NA for either flowering or veg. count
  mutate(total.plants = spring_flowering_recruit_count_per_subplot + spring_vegetative_recruit_count_per_subplot ) %>% 
  filter(year_t != 2016) ## DROP 2016 flowering since it doesn't contribute to a paired vital rate and will be easier to index everything in the combined model 


## prep the columns for the JAGS model
flowering_dat <- flowering %>%
  mutate(year = year_t -(min(year_t)-1),
         plot = newplot - (min(newplot)-1),
         water = as.integer(as.factor(water)))

## flowering from demo individuals 
flowering_individuals <- AGHY.demo.final %>% 
  select(newplot,water, subplot, ID, spring_survival_t1, tiller_number_t, inf_number_t, year_t, year_t1, tiller_number_t1, inf_number_t1, endo.stat) %>% 
  mutate(flow_t = if_else(inf_number_t > 0, 1, 0)) %>% 
  filter(endo.stat == 1 | endo.stat == 0,
         flow_t == 1 | flow_t == 0, 
         spring_survival_t1 == 1) %>% 
  mutate(year = year_t -(min(year_t)-1),
         plot = newplot - (min(newplot)-1),
         endo = as.integer(endo.stat)+1,
         water = as.integer(as.factor(water)),
         N.samples = 1)

flowering_sbplts <- flowering_individuals %>%  ## select the individuals for which we know endo stat and have info on whether or not they were flowering, and that were alive in t1
  group_by(newplot, subplot, year_t, flow_t)%>%
  summarize(n(), count = length(newplot),
            num_ep = sum(endo.stat),
            total_known = length(endo.stat),
            endo_flowering =  sum(flow_t[endo != "NA"])) %>% 
  ## endo_flowering is the number of plants of known endo stat that are flowering, typically this is equal to the number of known endo individuals since we get endo from seeds 
  full_join(flowering_dat)

df2_14 <- flowering_sbplts %>% 
  filter(year == 1)

df_15 <- flowering_sbplts %>% 
  filter(year == 2) %>% ## flowering and vegetative recruit counts were switched in 2015 -- fix this here.
  dplyr::rename(spring_flowering_recruit_count_per_subplot = spring_vegetative_recruit_count_per_subplot,
                spring_vegetative_recruit_count_per_subplot = spring_flowering_recruit_count_per_subplot)

flowering_sbplts2 <- df2_14 %>% 
  bind_rows(df_15)

## check for disagreement between number of demography individuals and subplot counts
issues_flow <- flowering_sbplts %>% 
  mutate(issue = ifelse(count > total.plants, 1, 0)) %>%  ## There's one issue where the subplot counts are less than the number of demo plants (plot 1)
filter(issue == 1)

flowering_sbplts = flowering_sbplts2

## going in and making the subplot count for plot 1, subplot 1 in 2014 4 instead of 3
flowering_sbplts <- flowering_sbplts %>% 
  mutate(total.plants = ifelse(newplot == 1 & subplot == 1 & year_t == 2014, 4, total.plants),
         spring_flowering_recruit_count_per_subplot = ifelse(newplot == 1 & subplot == 1 & year_t == 2014, 4, spring_flowering_recruit_count_per_subplot)) %>% 
  ## convert the NAs for total_known to be 0s to make one large df
  mutate(total_known = ifelse(is.na(total_known), 0, total_known),
         total_unknown = total.plants - total_known,
         total_unknown_flowered = spring_flowering_recruit_count_per_subplot - total_known,
         total_unknown_flowered = ifelse(total_unknown_flowered < 0, 0, total_unknown_flowered)) %>% 
  ungroup() %>% 
  mutate(year_t = ifelse(year_t == 2014, 1 | year_t == 2015, 2)) ## change year to be 1 for 2014 and 2 for 2015

## flowering_sbplts df has the information for E+ and E- known individuals that were alive in a given census year and have info on whether or not they flowered
## this df also has the subplot totals for all plants in each subplot when we don't have info on the individuals

## Condense this to the plot level for adjusting pop prev
flowering_plts <- flowering_sbplts %>% 
  group_by(newplot, year, water, plot) %>% 
  replace_na(list(total_known = 0, num_ep = 0, total_known_f = 0, total_known_ep_f = 0, total_plants_f = 0, endo_flowering = 0)) %>% 
  summarize(total_known_f = sum(total_known),
            total_known_ep_f = sum(num_ep),
            total_known_em_f = total_known_f - total_known_ep_f,
            total_plants_f = sum(total.plants),
            total_unknown_f = total_plants_f - total_known_f, 
            total_flowering_f = sum(spring_flowering_recruit_count_per_subplot),
            total_known_flowering_f = sum(endo_flowering),
            total_unknown_flowering_f = total_flowering_f - total_known_flowering_f)

# Bundle data for seed mass to seed count conversion ----------------------
linreg_dat <- AGHY.linreg %>% 
  dplyr::rename(plot = Plot) %>% 
  left_join(AGHY.plots %>% 
              select(c(plot, newplot, water))) %>% 
  mutate(water = as.integer(as.factor(water)),
         endo = Endo + 1)



# Bundle data for recruitment vital rate estimate -------------------------
plot_info_r <- AGHY_plot_dat %>% 
  mutate(plot = newplot) %>% 
  arrange(plot) %>% 
  distinct()

rec_final <- recruitment_raw %>% 
  mutate(plot = as.numeric(as.factor(plot)),
         year = as.numeric(as.factor(year))) %>% 
  mutate(total_plants_in_sbplt = spring_flowering_recruit_count_per_subplot + spring_vegetative_recruit_count_per_subplot) %>% 
  filter(total_plants_in_sbplt != "NA") %>%  ## drop rows for which we don't have an observation of counts
  ungroup() %>% 
  dplyr::left_join(plot_info_r, by = "plot") %>% 
  arrange(year, water) 


rec_plot_final <- rec_final %>% 
  group_by(plot,water,year) %>% 
  summarize(total_plants_sbplt = sum(total_plants_in_sbplt)) %>% 
  arrange(year, water)

### df that's the correct length for the rec estimate year t+1
rec_plot_final_yr1 <- rec_plot_final %>% 
  filter(year == 2 |
           year == 3)

rec_plot_final_yr <- rec_plot_final %>% 
  filter(year == 1 | 
           year == 2)

survival_dat2 <- survival_dat %>% ## get the plot, water, and year to link survival with recruitment for the estimate
  select(newplot, year_t, year_t1, water) %>%
  distinct() %>% 
  mutate(water = as.factor(case_when(water == "Add" ~ "Irrigated",
                                     water == "Control" ~ "Ambient"))) %>% 
  mutate(water = as.numeric(as.factor(water)),
         year = case_when(year_t1 == 2015 ~ 1,
                          year_t1 == 2016 ~ 2)) %>% 
  arrange(year_t1, water)



# Bundle data for vertical transmission estimate --------------------------
vtrans_dat <- AGHY.demo.final %>% 
  filter(year_t == 2014 | year_t == 2015, endo.stat == 1) %>% 
  select(plot, subplot, ID, year_t1, year_t, complete_seeds_scored, total_E_minus, total_E_plus, 
         e_plus_1, total.add.seeds, water, newplot, tot_seeds_scored_t, e_p_score_t, endo.stat) %>% 
  mutate(total_seeds_scored = case_when(complete_seeds_scored != "NA" ~ complete_seeds_scored,
                                        tot_seeds_scored_t != "NA" ~ tot_seeds_scored_t),
         total_ep_seeds = case_when(e_p_score_t != "NA" ~ e_p_score_t,
                                    total_E_plus != "NA" ~ total_E_plus)) %>% 
  select(plot, subplot, ID, water, newplot, year_t, total_seeds_scored, total_ep_seeds) %>% 
  mutate(water = case_when(water == "Add" ~ "Irrigated",
                           water == "Control" ~ "Ambient"),
         water = as.numeric(as.factor(water)),
         year = as.numeric(as.factor(year_t)))

## quick check to get an estimate of what the model outcome should be similar to
vtrans_check <- vtrans_dat %>% 
  dplyr::group_by(year, water) %>% 
  summarize(total_seeds_plt = sum(total_seeds_scored),
            total_ep_plt = sum(total_ep_seeds)) %>% 
  group_by(year,water) %>% 
  mutate(vtrans = total_ep_plt /total_seeds_plt)

#save.image("Data_for_JAGS_Model.RData") # save this as a data object to call in JAGS_parameterization_RScript.R
