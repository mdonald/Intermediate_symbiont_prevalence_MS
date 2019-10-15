## Title: AGHY Full Bayes model, prevalence and vital rates figures
## Purpose: Use the full distribution from pop prev estimates to estimate vital rates
## Notes: Adjusted pop prev works for survival, flowering, fertility and recruitment! DONE?!?! - just need to confirm that adj pop prev is correct
## Date Started: eons ago at this point
## Date Updated: 19 Sept 2019
## Updated recruitment model to the new version in the MS


rm(list=ls())
## Required packages
library(R2jags);load.module("mix")
library(mcmcplots)
library(tidyverse)
#install.packages("plyr")
#library(plyr)
library(dplyr)
library(Cairo)

## load workspace which has all of the data files already read in
load("AGHY_data_clean.RData") ## this has all of the raw data




## bundle data
AGHY_subplot_dat <- AGHY.merge%>%
  select(newplot,year_t,subplot,water,transmission,target_init_freq,E_plus_liberal,total)%>% 
  mutate(water = as.factor(case_when(water == "Add" ~ "Irrigated",
                                     water == "Control" ~ "Ambient"))) %>% 
  na.omit()%>%
  arrange(year_t,newplot,subplot)

AGHY_plot_dat <- AGHY.merge %>%
  select(newplot, water) %>% 
  distinct()%>% 
  mutate(water = case_when(water == "Add" ~ "Irrigated",
                           water == "Control" ~ "Ambient")) %>% 
  mutate(water = as.integer(as.factor(water))) %>% 
  arrange(newplot)

## levels
N.trt<-length(levels(AGHY_subplot_dat$water))
N.years<-length(unique(AGHY_subplot_dat$year_t))
N.plots<-length(unique(AGHY_subplot_dat$newplot))
N.obs<-nrow(AGHY_subplot_dat)

## data - predictor variables at the plot level
water<-AGHY_plot_dat$water
trans_reduce<-AGHY_plot_dat$trans_reduce
initial_prev<-AGHY_plot_dat$initial_prev

## data - response variable at subplot level
y.pos<-AGHY_subplot_dat$E_plus_liberal
N.samples<-AGHY_subplot_dat$total
plot<-AGHY_subplot_dat$newplot
year<-AGHY_subplot_dat$year_t-2013

## data for prediction
x.levels<-seq(0,1,0.01)
N.x.levels<-length(x.levels)

######################## SURVIVAL DATA PREP ###############################
mean_prev <-read.csv("mean_prev.csv")

## change mean_prev from wide to long 
mean_prev$newplot <-mean_prev$X
mean_prev_yr<- mean_prev %>% 
  group_by(newplot) %>% 
  gather(key = yearID, value = mean_prev, -newplot, -plot, -water,-X) %>% 
  mutate(year_t = case_when(yearID == "initial_prev" ~ "2013", ## add in years with which we can match and facet_grid 
                            yearID == "mean_2014" ~ "2014", ## the mean prev observed in 2014 
                            yearID == "mean_2015" ~ "2015",
                            yearID == "mean_2016" ~ "2016")) %>% 
  transform(year_t = as.numeric(year_t))




## survival data for individuals of known endo status

survival_dat <- AGHY.demo.final %>%
  filter(spring_survival_t1 != "NA") %>%
  select(newplot, subplot, ID, year_t, year_t1, spring_survival_t1, endo.stat) %>%
  left_join(AGHY.plots %>% 
              select(newplot, water)) %>%
  left_join(mean_prev_yr)

## need to incorporate endo stat
total_plants <- survival_dat %>% 
  group_by(year_t, newplot, mean_prev) %>% 
  select(spring_survival_t1, year_t, newplot, mean_prev, endo.stat) %>% 
  mutate(Ep_alive = case_when(endo.stat == 1 & spring_survival_t1 == 1 ~ 1),
         Em_alive = case_when(endo.stat == 0 & spring_survival_t1 == 1 ~ 0)) %>% 
  
  summarize(total_plants = length(spring_survival_t1), alive_plants = sum(spring_survival_t1), 
            Ep_known = sum(endo.stat, na.rm = T),
            Ep_alive = sum(Ep_alive, na.rm = T),
            Em_alive = sum(Em_alive, na.rm = T)) %>% 
  mutate(expected_Ep = total_plants * mean_prev,
         error_Ep = as.numeric(Ep_known > expected_Ep)) ## pretty good congruence between obsevered and expected Ep - only 1 that's observed as more Ep than expected and unrealistic


endo.unk <- survival_dat %>% 
  group_by(year_t, newplot, mean_prev) %>% 
  filter(is.na(endo.stat)) %>% 
  summarize(total_unknown = length(endo.stat),
            unk_alive = sum(spring_survival_t1))

surv_endo_unk_summarized <- endo.unk %>% 
  group_by(year_t) %>% 
  summarize(total_unknown = sum(total_unknown))

endo.known <- survival_dat %>% 
  group_by(year_t, newplot, mean_prev) %>% 
  filter(!is.na(endo.stat)) %>% 
  summarize(total_known = length(endo.stat),
            known_alive = sum(spring_survival_t1))

surv_endo_known_summarized <- endo.known %>% 
  group_by(year_t) %>% 
  summarize(total_known = sum(total_known))

surv_endo_summarized <- surv_endo_unk_summarized %>% 
  left_join(surv_endo_known_summarized) %>% 
  mutate(total_plants = total_known+ total_unknown,
         prop_unknown = total_unknown/total_plants,
         prop_known = total_known/total_plants) %>% 
  mutate(type = "Survival",
         year = case_when(year_t == 2014 ~ 1,
                         year_t == 2015 ~ 2))

survival_df_summarized <- total_plants %>% 
  full_join(endo.unk) %>% 
  full_join(endo.known) %>%
  mutate_at(vars(total_known, total_unknown), funs(replace(., is.na(.), 0))) %>% ## change the NAs to 0s in total_known and unknown -- there were a few instances in which all were known or unknown
  mutate(total_check = as.numeric(total_plants != (total_known + total_unknown)),
         expected_Em = total_plants - expected_Ep,
         Em_known = total_known - Ep_known)

surv_again <- survival_df_summarized %>% 
  mutate(expected_Ep = round(expected_Ep),
         expected_Em = round(expected_Em),
         diff_Ep = expected_Ep - Ep_known, ## how many Ep plants are "missing" -- expected number - number of known Ep
         diff_Em = expected_Em - Em_known) ## how many Em plants are "missing" -- expected number - number of known Em

surv_pop_prev_adj <- surv_again %>% 
  mutate(mean_prev_adj = diff_Ep/(diff_Ep + diff_Em)) %>% 
  select(newplot, year_t, mean_prev, mean_prev_adj ) %>% 
  ungroup()



### survival plot info for survival est for recruitment equation
### September 16, 2019

survival_dat2 <- survival_dat %>% 
  select(newplot, year_t, year_t1, water) %>%
  distinct() %>% 
  mutate(water = as.numeric(as.factor(water)),
         year = case_when(year_t1 == 2015 ~ 1,
                          year_t1 == 2016 ~ 2)) %>% 
  arrange(year_t1, newplot)

###### ^^ for the recruitment equation

## NUMBERS for UNKNOWN ENDO AND SURVIVAL
## in 2014 there are 189 plants in water and 191 in control &
## in 2015 there are 141 plants in water and 137 in control 

survival_dat %>% 
  group_by(year_t,endo.stat,water) %>% 
  select(spring_survival_t1) %>% 
  summarise(n()) 




## select the plants have have BOTH survival and endo data
survival_known <- survival_dat %>%
  filter(!is.na(endo.stat))%>%
  mutate(year = year_t -(min(year_t)-1),
         plot = newplot - (min(newplot)-1),
         endo = as.integer(endo.stat)+1,
         water = as.integer(water)) %>%
  select(year, plot, endo, water, spring_survival_t1, mean_prev, subplot, ID)%>%
  filter(!is.na(spring_survival_t1))

## select the plants for which there is survival info but we don't know their endo stat
## select the unknowns
survival_unk <- survival_dat %>%
  filter(is.na(endo.stat))%>%
  left_join(surv_pop_prev_adj) %>% 
  mutate(year = year_t - (min(year_t)-1),
         plot = newplot - (min(newplot)-1),
         water = as.integer(water)) %>%
  select(year, plot, endo.stat, water, spring_survival_t1, mean_prev, mean_prev_adj, subplot, ID)%>%
  filter(!is.na(spring_survival_t1))%>%
  unite(unique_ID, plot, subplot,year, sep = "_", remove = F) 


## combine to the subplot level so that we have number of individuals that survived and total alive last year
## these are just individuals for which we don't know their endo stat  
survival_unk_sbplt <- survival_unk %>% 
  group_by(unique_ID) %>% 
  summarize(alive.t1 = sum(spring_survival_t1),
            alive.t = length(unique_ID),
            plot.prev = unique(mean_prev),
            plot.prev.adj = unique(mean_prev_adj),
            water = unique(water)) %>% 
  mutate(plot.prev.adj_bound = ifelse(plot.prev.adj > 1, 1,
                                      ifelse(plot.prev.adj < 0, 0, plot.prev.adj))) %>% ## bound the adj plot prev by 0 and 1
  separate(unique_ID, c("plot", "subplot", "year"), sep = "_", remove = F) %>% 
  mutate(year = (as.numeric(year))) 


## adjust the population prev

survival_known_plot <- survival_dat %>%
  filter(!is.na(endo.stat))%>%
  mutate(year = year_t -(min(year_t)-1),
         plot = newplot - (min(newplot)-1),
         endo = as.integer(endo.stat)+1,
         water = as.integer(water)) %>%
  select(year, plot, endo, water, spring_survival_t1, mean_prev, subplot, ID, endo.stat)%>%
  filter(!is.na(spring_survival_t1)) %>% 
  group_by(year,plot,water) %>% 
  summarize(total_known = length(spring_survival_t1),
            total_known_ep = sum(endo.stat),
            total_known_em = total_known - total_known_ep)

## select the plants for which there is survival info but we don't know their endo stat
## select the unknowns
survival_unk_plot <- survival_dat %>%
  filter(is.na(endo.stat))%>%
  left_join(surv_pop_prev_adj) %>% 
  mutate(year = year_t - (min(year_t)-1),
         plot = newplot - (min(newplot)-1),
         water = as.integer(water)) %>%
  select(year, plot, endo.stat, water, spring_survival_t1, mean_prev, mean_prev_adj, subplot, ID)%>%
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



#####################################################################
################### FERTILITY AND RECRUITMENT #######################
#####################################################################

###### SEED MASS FERTILITY DATA

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
  select(-water) %>% 
  anti_join(seedmass_dat_issues) %>% ## drop these rows for which the number of infs collected is greater than observed
  left_join(mean_prev_yr) ## add the estiamted mean pop prevs -- ultimately this will be done in the bayes model

seedmass_dat %>% 
  group_by(year_t,endo.stat,water) %>% 
  select(seed_mass_t) %>% 
  summarise(n())

## There are 44 individuals for which we have seed mass but no endo score
seedmass_unk <- seedmass_dat %>% 
  filter(is.na(endo.stat)) %>% 
  mutate(year = year_t - (min(year_t)-1),
         plot = newplot - (min(newplot)-1),
         water = as.integer(water)) %>% 
  select(year,plot,endo.stat,water,seed_mass_t,inf_collected_t, mean_prev) %>% 
  filter(!is.na(seed_mass_t)) ## there are 14 individuals for which we have neither seed mass nor endo stat (can't do much with these)

seedmass_dat <- seedmass_dat%>% 
  mutate(year = year_t - (min(year_t)-1),
         plot = newplot - (min(newplot)-1),
         endo = as.integer(endo.stat+1),
         water = as.integer(water)) %>% 
  select(year,plot,endo,water,seed_mass_t,inf_collected_t, mean_prev) %>% 
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
  select(-water) %>% 
  left_join(mean_prev_yr) %>% 
  mutate(flowering_t = if_else(inf_number_t >= 1, 1, 0)) ## indicate whether the plant flowered (1) or not (0)

ggplot(infs_dat)+
  geom_histogram((aes(x=inf_number_t)))+
  facet_grid(endo.stat~water*year_t)

infs_dat %>% 
  group_by(year_t,endo.stat) %>% 
  select(inf_number_t) %>% 
  summarise(n())



## data for adjusting population prevalence for inf production 

total_plants_infs <- infs_dat_pop_prev %>% 
  group_by(year_t, newplot, mean_prev) %>% 
  select(inf_number_t, inf_collected_t, year_t, newplot, mean_prev, endo.stat, subplot, ID, flowering_t) %>% 
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
 # mutate(expected_Ep_flow = total_plants_flow * mean_prev,
  #       error_Ep = as.numeric(Ep_known_flow > expected_Ep_flow)) ## pretty good congruence between obsevered and expected Ep - only 1 that's observed as more Ep than expected and unrealistic


## numbers for MS
total_infs_collected <- infs_dat_pop_prev %>% 
  filter(inf_collected_t > 1) %>% 
  summarize(total_plants = n(),
            total_infs = sum(inf_collected_t))


total_infs_counted <- infs_dat_pop_prev %>% 
  filter(inf_number_t > 1) %>% 
  summarize(total_plants = n(),
            total_infs = sum(inf_number_t))

## pull out the data for which we don't have endo scores
inf.unk <- infs_dat_pop_prev %>% 
  group_by(year_t, newplot, mean_prev) %>% 
  filter(is.na(endo.stat)) %>% 
  summarize(total_unknown_infs = length(endo.stat),
            unk_flow = sum(flowering_t))

inf.known <- infs_dat_pop_prev %>% 
  group_by(year_t, newplot, mean_prev) %>% 
  filter(!is.na(endo.stat)) %>% 
  summarize(total_known_infs = length(endo.stat),
            known_flow = sum(flowering_t))

flow_df_summarized <- total_plants_infs %>% 
  full_join(inf.unk) %>% 
  full_join(inf.known) %>%
  mutate_at(vars(total_known_infs, total_unknown_infs), funs(replace(., is.na(.), 0))) %>% ## change the NAs to 0s in total_known and unknown -- there were a few instances in which all were known or unknown
  mutate(total_check = as.numeric(total_plants_flow != (total_known_infs + total_unknown_infs)))







## there are 681 flowering individuals for which we do not know their endo stat
infs_unk <- infs_dat_pop_prev %>% 
  filter(is.na(endo.stat),
         inf_number_t >= 1) %>% ## select only plants that flowered 
  mutate(year = year_t - (min(year_t)-1),
         plot = newplot - (min(newplot)-1),
         endo = 9999,
         water = as.integer(water)) %>% 
  select(year,plot,endo,water,inf_number_t, mean_prev) 

infs_dat <- infs_dat_pop_prev %>% 
  filter(inf_number_t >= 1) %>% 
  mutate(year = year_t - (min(year_t)-1),
         plot = newplot - (min(newplot)-1),
         endo = as.integer(endo.stat+1),
         water = as.integer(water)) %>% 
  select(year,plot,endo,water,inf_number_t, mean_prev) %>% 
  na.omit()


infs_endo_unk_summarized <- infs_unk %>% 
  group_by(year) %>% 
  summarize(total_unknown = sum(inf_number_t))


infs_endo_known_summarized <- inf.known %>% 
  group_by(year_t) %>% 
  summarize(total_known = sum(total_known_infs)) %>% 
  mutate(year = case_when(year_t == "2014" ~ 1,
                          year_t == "2015" ~ 2))

infs_endo_summarized <- infs_endo_unk_summarized %>% 
  left_join(infs_endo_known_summarized) %>% 
  mutate(total_plants = total_known+ total_unknown,
         prop_unknown = total_unknown/total_plants,
         prop_known = total_known/total_plants) %>% 
  mutate(type = "Infs")


## because we want the contribution of eps.infs.overdisp in the weighted mean
## we need to make the vector length start at the end of the N.infs.unk.vector
## an easy way to do this is to stack the dfs on top of one another
## then tell the for loop to start runing at the term that starts the unk df

infs_unk_stacked <- infs_dat %>% 
  bind_rows(infs_unk) 

## previously added in the population prev - but now this is all done in a single model
#%>% 
#  left_join(flow_pop_prev_adj) %>%  ## add in the adjusted population prevs
#  mutate(mean_prev_adj_bound = ifelse(mean_prev_adj > 1, 1,
#                                      ifelse(mean_prev_adj < 0, 0, mean_prev))) %>%   ## bound the adj plot prev by 0 and 1
#  mutate_at(vars(mean_prev_adj_bound), funs(replace(., is.na(.), 1))) ## remove the NAs -- replacing them with 1, this shouldn't matter since we start the model at the point where the endo stat is unknown (at 461)
##########
#########




flowering<- AGHY.merge %>% 
  select(year_t, water, newplot, subplot, spring_flowering_recruit_count_per_subplot,spring_vegetative_recruit_count_per_subplot) %>% 
  filter(!is.na(spring_vegetative_recruit_count_per_subplot), !is.na(spring_flowering_recruit_count_per_subplot)) %>%  ## drop subplots for which there is an NA for either flowering or veg. count
  mutate(total.plants = spring_flowering_recruit_count_per_subplot + spring_vegetative_recruit_count_per_subplot ) %>% 
  filter(year_t != 2016) ## DROP 2016 flowering since it doesn't contribute to a paired vital rate and will be easier to index everything in the combined model 


## select the plants have have BOTH survival and endo data
flowering_dat <- flowering %>%
  mutate(year = year_t -(min(year_t)-1),
         plot = newplot - (min(newplot)-1),
         water = as.integer(water))

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
         water = as.integer(water),
         N.samples = 1)

flowering_sbplts <- flowering_individuals %>%  ## select the individuals for which we KNOW endo stat and have info on whether or not they were flowering, AND that were ALIVE in t1
  group_by(newplot, subplot, year_t, flow_t)%>%
  summarize(n(), count = length(newplot),
            num_ep = sum(endo.stat),
            total_known = length(endo.stat),
            endo_flowering =  sum(flow_t[endo != "NA"])) %>% 
  ## endo_flowering is the number of plants of known endo stat that are flowering, typically this is equal to the number of known endo individuals since we get endo from seeds (there's one case in which )
  full_join(flowering_dat)

## check how many E+ and E- plants are known from each subplot
flowering_check <- flowering_sbplts %>% 
  mutate(num_em = total_known - num_ep)

total_known <- flowering_check %>%
  group_by(year, water) %>% 
  filter(!is.na(num_ep), !is.na(num_em)) %>% 
  summarize(num_ep = sum(num_ep), 
            num_em = sum(num_em))

totals_flowering <- flowering_check %>% 
  mutate(total_unknown = spring_flowering_recruit_count_per_subplot - total_known) %>% 
  group_by(year, water) %>% 
  filter(!is.na(total_unknown),
         !is.na(total_known)) %>% 
  summarize(total_unk = sum(total_unknown),
            total_known = sum(total_known),
            total_plants = sum(spring_flowering_recruit_count_per_subplot),
            total_ep = sum(num_ep),
            total_em = sum(num_em))

ggplot(flowering_check)+
  geom_histogram((aes(x=num_ep)))+
  facet_grid(water~year_t)

ggplot(flowering_check)+
  geom_histogram((aes(x=num_em)))+
  facet_grid(water~year_t)

## check for disagreement between number of demography individuals and subplot counts
issues_flow <- flowering_sbplts %>% 
  mutate(issue = ifelse(count > total.plants, 1, 0)) ## There's one issue where the subplot counts are less than the number of demo plants (plot 1)

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
## Jan 9, 2019
## get all the plots together 

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

## pull out the unknowns to for weighted mean approach

flow_unk_plot <- flowering_plts %>%
  filter(total_known_f == 0)


flow_endo_summarized <- flowering_plts %>% 
  group_by(year) %>% 
  summarize(total_plants = sum(total_plants_f),
            total_unknown = sum(total_unknown_f),
            total_known = sum(total_known_f),
            prop_known = total_known/total_plants,
            prop_unknown = total_unknown/total_plants) %>% 
  mutate(type = "Flowering")


summary_endo_status_MS <- surv_endo_summarized %>% 
  bind_rows(infs_endo_summarized) %>% 
  bind_rows(flow_endo_summarized)


## DATA FRAME FOR LINEAR REGRESSION CONVERSION FROM SEED MASS TO SEED COUNTS

linreg_dat <- AGHY.linreg %>% 
  dplyr::rename(plot = Plot) %>% 
  left_join(AGHY.plots) %>% 
  mutate(water = as.integer(water),
         endo = Endo + 1)


### RECRUITMENT IS OFF BY A YEAR -- it should be in 2015 and 2016 not 2014 and 2015
## data for recruitment - 
#rec_dat <- flowering_dat %>% 
#  group_by(water, plot, year) %>% 
#  summarize(flowering_plot = sum(spring_flowering_recruit_count_per_subplot),
#            veg_plot = sum(spring_vegetative_recruit_count_per_subplot),
#            total_plants_3sbplts = sum(total.plants),
#            n_sbplts = n()) %>% 
#  mutate(plant_density_sbplt = total_plants_3sbplts / n_sbplts) %>% 
#  unite(year_water, c("year", "water"), remove = FALSE) %>% 
#  arrange(year_water)

flow_15 <- flow.15 %>% 
  select(year_t, water.x, newplot, subplot, spring_flowering_recruit_count_per_subplot, spring_vegetative_recruit_count_per_subplot, total.plants) %>% 
  dplyr::rename(water = water.x)

flow_16 <- flow.16 %>% 
  select(year_t, water.x, newplot, subplot, spring_flowering_recruit_count_per_subplot, spring_vegetative_recruit_count_per_subplot, total.plants) %>% 
  dplyr::rename(water = water.x)

flow_15_16 <- flow_15 %>% 
  bind_rows(flow_16) %>% 
  filter(year_t != 2014) %>% 
  mutate(year = year_t -(min(year_t)-1),
         plot = newplot - (min(newplot)-1),
         water = as.integer(water))

rec_dat <- flow_15_16 %>% 
  group_by(water, plot, year) %>% 
  summarize(flowering_plot = sum(spring_flowering_recruit_count_per_subplot),
            veg_plot = sum(spring_vegetative_recruit_count_per_subplot),
            total_plants_3sbplts = sum(total.plants),
            n_sbplts = n()) %>% 
  mutate(plant_density_sbplt = total_plants_3sbplts / n_sbplts) %>% 
  unite(year_water, c("year", "water"), remove = FALSE) %>% 
  arrange(year_water)



### for new version of recruiment 
library(readxl)
recruitment_raw <- read_excel("AGHY_SFAEF_life_history_expt.xlsx", sheet = "Subplot-level data")

plot_info_r <- AGHY_plot_dat %>% 
  mutate(plot = newplot) %>% 
  arrange(plot) %>% 
  distinct()

rec_final <- recruitment_raw %>% 
  mutate(plot = as.numeric(as.factor(plot)),
         year = as.numeric(as.factor(year))) %>% 
  mutate(total_plants_in_sbplt = spring_flowering_recruit_count_per_subplot + spring_vegetative_recruit_count_per_subplot) %>% 
  filter(total_plants_in_sbplt != "NA") %>%  ## drop the 6 rows for which we don't have an observation of counts
  ungroup() %>% 
  dplyr::left_join(plot_info_r, by = "plot") %>% 
  arrange(year, plot) 


rec_plot_final <- rec_final %>% 
  group_by(plot,water,year) %>% 
  summarize(total_plants_sbplt = sum(total_plants_in_sbplt)) %>% 
  arrange(year, plot)

### df that's the correct length for the rec estimate year t+1
rec_plot_final_yr1 <- rec_plot_final %>% 
  filter(year == 2 |
           year == 3)

rec_plot_final_yr <- rec_plot_final %>% 
  filter(year == 1 | 
           year == 2)

### for neg binom model of recruitment 
x <- rnorm(n=N.sbplts.avg.rec)
dat <- data.frame(x)
X=cbind(1,dat$x)


###########################################################################
### data for vertical transmission ###
######################################

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


vtrans_check <- vtrans_dat %>% 
  dplyr::group_by(year, water) %>% 
  summarize(total_seeds_plt = sum(total_seeds_scored),
            total_ep_plt = sum(total_ep_seeds)) %>% 
  group_by(year,water) %>% 
  mutate(vtrans = total_ep_plt /total_seeds_plt)



###########################################################################
###########################################################################


jag.data<-list(N.trt=N.trt,
               N.years=N.years,
               N.plots=N.plots,
               N.obs=N.obs,
               water=water,
               y.pos=y.pos,
               N.samples=N.samples,
               plot=plot,
               year=year,
               x.levels=x.levels,
               N.x.levels=N.x.levels,
               
               
               ### SURVIVAL DATA
               
               N.endo = max(survival_known$endo),
               N.water = max(survival_known$water),
               N.year = max(survival_known$year), 
               N.plot = max(survival_known$plot),
               
               N.obs.surv.known = nrow(survival_known),
               survival.k = survival_known$spring_survival_t1,
               # survival.known = survival_known$spring_survival_t1,
               endo.k = survival_known$endo,
               water.k = survival_known$water,
               year.k = survival_known$year,
               plot.k = survival_known$plot,
               
               N.obs.surv.unknown = nrow(survival_unk_sbplt),
               alive.unk.t1 = survival_unk_sbplt$alive.t1,
               alive.unk.t = survival_unk_sbplt$alive.t,
               #   alive.unk.t1.2 = survival_unk_sbplt$alive.t1,
               #   alive.unk.t.2 = survival_unk_sbplt$alive.t,
               #  plot.prev.surv = survival_unk_sbplt$plot.prev.adj, ## using the adj. pop prev (** decide if it should be bound between 0 and 1)
               water.unk.surv = survival_unk_sbplt$water,
               year.unk.surv = as.integer(survival_unk_sbplt$year),
               plot.unk.surv = as.integer(survival_unk_sbplt$plot),
               
               
               ## ADJUST POP PREV FOR SURVIVAL
               N.plots.known.endo.surv = nrow(survival_plot),
               year.adj.surv = survival_plot$year,
               total.demo.plants.3.sbplts = survival_plot$total_demo_plants,
               plot.surv.adj = survival_plot$plot,
               known_ep = survival_plot$total_known_ep,
               total_unknown_surv = survival_plot$total_unknown,
               
               #### FERTILITY, FLOWERING AND RECRUITMENT
               
               ### FERTILITY DATA
               #N.endo = max(seedmass_dat$endo),
               #N.water = max(seedmass_dat$water),
               #N.year.fert = max(seedmass_dat$year),
               N.obs.seed = nrow(seedmass_dat),
               # N.plot = max(seedmass_dat$plot),
               
               #seedmass = seedmass_dat$seed_mass_t,
               #infs_massed = seedmass_dat$inf_collected_t,
               endo.seed = seedmass_dat$endo,
               water.seed = seedmass_dat$water,
               year.seed = seedmass_dat$year,
               plot.seed = seedmass_dat$plot,
               y.seed = log(seedmass_dat$seed_mass_t/seedmass_dat$inf_collected_t),
               
               N.obs.infs = nrow(infs_dat),
               endo.infs = infs_dat$endo,
               water.infs = infs_dat$water,
               year.infs = infs_dat$year,
               plot.infs = infs_dat$plot,
               y.infs = infs_dat$inf_number_t,
               
               # ADJUST POP Prev by known endo stat for inf production at the plot level
               
               N.plots.known.endo.inf = nrow(total_plants_infs),
               year_unk_adj = total_plants_infs$year_t,
               total_plants_infs = total_plants_infs$total_plants_flow,
               plot.inf.adj = total_plants_infs$newplot,
               known_ep_inf = total_plants_infs$Ep_known_flow,
               total_unknown_inf = total_plants_infs$total_flowering_plants_unknown,
               
               ## data for weighted mean approach for unknown infs
               
               infs.unk = infs_unk_stacked$inf_number_t,
               water_unk = infs_unk_stacked$water,
               plot_unk = infs_unk_stacked$plot,
               year_unk = infs_unk_stacked$year,
               N.obs.infs.unk = nrow(infs_unk_stacked),
               
               ### FLOWERING DATA
               # N.water = max(flowering_dat$water),
               # N.endo = 2,
               # N.year = max(flowering_dat$year),
               # N.plot = max(flowering_dat$plot),
               
               ### FLOWERING DATA FOR KNOWN ENDO 
               N.obs.flow.known = nrow(flowering_individuals),
               endo.flow.k = flowering_individuals$endo,
               water.flow.k = flowering_individuals$water,
               year.flow.k = flowering_individuals$year,
              # plot.flow.k = flowering_individuals$plot,
               flowering.k = flowering_individuals$flow_t,
             # flowering.k.new = rbinom(82, 1, 0.5), ## fake data to try and see what's wrong with the model
               N.samples.flow = flowering_individuals$N.samples, ## for betabin - N is 1
               
               ### FLOWERING DATA FOR ADJUSTING POP PREV
               N.plots.known.endo.flow = nrow(flowering_plts),
               plot.adj.flow = flowering_plts$plot,
               year.adj.flow = flowering_plts$year,
               total.plants.3.sbplts.f = flowering_plts$total_plants_f,
               known_ep_f = flowering_plts$total_known_ep_f,
               total_unknown_flow = flowering_plts$total_unknown_f,
               
               ### FLOWERING DATA FOR WEIGHTED MEAN APPROACH
               
               N.obs.flow.unknown = nrow(flowering_sbplts),
               year.f = flowering_sbplts$year_t,
               plot.f = flowering_sbplts$newplot,
               water.f = flowering_sbplts$water,
               total.unk.plants.f = flowering_sbplts$total_unknown,##########################HERE
               flow.unk = flowering_sbplts$total_unknown_flowered,
             
               
               ## linear regression seed mass and seed count data
               
               N.obs.seed.lin.reg = nrow(linreg_dat),
               #endo.lin.reg = linreg_dat$endo, ## currently not splitting by endo and water
               #water.lin.reg = linreg_dat$water,
               seed.count.linreg = (linreg_dat$seed_count),
               seed.mass.linreg = (linreg_dat$seed_mass),
               
               
               # ## recruitment data for determing endo stat 
               # N.sbplts.avg.rec = nrow(flowering_dat),
               # rec.sbplt = flowering_dat$total.plants + 1,
               # water.r = flowering_dat$water,
               # plot.r = flowering_dat$plot,
               # year.r = flowering_dat$year,
               # observed.plants = rec_dat$total_plants_3sbplts,
               # # plot.prev.rec = rec_dat$mean_prev,
               # plots.rec = rec_dat$plot,
               # N.plots.rec = nrow(rec_dat),
               # water.rec = rec_dat$water,
               # year.rec = rec_dat$year,
             
             ## NEW rec
             N.sbplts.avg.rec = nrow(rec_final),
             rec.sbplt = rec_final$total_plants_in_sbplt,
             water.r = rec_final$water,
             plot.r = rec_final$plot,
             year.r = rec_final$year,
             
             plots.rec_t = rec_plot_final_yr$plot,
             N.plots.rec_t = nrow(rec_plot_final_yr),
             water.rec_t = rec_plot_final_yr$water,
             year.rec_t = rec_plot_final_yr$year,
             
             surv.rec.plot = survival_dat2$newplot,
             N.plots.surv.rec = nrow(survival_dat2),
             water.surv.rec = survival_dat2$water,
             year.surv.rec = survival_dat2$year,
             
             
             plots.rec_t1 = rec_plot_final_yr1$plot,
             N.plots.rec_t1 = nrow(rec_plot_final_yr1),
             water.rec_t1 = rec_plot_final_yr1$water,
             year.rec_t1 = rec_plot_final_yr1$year,
             ### neg binom rec
              X = X,
              mu.beta=rep(0,2),
              tau.beta=diag(.0001,2),
         
               
               ##### VERTICAL TRANSMISSION
               water.vt = vtrans_dat$water,
               year.vt = vtrans_dat$year,
               plot.vt = vtrans_dat$newplot,
               total.seeds = vtrans_dat$total_seeds_scored,
               ep.seeds = vtrans_dat$total_ep_seeds,
               N.obs.vt = length(vtrans_dat$total_ep_seeds))



## Inits function
inits<-function(){list(beta0.mean=matrix(rnorm(N.trt*(N.years-1),0,2),N.trt,(N.years-1)),
                       beta1.mean=matrix(rnorm(N.trt*(N.years-1),0,2),N.trt,(N.years-1)),
                       sigma0.plot=runif(1,0,10),
                       sigma1.plot=runif(1,0,10),
                       a=runif(1,0,10),
                       
                       ## survival inits
                       mu_surv = array(rnorm(jag.data$N.endo*jag.data$N.water*jag.data$N.year),
                                       dim=c(jag.data$N.endo,jag.data$N.water,jag.data$N.year)),
                       
                       #sigma.surv = rlnorm(jag.data$N.year),
                       sigma.surv = rlnorm(1),
                       ## FLOWERING AND RECRUITMENT
                       
                       mu.seed = array(rnorm(jag.data$N.endo*jag.data$N.water*jag.data$N.year),
                                       dim=c(jag.data$N.endo,jag.data$N.water,jag.data$N.year)),
                       mu.infs = array(rnorm(jag.data$N.endo*jag.data$N.water*jag.data$N.year),
                                       dim=c(jag.data$N.endo,jag.data$N.water,jag.data$N.year)),
                       sigma.res = rlnorm(1),
                       sigma.infs.overdisp = rlnorm(1),
                    #   sigma.plot = rlnorm(jag.data$N.year),
                       sigma.plot = rlnorm(1),
                      # sigma.plot.infs = rlnorm(jag.data$N.year),
                       sigma.plot.infs = rlnorm(1),
                       mu_flow = array(rnorm(jag.data$N.endo*jag.data$N.water*jag.data$N.year),
                                       dim=c(jag.data$N.endo,jag.data$N.water,jag.data$N.year)),
                       #sigma.flow = rlnorm(jag.data$N.year),
                       sigma.flow = rlnorm(1),
                       b=runif(1,0,10),
                      # c = runif(1,0,10),
                       
                       beta0_seed_linreg = rnorm(1,0,1),
                       slope_lin_reg = rnorm(1,0,1),
                       
                       sigma.count = rlnorm(1),
                       
                       ## recruitment
                       mu_rec = array(rnorm(jag.data$N.water*jag.data$N.years),
                                      dim=c(jag.data$N.water,jag.data$N.years)),
                    
                      # mu.beta=rep(0,2), ## thinking of using for flowering neg binom
                      # tau.beta=diag(.0001,2), ## thinking of using for flowering neg binom  
                       #sigma.rec = rlnorm(jag.data$N.year),
                       sigma.rec = rlnorm(1),
                       ## VERTICAL TRANSMISSION
                       mu_vtrans = array(rnorm(jag.data$N.water*jag.data$N.year),
                                         dim=c(jag.data$N.water,jag.data$N.year)),
                      # sigma.vt = rlnorm(jag.data$N.year)
                      sigma.vt = rlnorm(1)
)
}

## Params to estimate
parameters<-c("beta0.mean","beta1.mean", ## Population prev params
               "sigma0.plot","sigma1.plot","a",
               "Eplus.add.pred","Eplus.control.pred",
               "prev","fit","fit.new",
               
               ## SURVIVAL PARAMS
               
               "mu_surv", "sigma.surv", "weighted.mean.surv",
               "survival.vital.rate",
               "fit.surv", "fit.surv.new",
              
               ## FERTILITY AND RECRUITMENT
               "mu.seed",
               "mu.infs",
               "prob.ep.f",
               "prob.em.f",
             #  "ep_seeds",
            #   "em_seeds",
              # "ep_flowering",
              # "em_flowering",

#               ## fertility vital rate
#               # "em_seeds_ambient.14",
#               # "em_seeds_ambient.15",
#               # "em_seeds_irrigated.14",
#               # "em_seeds_irrigated.15",
#               #
#               # "ep_seeds_ambient.14",
#               # "ep_seeds_ambient.15",
#               # "ep_seeds_irrigated.14",
#               # "ep_seeds_irrigated.15",
#
               "fit.seed",
               "fit.seed.new",

               "fit.infs",
               "fit.infs.new",

               "fit.infs.unk",
               "fit.infs.unk.new",

               "fit.flow",
               "fit.flow.new",

               "fit.linreg",
               "fit.linreg.new",

               #"prob_ep_rec",
               #"prob_em_rec",
               "prob.ep.surv", #### new version
               "prob.em.surv", ### new version
               "ep_rec_plt",
              "em_rec_plt",
               "flowering.vital.rate",
               "seedmass.per.cap",
               #"avg.rec",
               "fit.rec",
               "fit.rec.new",
#
#               ## VERTICAL TRANSMISSION
               "vtrans",
#
               "fit.vtrans",
               "fit.vtrans.new",
#
#               ## adjusted pop prev
#               "prev_adj",
#               "prev_adj_f",

              ## Recruitment probability (derived)

              # "new_prob_ep_rec.ambient.14",
              # "new_prob_ep_rec.irrigated.14",
              #
              # "new_prob_em_rec.ambient.14",
              # "new_prob_em_rec.irrigated.14",
              #
              # "new_prob_ep_rec.ambient.15",
              # "new_prob_ep_rec.irrigated.15",
              #
              # "new_prob_em_rec.ambient.15",
              # "new_prob_em_rec.irrigated.15",
              "p.rec",
              "ep_rec_plt_t1",
              "ep_rec_plt_t",
              "ep_rec_flow_t",
              "ep_seeds_plt_t",
              
              "ep_seeds_rec",
              "ep_rec_plt",

             "mu_vtrans",
              "mean.vtrans",

              

              # "ep_rec_plt",
              # "em_rec_plt",
               "new_prob_ep_rec",
              # "ep_seeds_rec",
               "prob_ep_rec_new.ambient.14",
               "prob_ep_rec_new.ambient.15",
              #
             "prob_ep_rec_new.irrigated.14",
              "prob_ep_rec_new.irrigated.15",

              "prob_em_rec_new.ambient.14",
              "prob_em_rec_new.ambient.15",
              #
              "prob_em_rec_new.irrigated.14",
              "prob_em_rec_new.irrigated.15",
              

              #### random effects to check fequency dependence
               "eps.surv",
               "eps.flow",
               "eps.rec",
               "eps.seed",
               "eps.infs",
               "eps.vt"

)


## MCMC settings
ni<-10000
nb<-2000
nt<-20
nc<-3

## run JAGS
AGHY.endochange.20142016<-jags(data=jag.data,inits=inits,parameters.to.save=parameters,
                               model.file="./vital_rate_scripts/Bayes models/AGHY_bayes_full_updated_MD_NM2.txt",
                               n.thin=nt,n.chains=nc,n.burnin=nb,
                               n.iter=ni,working.directory=getwd())
#mcmcplot(AGHY.endochange.20142016)

## save the workspace to load for figure creation
#save.image("AGHY_model_Bayes_full_output_update.RData")

#load("AGHY_model_Bayes_full_output_update.RData")

bayes.data.summary.surv<-as.data.frame(AGHY.endochange.20142016$BUGSoutput$summary, keep.rownames=TRUE)[]

bayes.data.summary.surv<-cbind(rownames(bayes.data.summary.surv),bayes.data.summary.surv)
names(bayes.data.summary.surv)[names(bayes.data.summary.surv) =="rownames(bayes.data.summary.surv)"] <- "info"



plot(AGHY.endochange.20142016$BUGSoutput$sims.list$fit, AGHY.endochange.20142016$BUGSoutput$sims.list$fit.new,
     xlab="Discrepancy for actual data",
     ylab="Discrepancy for new data",xlim=c(600,1200),ylim=c(600,1200))
abline(0,1, col='darkgray',lwd=3)

## save posterior predictive check
Cairo(file="PP_pop_prev.png", type="png", bg="white")
plot(AGHY.endochange.20142016$BUGSoutput$sims.list$fit, AGHY.endochange.20142016$BUGSoutput$sims.list$fit.new,
     xlab="Discrepancy for actual data",
     ylab="Discrepancy for new data",xlim=c(600,1200),ylim=c(600,1200))
abline(0,1, col='darkgray',lwd=3)
dev.off() # creates a file "plot.png" with the above plot



prev_hold<-data.frame(rbind(AGHY.endochange.20142016$BUGSoutput$sims.list$prev[,,1],
                            AGHY.endochange.20142016$BUGSoutput$sims.list$prev[,,2],
                            AGHY.endochange.20142016$BUGSoutput$sims.list$prev[,,3]))
colnames(prev_hold)<-unique(plot)

prev_dat <- prev_hold %>%
  mutate(year = rep(2014:2016,each=dim(AGHY.endochange.20142016$BUGSoutput$sims.list$prev)[1]))%>%
  gather('1':'47',key="newplot",value="prevalence")%>%
  mutate(plot = as.integer(newplot),
         newplot = as.integer(newplot))%>%
  dplyr::group_by(newplot,year)%>%
  mutate(mean.prev = mean(prevalence),
         low.CI = quantile(prevalence,probs=0.05),
         high.CI = quantile(prevalence,probs=0.95))%>%
  select(plot, newplot, water, mean.prev, low.CI, high.CI) %>% 
  distinct() %>% 
  left_join(AGHY_plot_dat)%>%
  mutate(water = ifelse(water==1, "Ambient", "Irrigated")) %>% 
  ungroup()

#%>%
#mutate(trans_reduce = ifelse(year==2014,trans_reduce,NA))%>%
#filter(trans_reduce == 0 | is.na(trans_reduce))

mean_prev<-prev_dat%>%
  mutate(year=paste('mean', year,sep="_"))%>%
  select(plot,water,mean.prev,year,initial_prev)%>%
  spread(year,mean.prev)
lowCI_prev<-prev_dat%>%
  mutate(year=paste('lowCI', year,sep="_"))%>%
  select(plot,water,low.CI,year)%>%
  spread(year,low.CI)
highCI_prev<-prev_dat%>%
  mutate(year=paste('highCI', year,sep="_"))%>%
  select(plot,water,high.CI,year)%>%
  spread(year,high.CI)  
prev_wide <- left_join(left_join(mean_prev,lowCI_prev,by="plot"),highCI_prev,by="plot") 


## APRIL 22, 2018 -- write the mean_prev df to csv, to that we can use it to build the vital rates script
## then the WORKING VR scripts will be merged with this one and run as a whole, so that all of the uncertainty will be estimated together
#write.csv(mean_prev, "mean_prev.csv")

##get the CI for the Eplus.control.pred and add predictions 

CIS <- as.data.frame(prev_hold) %>% 
  filter(str_detect(info, "^Eplus.")) %>% 
  rename(low_CI = '2.5%',
         high_CI = '97.5%')

predictions <- tibble(x.prev = rep(x.levels,(N.years-1)*N.trt),
                      pred_prev = c(as.vector(AGHY.endochange.20142016$BUGSoutput$mean$Eplus.control.pred),
                                    as.vector(AGHY.endochange.20142016$BUGSoutput$mean$Eplus.add.pred)),
                      year = rep(2015:2016,each=N.x.levels,times=N.trt),
                      water = c(rep("Irrigated",N.x.levels*(N.years-1)),rep("Ambient",N.x.levels*(N.years-1)))) %>% 
  mutate(year_previous = year - 1) %>% 
  unite(water_year, c(water, year_previous), remove = FALSE) 



ggplot(prev_wide)+
  geom_point(aes(x=mean_2014,y=mean_2015,size=1),color="tomato")+
  geom_point(aes(x=mean_2015,y=mean_2016,size=1),color="aquamarine3",inherit.aes = FALSE)+
  facet_wrap(~water)+
  geom_line(data=predictions,aes(x=x.prev,y=pred_prev,color=as.factor(year)),size=1.5, inherit.aes = FALSE)+
  geom_line(data=predictions,aes(x=x.prev,y=x.prev),inherit.aes = FALSE)+
  geom_segment(aes(x=mean_2014,y=mean_2015,xend=mean_2015,yend=mean_2016),arrow=arrow(length=unit(0.35,"cm")))+
  theme_bw()



mean_prev_yr<- prev_wide %>% 
  select(plot, mean_2014, mean_2015, mean_2016, water) %>% 
  group_by(plot, water) %>% 
  gather(key = yearID, value = mean_prev, -plot, -water) %>% 
  mutate(year_t = case_when(yearID == "initial_prev" ~ "2013", ## add in years with which we can match and facet_grid 
                            yearID == "mean_2014" ~ "2014", ## the mean prev observed in 2014 
                            yearID == "mean_2015" ~ "2015",
                            yearID == "mean_2016" ~ "2016")) %>% 
  transform(year_t = as.numeric(year_t))

prev_14_15 <- mean_prev_yr %>% 
  filter(year_t != 2016)

yrt1<- mean_prev_yr %>% 
  filter(year_t != 2014)%>%
  mutate(year_t1 = year_t, mean_prev_t1 = as.numeric(mean_prev)) %>% 
  select(year_t1, mean_prev_t1) %>% 
  bind_cols(prev_14_15) %>%
  unite(water_year, c(water, year_t), remove = FALSE) %>% 
  mutate(water_year = as.factor(water_year))


pop_prev_vectors <-ggplot(yrt1, aes(x = mean_prev, y = mean_prev_t1))+
  geom_point(aes(fill = (water_year)), shape=21, size = 3) + 
  scale_fill_manual(values=c("#ffae19","#dd7509", "#68BAE9", "#0338a8")) + 
  scale_colour_manual(values=c("#ffae19","#dd7509","#68BAE9", "#0338a8"))+
  geom_line(data=predictions,aes(x=x.prev,y=pred_prev,color=(water_year)),size=1.5, alpha=0.75, inherit.aes = FALSE)+
  geom_line(data=predictions,aes(x=x.prev,y=x.prev),inherit.aes = FALSE)+
  geom_segment(data =prev_wide, aes(x=mean_2014,y=mean_2015,xend=mean_2015,yend=mean_2016),arrow=arrow(length=unit(0.35,"cm")))+
  facet_wrap(~water)+
  xlab("Mean symbiont prevalence year_t")+
  ylab("Mean symbiont prevalence year_t+1")+
  theme_classic()+
  theme(legend.title = element_blank())

pop_prev_vectors

#ggsave("../AGHY_Final/vital_rate_scripts/Figures/Figures_full_model/pop_prev_vectors.png", pop_prev_vectors)



## find the equilibria points -- where the predictions cross the 1:1 line -- pred_prev_r == x.prev (probably a better way to do this without rounding)
equilibria <- predictions %>% 
  mutate(pred_prev_r = round(pred_prev, digits = 2)) %>% 
  filter(x.prev == pred_prev_r)




eqm_prev <- data.frame(matrix(NA,(AGHY.endochange.20142016$BUGSoutput$n.sims)*2,(N.years-1)))
names(eqm_prev) <- c("y1","y2")
eqm_prev$water <- rep(c("Irrigated","Ambient"),each=AGHY.endochange.20142016$BUGSoutput$n.sims)
#eqm_prev$eqm <- NA

for(i in 1:AGHY.endochange.20142016$BUGSoutput$n.sims){
  for(j in 1:(N.years-1)){
    eqm_prev[i,j] <- x.levels[which.min(abs(x.levels - AGHY.endochange.20142016$BUGSoutput$sims.list$Eplus.control.pred[i,,j]))]
    eqm_prev[i + AGHY.endochange.20142016$BUGSoutput$n.sims,j] <- x.levels[which.min(abs(x.levels - AGHY.endochange.20142016$BUGSoutput$sims.list$Eplus.add.pred[i,,j]))]
  }
}

pop_prev_yr <-eqm_prev%>%
  gather(c("y1","y2"),key="year",value="prevalence")%>%
  mutate(year = as.integer(ifelse(year=="y1",2015,2016)))%>%
  ggplot()+
  geom_histogram(aes(x=prevalence,fill=as.factor(year)),position="dodge",binwidth = 0.01)+
  scale_fill_manual(values = c("firebrick", "black"))+
  facet_grid(water~.)+
  theme_bw()

pop_prev_yr
#ggsave("../AGHY_Final/vital_rate_scripts/Figures/Figures_full_model/pop_prev_yr.png", pop_prev_yr)

pop_prev_water<-eqm_prev%>%
  gather(c("y1","y2"),key="year",value="prevalence")%>%
  mutate(year = as.integer(ifelse(year=="y1",2015,2016)))%>%
  ggplot()+
  geom_histogram(aes(x=prevalence,fill=water),position="dodge",binwidth = 0.01)+
  scale_fill_manual(values = c("orange3", "steelblue4"))+
  facet_grid(as.factor(year)~., scales = "free")+
  theme_bw()

pop_prev_water

#ggsave("../AGHY_Final/vital_rate_scripts/Figures/Figures_full_model/pop_prev_water.png", pop_prev_water)

## do them as densities 
pop_prev_density <-eqm_prev%>%
  gather(c("y1","y2"),key="year",value="prevalence")%>%
  mutate(year = as.integer(ifelse(year=="y1",2015,2016)))%>%
  ggplot()+
  geom_density(aes(x=prevalence,fill=water),alpha = .85)+
  scale_fill_manual(values = c("orange3", "steelblue4"))+
  facet_grid(as.factor(year)~., scales = "free")+
  theme_bw()

pop_prev_density

#ggsave("../AGHY_Final/vital_rate_scripts/Figures/Figures_full_model/pop_prev_density.png", pop_prev_density)



## SURVIVAL 


plot(AGHY.endochange.20142016$BUGSoutput$sims.list$fit.surv,
     AGHY.endochange.20142016$BUGSoutput$sims.list$fit.surv.new,
     xlab="SSQ for actual data",
     ylab="SSQ for perfect (new) data")
abline(0,1)

AGHY.surv.out$BUGSoutput$sims.list$survival.vital.rate ## read about what this output is

bayes.rec2 <-  bayes.data.summary.surv %>% 
  filter(str_detect(info, "^p.rec"))


bayes.surv.rec <- bayes.data.summary.surv %>% 
  filter(str_detect(info, "^ep_rec_plt_t1"))

bayes.surv.rec2 <- bayes.data.summary.surv %>% 
  filter(str_detect(info, "ep_seeds_plt_t"))

bayes.surv.rec21 <- bayes.data.summary.surv %>% 
  filter(str_detect(info, "^ep_rec_flow_t"))

bayes.rec_rate <- bayes.data.summary.surv %>% 
  filter(str_detect(info, "^prob_ep_rec_new"))

bayes.rec_rate_em <- bayes.data.summary.surv %>% 
  filter(str_detect(info, "^prob_em_rec_new"))

bayes.rec2 <- bayes.rec_rate %>% 
  bind_rows(bayes.rec_rate_em)

lambda <- bayes.data.summary.surv %>% 
  filter(str_detect(info, "^mean.vtrans")) 


ggplot(bayes.rec2, aes(info, mean, colour = info))+
  geom_point(size=3, stat="identity")



bayes.surv <-  bayes.data.summary.surv %>% 
  filter(str_detect(info, "^survival.vital.rate")) %>% ## starts with seedmass.per.cap to get the derived seed mass per capita estimates, which are indexed over endo, water, and year
  mutate(endo = case_when(info == "survival.vital.rate[1,1,1]" ~ "Negative", ## [i,j,k], where i is endo status (1 = Negative, 2 = Positive)
                          info == "survival.vital.rate[2,1,1]" ~ "Positive" ,
                          info == "survival.vital.rate[1,2,1]" ~ "Negative",
                          info == "survival.vital.rate[2,2,1]" ~ "Positive" ,
                          info == "survival.vital.rate[1,1,2]" ~ "Negative",
                          info == "survival.vital.rate[2,1,2]" ~ "Positive",
                          info == "survival.vital.rate[1,2,2]" ~ "Negative",
                          info == "survival.vital.rate[2,2,2]" ~ "Positive"),
         
         
         water = case_when(info == "survival.vital.rate[1,1,1]" ~ "Irrigated", ## [i,j,k], where j is water status (1 = Add, 2 = Control) - going to call it 1 = Irrigated, 2 = Ambient
                           info == "survival.vital.rate[2,1,1]" ~ "Irrigated" ,
                           info == "survival.vital.rate[1,2,1]" ~ "Ambient",
                           info == "survival.vital.rate[2,2,1]" ~ "Ambient" ,
                           info == "survival.vital.rate[1,1,2]" ~ "Irrigated",
                           info == "survival.vital.rate[2,1,2]" ~ "Irrigated",
                           info == "survival.vital.rate[1,2,2]" ~ "Ambient",
                           info == "survival.vital.rate[2,2,2]" ~ "Ambient"),
         
         
         year = case_when(info == "survival.vital.rate[1,1,1]" ~ "2015", ## [i,j,k], where k is year (1 = 2014, 2 = 2015)
                          info == "survival.vital.rate[2,1,1]" ~ "2015",
                          info == "survival.vital.rate[1,2,1]" ~ "2015",
                          info == "survival.vital.rate[2,2,1]" ~ "2015",
                          info == "survival.vital.rate[1,1,2]" ~ "2016",
                          info == "survival.vital.rate[2,1,2]" ~ "2016",
                          info == "survival.vital.rate[1,2,2]" ~ "2016",
                          info == "survival.vital.rate[2,2,2]" ~ "2016")) %>% 
  unite(water_endo, water, endo, remove = FALSE) 



survival_plot_scales_free<- ggplot(bayes.surv, aes(endo, mean, colour = water))+
  geom_point(size=3, stat="identity")+
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=.5, size =1)+
  scale_colour_manual(values = c("orange3","steelblue4"))+
  #scale_x_discrete("Endophyte Status", labels = c("Irrigated _ 0" = "E-","Ambient _ 0" = "E-",
  #                                               "Irrigated _ 1" = "E+","Ambient _ 1" = "E+"))+
  labs(title = "Survival Probability", y = "Probability of survival", x = "Endophyte status")+
  guides(color=guide_legend("Watering Regime"))+
  theme(text = element_text(size=18))+
  theme_bw() + #+ theme(panel.border = element_blank()) +
  facet_grid(year~water, scales = "free")

survival_prob_single_model <- ggplot(bayes.surv, aes(water, mean, colour = water, shape = endo))+
  geom_point(size=3, stat="identity", position = position_dodge(width = 1))+
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=.5, size =1,  position = position_dodge(width = 1))+
  scale_colour_manual(values = c("orange3","steelblue4"))+
  #scale_x_discrete("Endophyte Status", labels = c("Irrigated _ 0" = "E-","Ambient _ 0" = "E-",
  #                                               "Irrigated _ 1" = "E+","Ambient _ 1" = "E+"))+
  labs(title = "Survival Probability", y = "Probability of survival", x = "Watering regime")+
  guides(color=guide_legend("Watering Regime"))+
  theme(text = element_text(size=18))+
  theme_bw() + #+ theme(panel.border = element_blank()) +
  facet_grid(~year)+
  ylim(0,1)


survival_prob_single_model




bayes.surv_matched <- (bayes.surv) %>% 
  unite(water_year_endo, water, year, endo, sep = "_", remove = F) %>% 
  mutate(water_year_endo_for_legend = (case_when(water_year_endo == "Ambient_2015_Negative" ~ "Ambient E- '14 - '15",
                                                 water_year_endo == "Ambient_2015_Positive" ~ "Ambient E+ '14 - '15",
                                                 water_year_endo == "Ambient_2016_Negative" ~ "Ambient E- '15 - '16",
                                                 water_year_endo == "Ambient_2016_Positive" ~ "Ambient E+ '15 - '16",
                                                 water_year_endo == "Irrigated_2015_Negative" ~ "Irrigated E- '14 - '15",
                                                 water_year_endo == "Irrigated_2015_Positive" ~ "Irrigated E+ '14 - '15",
                                                 water_year_endo == "Irrigated_2016_Negative" ~ "Irrigated E- '15 - '16",
                                                 water_year_endo == "Irrigated_2016_Positive" ~ "Irrigated E+ '15 - '16")),
         
         water_endo_for_legend = case_when(water_year_endo == "Ambient_2015_Negative" ~ "Ambient E-",
                                           water_year_endo == "Ambient_2015_Positive" ~ "Ambient E+",
                                           water_year_endo == "Ambient_2016_Negative" ~ "Ambient E-",
                                           water_year_endo == "Ambient_2016_Positive" ~ "Ambient E+",
                                           water_year_endo == "Irrigated_2015_Negative" ~ "Irrigated E-",
                                           water_year_endo == "Irrigated_2015_Positive" ~ "Irrigated E+",
                                           water_year_endo == "Irrigated_2016_Negative" ~ "Irrigated E-",
                                           water_year_endo == "Irrigated_2016_Positive" ~ "Irrigated E+"),
        # water = case_when(water == "1" ~ "Ambient",
         #                  water == "2" ~ "Irrigated"),
         trans_year = case_when(year == "2015" ~ "'14 - '15",
                                year == "2016" ~ "'15 - '16"))


survival_match_pop_prev <- ggplot(bayes.surv_matched, aes(trans_year, mean, colour = water_year_endo_for_legend, shape = water_year_endo_for_legend))+
  geom_point(size=3, stat="identity", position = position_dodge(width = 1))+
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=.5, size =1,  position = position_dodge(width = 1))+
  scale_colour_manual(name = "",values = c("#ffae19","#dd7509",
                                           "#ffae19",  "#dd7509",
                                           "#68BAE9", "#0338a8", 
                                           "#68BAE9", "#0338a8"))+
  scale_shape_manual("", values = c(17, 17,
                                    16,16,
                                    17,17,
                                    16,16))+
  #scale_x_discrete("Endophyte Status", labels = c("Irrigated _ 0" = "E-","Ambient _ 0" = "E-",
  #                                               "Irrigated _ 1" = "E+","Ambient _ 1" = "E+"))+
  labs(title = "Survival Probability", y = "Probability of survival", x = "Transition year")+
  # guides(color=guide_legend("Watering Regime and Year"))+
  theme(text = element_text(size=18))+
  theme_classic() + #+ theme(panel.border = element_blank()) +
  facet_grid(~water)#+
#ylim(0,1)



#ggsave("../AGHY_Final/vital_rate_scripts/Figures/Figures_full_model/survival_prob_single_model.png", survival_prob_single_model, width = 8.5, height = 9.5, units = "in")



## flowering


bayes.flow <-  bayes.data.summary.surv %>% 
  filter(str_detect(info, "^flowering.vital.rate")) %>% ## starts with seedmass.per.cap to get the derived seed mass per capita estimates, which are indexed over endo, water, and year
  mutate(endo = case_when(info == "flowering.vital.rate[1,1,1]" ~ "Negative", ## [i,j,k], where i is endo status (1 = Negative, 2 = Positive)
                          info == "flowering.vital.rate[2,1,1]" ~ "Positive" ,
                          info == "flowering.vital.rate[1,2,1]" ~ "Negative",
                          info == "flowering.vital.rate[2,2,1]" ~ "Positive" ,
                          info == "flowering.vital.rate[1,1,2]" ~ "Negative",
                          info == "flowering.vital.rate[2,1,2]" ~ "Positive",
                          info == "flowering.vital.rate[1,2,2]" ~ "Negative",
                          info == "flowering.vital.rate[2,2,2]" ~ "Positive",
                          info == "flowering.vital.rate[1,1,3]" ~ "Negative",
                          info == "flowering.vital.rate[2,1,3]" ~ "Positive",
                          info == "flowering.vital.rate[1,2,3]" ~ "Negative",
                          info == "flowering.vital.rate[2,2,3]" ~ "Positive"),
         
         water = case_when(info == "flowering.vital.rate[1,1,1]" ~ "Irrigated", ## [i,j,k], where j is water status (1 = Add, 2 = Control) - going to call it 1 = Irrigated, 2 = Ambient
                           info == "flowering.vital.rate[2,1,1]" ~ "Irrigated" ,
                           info == "flowering.vital.rate[1,2,1]" ~ "Ambient",
                           info == "flowering.vital.rate[2,2,1]" ~ "Ambient" ,
                           info == "flowering.vital.rate[1,1,2]" ~ "Irrigated",
                           info == "flowering.vital.rate[2,1,2]" ~ "Irrigated",
                           info == "flowering.vital.rate[1,2,2]" ~ "Ambient",
                           info == "flowering.vital.rate[2,2,2]" ~ "Ambient",
                           info == "flowering.vital.rate[1,1,3]" ~ "Irrigated",
                           info == "flowering.vital.rate[2,1,3]" ~ "Irrigated",
                           info == "flowering.vital.rate[1,2,3]" ~ "Ambient",
                           info == "flowering.vital.rate[2,2,3]" ~ "Ambient"),
         
         year = case_when(info == "flowering.vital.rate[1,1,1]" ~ "2014", ## [i,j,k], where k is year (1 = 2014, 2 = 2015)
                          info == "flowering.vital.rate[2,1,1]" ~ "2014",
                          info == "flowering.vital.rate[1,2,1]" ~ "2014",
                          info == "flowering.vital.rate[2,2,1]" ~ "2014",
                          info == "flowering.vital.rate[1,1,2]" ~ "2015",
                          info == "flowering.vital.rate[2,1,2]" ~ "2015",
                          info == "flowering.vital.rate[1,2,2]" ~ "2015",
                          info == "flowering.vital.rate[2,2,2]" ~ "2015",
                          info == "flowering.vital.rate[1,1,3]" ~ "2016",
                          info == "flowering.vital.rate[2,1,3]" ~ "2016",
                          info == "flowering.vital.rate[1,2,3]" ~ "2016",
                          info == "flowering.vital.rate[2,2,3]" ~ "2016")) %>% 
  unite(water_endo, water, endo, remove = FALSE) 



flowering_plot_scales_free<- ggplot(bayes.flow, aes(endo, mean, colour = water))+
  geom_point(size=3, stat="identity")+
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=.5, size =1)+
  scale_colour_manual(values = c("orange3","steelblue4"))+
  #scale_x_discrete("Endophyte Status", labels = c("Irrigated _ 0" = "E-","Ambient _ 0" = "E-",
  #                                               "Irrigated _ 1" = "E+","Ambient _ 1" = "E+"))+
  labs(title = "Flowering Probability", y = "Probability of flowering", x = "Endophyte status")+
  guides(color=guide_legend("Watering Regime"))+
  theme(text = element_text(size=18))+
  theme_bw() + #+ theme(panel.border = element_blank()) +
  facet_grid(year~water, scales = "free")

flowering_prob_single_model <- ggplot(bayes.flow, aes(water, mean, colour = water, shape = endo))+
  geom_point(size=3, stat="identity", position = position_dodge(width = 1))+
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=.5, size =1,  position = position_dodge(width = 1))+
  scale_colour_manual(values = c("orange3","steelblue4"))+
  #scale_x_discrete("Endophyte Status", labels = c("Irrigated _ 0" = "E-","Ambient _ 0" = "E-",
  #                                               "Irrigated _ 1" = "E+","Ambient _ 1" = "E+"))+
  labs(title = "Flowering Probability", y = "Probability of flowering", x = "Watering regime")+
  guides(color=guide_legend("Watering Regime"))+
  theme(text = element_text(size=18))+
  theme_bw() + #+ theme(panel.border = element_blank()) +
  facet_grid(~year)+
  ylim(0,1)


flowering_prob_single_model



bayes.flow_matched <- (bayes.flow) %>% 
  unite(water_year_endo, water, year, endo, sep = "_", remove = F) %>% 
  mutate(water_year_endo_for_legend = (case_when(water_year_endo == "Ambient_2014_Negative" ~ "Ambient E- '14 - '15",
                                                 water_year_endo == "Ambient_2014_Positive" ~ "Ambient E+ '14 - '15",
                                                 water_year_endo == "Ambient_2015_Negative" ~ "Ambient E- '15 - '16",
                                                 water_year_endo == "Ambient_2015_Positive" ~ "Ambient E+ '15 - '16",
                                                 water_year_endo == "Irrigated_2014_Negative" ~ "Irrigated E- '14 - '15",
                                                 water_year_endo == "Irrigated_2014_Positive" ~ "Irrigated E+ '14 - '15",
                                                 water_year_endo == "Irrigated_2015_Negative" ~ "Irrigated E- '15 - '16",
                                                 water_year_endo == "Irrigated_2015_Positive" ~ "Irrigated E+ '15 - '16")),
         
         water_endo_for_legend = case_when(water_year_endo == "Ambient_2014_Negative" ~ "Ambient E-",
                                           water_year_endo == "Ambient_2014_Positive" ~ "Ambient E+",
                                           water_year_endo == "Ambient_2015_Negative" ~ "Ambient E-",
                                           water_year_endo == "Ambient_2015_Positive" ~ "Ambient E+",
                                           water_year_endo == "Irrigated_2014_Negative" ~ "Irrigated E-",
                                           water_year_endo == "Irrigated_2014_Positive" ~ "Irrigated E+",
                                           water_year_endo == "Irrigated_2015_Negative" ~ "Irrigated E-",
                                           water_year_endo == "Irrigated_2015_Positive" ~ "Irrigated E+"),
         # water = case_when(water == "1" ~ "Ambient",
         #                  water == "2" ~ "Irrigated"),
         trans_year = case_when(year == "2014" ~ "'14 - '15",
                                year == "2015" ~ "'15 - '16"))




flow_match_pop_prev <- ggplot(bayes.flow_matched, aes(trans_year, mean, colour = water_year_endo_for_legend, shape = water_year_endo_for_legend))+
  geom_point(size=3, stat="identity", position = position_dodge(width = 1))+
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=.5, size =1,  position = position_dodge(width = 1))+
  scale_colour_manual(name = "",values = c("#ffae19","#dd7509",
                                           "#ffae19",  "#dd7509",
                                           "#68BAE9", "#0338a8", 
                                           "#68BAE9", "#0338a8"))+
  scale_shape_manual("", values = c(17, 17,
                                    16,16,
                                    17,17,
                                    16,16))+
  #scale_x_discrete("Endophyte Status", labels = c("Irrigated _ 0" = "E-","Ambient _ 0" = "E-",
  #                                               "Irrigated _ 1" = "E+","Ambient _ 1" = "E+"))+
  labs(title = "Flowering Probability", y = "Probability of flowering", x = "Transition year")+
  # guides(color=guide_legend("Watering Regime and Year"))+
  theme(text = element_text(size=18))+
  theme_classic() + #+ theme(panel.border = element_blank()) +
  facet_grid(~water)+
  theme(legend.position = "bottom")

flow_match_pop_prev


#ggsave("../AGHY_Final/vital_rate_scripts/Figures/Figures_full_model/flowering_prob_single_model.png", flowering_prob_single_model)  

#bayes.linreg <- bayes.linreg %>% 
#  mutate(exp = exp(mean))

bayes.linreg <-  bayes.data.summary.surv %>% 
  filter (str_detect(info, "^slope_lin_reg"))


bayes.flow.percap <-  bayes.data.summary.surv %>% 
  filter(str_detect(info, "^seedmass.per.cap")) %>% ## starts with seedmass.per.cap to get the derived seed mass per capita estimates, which are indexed over endo, water, and year
  mutate(endo = case_when(info == "seedmass.per.cap[1,1,1]" ~ "Negative", ## [i,j,k], where i is endo status (1 = Negative, 2 = Positive)
                          info == "seedmass.per.cap[2,1,1]" ~ "Positive" ,
                          info == "seedmass.per.cap[1,2,1]" ~ "Negative",
                          info == "seedmass.per.cap[2,2,1]" ~ "Positive" ,
                          info == "seedmass.per.cap[1,1,2]" ~ "Negative",
                          info == "seedmass.per.cap[2,1,2]" ~ "Positive",
                          info == "seedmass.per.cap[1,2,2]" ~ "Negative",
                          info == "seedmass.per.cap[2,2,2]" ~ "Positive"),
         
         water = case_when(info == "seedmass.per.cap[1,1,1]" ~ "Irrigated", ## [i,j,k], where j is water status (1 = Add, 2 = Control) - going to call it 1 = Irrigated, 2 = Ambient
                           info == "seedmass.per.cap[2,1,1]" ~ "Irrigated" ,
                           info == "seedmass.per.cap[1,2,1]" ~ "Ambient",
                           info == "seedmass.per.cap[2,2,1]" ~ "Ambient" ,
                           info == "seedmass.per.cap[1,1,2]" ~ "Irrigated",
                           info == "seedmass.per.cap[2,1,2]" ~ "Irrigated",
                           info == "seedmass.per.cap[1,2,2]" ~ "Ambient",
                           info == "seedmass.per.cap[2,2,2]" ~ "Ambient"),
         
         year = case_when(info == "seedmass.per.cap[1,1,1]" ~ "2014", ## [i,j,k], where k is year (1 = 2014, 2 = 2015)
                          info == "seedmass.per.cap[2,1,1]" ~ "2014",
                          info == "seedmass.per.cap[1,2,1]" ~ "2014",
                          info == "seedmass.per.cap[2,2,1]" ~ "2014",
                          info == "seedmass.per.cap[1,1,2]" ~ "2015",
                          info == "seedmass.per.cap[2,1,2]" ~ "2015",
                          info == "seedmass.per.cap[1,2,2]" ~ "2015",
                          info == "seedmass.per.cap[2,2,2]" ~ "2015")) %>% 
  unite(water_endo, water, endo, remove = FALSE) 


per_capita_seedmass_plot_scales_free<- ggplot(bayes.flow.percap, aes(endo, mean, colour = water))+
  geom_point(size=3, stat="identity")+
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=.5, size =1)+
  scale_colour_manual(values = c("orange3","steelblue4"))+
  #scale_x_discrete("Endophyte Status", labels = c("Irrigated _ 0" = "E-","Ambient _ 0" = "E-",
  #                                               "Irrigated _ 1" = "E+","Ambient _ 1" = "E+"))+
  labs(title = "Per capita seed mass (g)", y = "Per capita seed mass", x = "Endophyte status")+
  guides(color=guide_legend("Watering Regime"))+
  theme(text = element_text(size=18))+
  theme_bw() + #+ theme(panel.border = element_blank()) +
  facet_grid(year~water, scales = "free")

seedmass_percap_single_model <- ggplot(bayes.flow.percap, aes(endo, mean, colour = water))+
  geom_point(size=3, stat="identity")+
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=.5, size =1)+
  scale_colour_manual(values = c("orange3","steelblue4"))+
  #scale_x_discrete("Endophyte Status", labels = c("Irrigated _ 0" = "E-","Ambient _ 0" = "E-",
  #                                               "Irrigated _ 1" = "E+","Ambient _ 1" = "E+"))+
  labs(title = "Per capita seed mass (g)", y = "Per capita seed mass", x = "Endophyte status")+
  guides(color=guide_legend("Watering Regime"))+
  theme(text = element_text(size=18))+
  theme_bw() + #+ theme(panel.border = element_blank()) +
  facet_grid(year~water)


## FIGURE: Per capita seed mass 
seedmass_percap_single_model

#ggsave("../AGHY_Final/vital_rate_scripts/Figures/Figures_full_model/seedmass_percap_single_model.png", seedmass_percap_single_model)  


bayes.percap_matched <- (bayes.flow.percap) %>% 
  unite(water_year_endo, water, year, endo, sep = "_", remove = F) %>% 
  mutate(water_year_endo_for_legend = (case_when(water_year_endo == "Ambient_2014_Negative" ~ "Ambient E- '14 - '15",
                                                 water_year_endo == "Ambient_2014_Positive" ~ "Ambient E+ '14 - '15",
                                                 water_year_endo == "Ambient_2015_Negative" ~ "Ambient E- '15 - '16",
                                                 water_year_endo == "Ambient_2015_Positive" ~ "Ambient E+ '15 - '16",
                                                 water_year_endo == "Irrigated_2014_Negative" ~ "Irrigated E- '14 - '15",
                                                 water_year_endo == "Irrigated_2014_Positive" ~ "Irrigated E+ '14 - '15",
                                                 water_year_endo == "Irrigated_2015_Negative" ~ "Irrigated E- '15 - '16",
                                                 water_year_endo == "Irrigated_2015_Positive" ~ "Irrigated E+ '15 - '16")),
         
         water_endo_for_legend = case_when(water_year_endo == "Ambient_2014_Negative" ~ "Ambient E-",
                                           water_year_endo == "Ambient_2014_Positive" ~ "Ambient E+",
                                           water_year_endo == "Ambient_2015_Negative" ~ "Ambient E-",
                                           water_year_endo == "Ambient_2015_Positive" ~ "Ambient E+",
                                           water_year_endo == "Irrigated_2014_Negative" ~ "Irrigated E-",
                                           water_year_endo == "Irrigated_2014_Positive" ~ "Irrigated E+",
                                           water_year_endo == "Irrigated_2015_Negative" ~ "Irrigated E-",
                                           water_year_endo == "Irrigated_2015_Positive" ~ "Irrigated E+"),
         # water = case_when(water == "1" ~ "Ambient",
         #                  water == "2" ~ "Irrigated"),
         trans_year = case_when(year == "2014" ~ "'14 - '15",
                                year == "2015" ~ "'15 - '16"))


flow_percap_match_pop_prev <- ggplot(bayes.percap_matched, aes(trans_year, mean, colour = water_year_endo_for_legend, shape = water_year_endo_for_legend))+
  geom_point(size=3, stat="identity", position = position_dodge(width = 1))+
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=.5, size =1,  position = position_dodge(width = 1))+
  scale_colour_manual(name = "",values = c("#ffae19","#dd7509",
                                           "#ffae19",  "#dd7509",
                                           "#68BAE9", "#0338a8", 
                                           "#68BAE9", "#0338a8"))+
  scale_shape_manual("", values = c(17, 17,
                                    16,16,
                                    17,17,
                                    16,16))+
  #scale_x_discrete("Endophyte Status", labels = c("Irrigated _ 0" = "E-","Ambient _ 0" = "E-",
  #                                               "Irrigated _ 1" = "E+","Ambient _ 1" = "E+"))+
  labs(title = "Per capita seed mass", y = "Per capita seed mass (g)", x = "Transition year")+
  # guides(color=guide_legend("Watering Regime and Year"))+
  theme(text = element_text(size=18))+
  theme_classic() + #+ theme(panel.border = element_blank()) +
  facet_grid(~water, scales= "free")+
  theme(legend.position = "bottom")

flow_percap_match_pop_prev

### components of fertility


bayes.flow.percap.mu <-  bayes.data.summary.surv %>% 
  filter(str_detect(info, "^mu.")) %>% ## starts with seedmass.per.cap to get the derived seed mass per capita estimates, which are indexed over endo, water, and year
  mutate(endo = case_when(info == "mu.infs[1,1,1]" ~ "Negative", ## [i,j,k], where k is year (1 = 2014, 2 = 2015)
                          info == "mu.infs[2,1,1]" ~ "Positive",
                          info == "mu.infs[1,2,1]" ~ "Negative",
                          info == "mu.infs[2,2,1]" ~ "Positive",
                          info == "mu.infs[1,1,2]" ~ "Negative", ## [i,j,k], where k is year (1 = 2014, 2 = 2015)
                          info == "mu.infs[2,1,2]" ~ "Positive",
                          info == "mu.infs[1,2,2]" ~ "Negative",
                          info == "mu.infs[2,2,2]" ~ "Positive",
                          info == "mu.seed[1,1,1]" ~ "Negative", ## [i,j,k], where k is year (1 = 2014, 2 = 2015)
                          info == "mu.seed[2,1,1]" ~ "Positive",
                          info == "mu.seed[1,2,1]" ~ "Negative",
                          info == "mu.seed[2,2,1]" ~ "Positive",
                          info == "mu.seed[1,1,2]" ~ "Negative", ## [i,j,k], where k is year (1 = 2014, 2 = 2015)
                          info == "mu.seed[2,1,2]" ~ "Positive",
                          info == "mu.seed[1,2,2]" ~ "Negative",
                          info == "mu.seed[2,2,2]" ~ "Positive"),
         
         water = case_when(info == "mu.infs[1,1,1]" ~ "Irrigated", ## [i,j,k], where k is year (1 = 2014, 2 = 2015)
                           info == "mu.infs[2,1,1]" ~ "Irrigated",
                           info == "mu.infs[1,2,1]" ~ "Ambient",
                           info == "mu.infs[2,2,1]" ~ "Ambient",
                           info == "mu.infs[1,1,2]" ~ "Irrigated", ## [i,j,k], where k is year (1 = 2014, 2 = 2015)
                           info == "mu.infs[2,1,2]" ~ "Irrigated",
                           info == "mu.infs[1,2,2]" ~ "Ambient",
                           info == "mu.infs[2,2,2]" ~ "Ambient",
                           info == "mu.seed[1,1,1]" ~ "Irrigated", ## [i,j,k], where k is year (1 = 2014, 2 = 2015)
                           info == "mu.seed[2,1,1]" ~ "Irrigated",
                           info == "mu.seed[1,2,1]" ~ "Ambient",
                           info == "mu.seed[2,2,1]" ~ "Ambient",
                           info == "mu.seed[1,1,2]" ~ "Irrigated", ## [i,j,k], where k is year (1 = 2014, 2 = 2015)
                           info == "mu.seed[2,1,2]" ~ "Irrigated",
                           info == "mu.seed[1,2,2]" ~ "Ambient",
                           info == "mu.seed[2,2,2]" ~ "Ambient"),
         
         type_mean = case_when(info == "mu.infs[1,1,1]" ~ "infs", ## [i,j,k], where k is year (1 = 2014, 2 = 2015)
                               info == "mu.infs[2,1,1]" ~ "infs",
                               info == "mu.infs[1,2,1]" ~ "infs",
                               info == "mu.infs[2,2,1]" ~ "infs",
                               info == "mu.infs[1,1,2]" ~ "infs", ## [i,j,k], where k is year (1 = 2014, 2 = 2015)
                               info == "mu.infs[2,1,2]" ~ "infs",
                               info == "mu.infs[1,2,2]" ~ "infs",
                               info == "mu.infs[2,2,2]" ~ "infs",
                               info == "mu.seed[1,1,1]" ~ "seeds", ## [i,j,k], where k is year (1 = 2014, 2 = 2015)
                               info == "mu.seed[2,1,1]" ~ "seeds",
                               info == "mu.seed[1,2,1]" ~ "seeds",
                               info == "mu.seed[2,2,1]" ~ "seeds",
                               info == "mu.seed[1,1,2]" ~ "seeds", ## [i,j,k], where k is year (1 = 2014, 2 = 2015)
                               info == "mu.seed[2,1,2]" ~ "seeds",
                               info == "mu.seed[1,2,2]" ~ "seeds",
                               info == "mu.seed[2,2,2]" ~ "seeds"),
         year = case_when(info == "mu.infs[1,1,1]" ~ "2014", ## [i,j,k], where k is year (1 = 2014, 2 = 2015)
                          info == "mu.infs[2,1,1]" ~ "2014",
                          info == "mu.infs[1,2,1]" ~ "2014",
                          info == "mu.infs[2,2,1]" ~ "2014",
                          info == "mu.infs[1,1,2]" ~ "2015", ## [i,j,k], where k is year (1 = 2014, 2 = 2014)
                          info == "mu.infs[2,1,2]" ~ "2015",
                          info == "mu.infs[1,2,2]" ~ "2015",
                          info == "mu.infs[2,2,2]" ~ "2015",
                          info == "mu.seed[1,1,1]" ~ "2014", ## [i,j,k], where k is year (1 = 2014, 2 = 2014)
                          info == "mu.seed[2,1,1]" ~ "2014",
                          info == "mu.seed[1,2,1]" ~ "2014",
                          info == "mu.seed[2,2,1]" ~ "2014",
                          info == "mu.seed[1,1,2]" ~ "2015", ## [i,j,k], where k is year (1 = 2014, 2 = 2015)
                          info == "mu.seed[2,1,2]" ~ "2015",
                          info == "mu.seed[1,2,2]" ~ "2015",
                          info == "mu.seed[2,2,2]" ~ "2015"))%>% 
  #unite(water_endo = c(water, endo),remove = FALSE) %>% 
  unite(type_year, type_mean, year,remove = FALSE) %>% 
  #rename(water_endo = names(.)[11]) %>% 
  dplyr::rename(type_year = names(.)[13]) %>% 
  mutate(means_real = exp(mean),
         CI_low = exp(`2.5%`),
         CI_high = exp(`97.5%`)) %>% 
  filter(type_year != "NA_NA")











components_fert<- ggplot(bayes.flow.percap.mu, aes(endo, means_real, colour = water))+
  geom_point(size=3, stat="identity", position = position_dodge(width = .75))+
  geom_errorbar(aes(ymin=CI_low, ymax=CI_high), width=.5, size =1, position = position_dodge(width = .75))+
  scale_colour_manual(values = c("orange3","steelblue4"))+
  #scale_x_discrete("Endophyte Status", labels = c("Irrigated _ 0" = "E-","Ambient _ 0" = "E-",
  #                                               "Irrigated _ 1" = "E+","Ambient _ 1" = "E+"))+
  labs(title = "Components of fertility", y = "Estimated per capita seed mass (g) and number of inflorescences", x = "Endophyte status")+
  guides(color=guide_legend("Watering Regime"))+
  theme(text = element_text(size=18))+
  theme_bw() + 
  facet_grid(type_mean ~ year, scales = "free")

components_fert


#ggsave("../AGHY_Final/vital_rate_scripts/Figures/components_fertNov8_YES.png", components_fert)  

bayes.comp_fert_matched <- (bayes.flow.percap.mu) %>% 
  unite(water_year_endo, water, year, endo, sep = "_", remove = F) %>% 
  mutate(water_year_endo_for_legend = (case_when(water_year_endo == "Ambient_2014_Negative" ~ "Ambient E- '14 - '15",
                                                 water_year_endo == "Ambient_2014_Positive" ~ "Ambient E+ '14 - '15",
                                                 water_year_endo == "Ambient_2015_Negative" ~ "Ambient E- '15 - '16",
                                                 water_year_endo == "Ambient_2015_Positive" ~ "Ambient E+ '15 - '16",
                                                 water_year_endo == "Irrigated_2014_Negative" ~ "Irrigated E- '14 - '15",
                                                 water_year_endo == "Irrigated_2014_Positive" ~ "Irrigated E+ '14 - '15",
                                                 water_year_endo == "Irrigated_2015_Negative" ~ "Irrigated E- '15 - '16",
                                                 water_year_endo == "Irrigated_2015_Positive" ~ "Irrigated E+ '15 - '16")),
         
         water_endo_for_legend = case_when(water_year_endo == "Ambient_2014_Negative" ~ "Ambient E-",
                                           water_year_endo == "Ambient_2014_Positive" ~ "Ambient E+",
                                           water_year_endo == "Ambient_2015_Negative" ~ "Ambient E-",
                                           water_year_endo == "Ambient_2015_Positive" ~ "Ambient E+",
                                           water_year_endo == "Irrigated_2014_Negative" ~ "Irrigated E-",
                                           water_year_endo == "Irrigated_2014_Positive" ~ "Irrigated E+",
                                           water_year_endo == "Irrigated_2015_Negative" ~ "Irrigated E-",
                                           water_year_endo == "Irrigated_2015_Positive" ~ "Irrigated E+"),
         # water = case_when(water == "1" ~ "Ambient",
         #                  water == "2" ~ "Irrigated"),
         trans_year = case_when(year == "2014" ~ "'14 - '15",
                                year == "2015" ~ "'15 - '16"))


comp_fert_match_pop_prev <- ggplot(bayes.comp_fert_matched, aes(trans_year, means_real, colour = water_year_endo_for_legend, shape = water_year_endo_for_legend))+
  geom_point(size=3, stat="identity", position = position_dodge(width = 1))+
  geom_errorbar(aes(ymin=CI_low, ymax=CI_high), width=.5, size =1,  position = position_dodge(width = 1))+
  scale_colour_manual(name = "",values = c("#ffae19","#dd7509",
                                           "#ffae19",  "#dd7509",
                                           "#68BAE9", "#0338a8", 
                                           "#68BAE9", "#0338a8"))+
  scale_shape_manual("", values = c(17, 17,
                                    16,16,
                                    17,17,
                                    16,16))+
  #scale_x_discrete("Endophyte Status", labels = c("Irrigated _ 0" = "E-","Ambient _ 0" = "E-",
  #                                               "Irrigated _ 1" = "E+","Ambient _ 1" = "E+"))+
  labs(title = "Per capita fertility components", y = "Estimated per capita number of inflorescences and seed mass (g)", x = "Transition year")+
  # guides(color=guide_legend("Watering Regime and Year"))+
  theme(text = element_text(size=18))+
  theme_classic() + #+ theme(panel.border = element_blank()) +
  facet_grid(type_mean~water, scales= "free")+
  theme(legend.position = "bottom")

comp_fert_match_pop_prev

### RECRUITMENT

## get plot info to match with probability of recruitment output
plot.info.ep <- rec_dat %>% 
  select(water, plot, year) %>% 
  mutate(endo = 2)

plot.info <- rec_dat %>% 
  select(water, plot, year) %>% 
  mutate(endo = 1) %>% 
  bind_rows(plot.info.ep)


recruits_EP <- bayes.data.summary.surv %>% 
  filter(str_detect(info, "^ep_"))


bayes.rec_new <- bayes.data.summary.surv %>% 
  filter(str_detect(info, "p.rec"))

bayes.rec <-  bayes.data.summary.surv %>% 
  filter(str_detect(info, "^prob_ep_rec")|
           str_detect(info, "^prob_em_rec")) %>% ## starts with seedmass.per.cap to get the derived seed mass per capita estimates, which are indexed over endo, water, and year
  rownames_to_column(var = "rowname") #%>%
  #bind_cols(plot.info) %>% 
#  mutate(endo = as.factor(endo),
 #        water = as.factor(water),
  #       year = case_when(year == "1" ~ "2014",
   #                       year == "2" ~ "2015"))
#rec_plot_dat <-  bayes.rec %>%
#  group_by(endo, year, water) %>%
#  rename(prob = mean) %>% 
##  summarise(mean = mean(prob, na.rm = TRUE),
#            sd = sd(prob, na.rm = TRUE),
#            n = n()) %>%
#  mutate(se = sd / sqrt(n),
#         lower.ci = mean - qt(1 - (0.05 / 2), n - 1) * se,
#         upper.ci = mean + qt(1 - (0.05 / 2), n - 1) * se)


#ep_rec_plt
#recruits <- bayes.data.summary.surv %>% 
#  filter(str_detect(info, "_rec_plt"))

#recruits1 <- bayes.data.summary.surv %>% 
#  filter(str_detect(info, "avg.rec"))

bayes.rec.matched <- (bayes.rec) %>% 
  separate(info, c("blah","endo", "type", "blah2","water", "year"), remove = F) %>% 
  unite(water_year_endo, c("water","year","endo"), sep = "_", remove = F) %>% 
  mutate(water_year_endo_for_legend = (case_when(water_year_endo == "ambient_14_em" ~ "Ambient E- '14 - '15",
                                                 water_year_endo == "ambient_14_ep" ~ "Ambient E+ '14 - '15",
                                                 water_year_endo == "ambient_15_em" ~ "Ambient E- '15 - '16",
                                                 water_year_endo == "ambient_15_ep" ~ "Ambient E+ '15 - '16",
                                                 water_year_endo == "irrigated_14_em" ~ "Irrigated E- '14 - '15",
                                                 water_year_endo == "irrigated_14_ep" ~ "Irrigated E+ '14 - '15",
                                                 water_year_endo == "irrigated_15_em" ~ "Irrigated E- '15 - '16",
                                                 water_year_endo == "irrigated_15_ep" ~ "Irrigated E+ '15 - '16")),
         
         water_endo_for_legend = case_when(water_year_endo == "ambient_14_em" ~ "Ambient E-",
                                           water_year_endo == "ambient_14_ep" ~ "Ambient E+",
                                           water_year_endo == "ambient_15_em" ~ "Ambient E-",
                                           water_year_endo == "ambient_15_ep" ~ "Ambient E+",
                                           water_year_endo == "irrigated_14_em" ~ "Irrigated E-",
                                           water_year_endo == "irrigated_14_ep" ~ "Irrigated E+",
                                           water_year_endo == "irrigated_15_em" ~ "Irrigated E-",
                                           water_year_endo == "irrigated_15_ep" ~ "Irrigated E+"),
         # water = case_when(water == "1" ~ "Ambient",
         #                  water == "2" ~ "Irrigated"),
         trans_year = case_when(year == "14" ~ "'14 - '15",
                                year == "15" ~ "'15 - '16")) %>% 
  mutate(CI_low = (`2.5%`),
         CI_high = (`97.5%`),
         water = (case_when(water == "ambient" ~ "Ambient",
                            water == "irrigated" ~ "Irrigated")))




rec_match_pop_prev <- ggplot(bayes.rec.matched, aes(trans_year, mean, colour = water_year_endo_for_legend, shape = water_year_endo_for_legend))+
  geom_point(size=3, stat="identity", position = position_dodge(width = 1))+
  geom_errorbar(aes(ymin=CI_low, ymax=CI_high), width=.5, size =1,  position = position_dodge(width = 1))+
  scale_colour_manual(name = "",values = c("#ffae19","#dd7509",
                                           "#ffae19",  "#dd7509",
                                           "#68BAE9", "#0338a8", 
                                           "#68BAE9", "#0338a8"))+
  scale_shape_manual("", values = c(17, 17,
                                    16,16,
                                    17,17,
                                    16,16))+
  #scale_x_discrete("Endophyte Status", labels = c("Irrigated _ 0" = "E-","Ambient _ 0" = "E-",
  #                                               "Irrigated _ 1" = "E+","Ambient _ 1" = "E+"))+
  labs(title = "Probability of recruitment", y = "Probability of recruitment", x = "Transition year")+
  # guides(color=guide_legend("Watering Regime and Year"))+
  theme(text = element_text(size=18))+
  theme_classic() + #+ theme(panel.border = element_blank()) +
  facet_grid(~water, scales= "free")+
  theme(legend.position = "bottom")

rec_match_pop_prev

  


############################
### VERTICAL TRANSMISSION ##
############################
bayes.vt <- bayes.data.summary.surv %>%
  filter(str_detect(info, "^vtrans")) %>% 
  mutate(water = case_when(info == "vtrans[1,1]" ~ "Ambient",
                           info == "vtrans[2,1]" ~ "Irrigated",
                           info == "vtrans[1,2]" ~ "Ambient",
                           info == "vtrans[2,2]" ~ "Irrigated"),
         
         year = case_when(info == "vtrans[1,1]" ~ "2014",
                          info == "vtrans[2,1]" ~ "2014",
                          info == "vtrans[1,2]" ~ "2015",
                          info == "vtrans[2,2]" ~ "2015")) %>% 
  unite(water_year, c(water, year), remove = FALSE)
  




## shows probabilty of vertical transmission for recruits in 2014 and 2015
vtrans_year_single_model<- ggplot(bayes.vt, aes(water, mean, colour = water_year))+
  geom_point(size=3, stat="identity")+
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=.5, size =1)+
  #scale_colour_manual(values = c("orange3","steelblue4"))+
  scale_colour_manual(name = "",values = c("#ffae19","#dd7509",
                                           #"#ffae19",  "#dd7509",
                                           #"#68BAE9", "#0338a8", 
                                           "#68BAE9", "#0338a8"))+
  #scale_x_discrete("Endophyte Status", labels = c("Irrigated _ 0" = "E-","Ambient _ 0" = "E-",
  #                                               "Irrigated _ 1" = "E+","Ambient _ 1" = "E+"))+
  labs(title = "Vertical transmission", y = "Vertical tranmsission probability", x = "Watering Regime")+
  scale_y_continuous(limits = c(0,1))+  
  guides(color=guide_legend("Watering Regime"))+
  theme(text = element_text(size=18))+
  theme_bw() + #+ theme(panel.border = element_blank()) +
  facet_grid(~year)

vtrans_year_single_model

#ggsave("../AGHY_Final/vital_rate_scripts/Figures/Figures_full_model/vtrans_year_single_model.png", vtrans_year_single_model) 


## good fit for vtrans posterior predictive check with the beta binomial
plot(AGHY.endochange.20142016$BUGSoutput$sims.list$fit.vtrans,
     AGHY.endochange.20142016$BUGSoutput$sims.list$fit.vtrans.new,
     xlab="SSQ for actual data",
     ylab="SSQ for perfect (new) data")
abline(0,1)





### make a combo figure for the vital rates
bayes.flow_matched <- bayes.flow_matched %>% 
                      mutate(vital_rate = "Flowering")

bayes.surv_matched <- bayes.surv_matched %>% 
                      mutate(vital_rate = "Survival")



rec_matched <- bayes.rec.matched %>% 
               select(water_year_endo, water_year_endo_for_legend, water_endo_for_legend, trans_year, mean, se, water, endo, year) %>% 
               mutate('2.5%' = mean - se,
                      '97.5%' = mean + se,
                      vital_rate = "Recruitment") %>%  ## this isn't equivalent to what the CIs are from the Bayes model - but in order to put them all on the same figure, I'll put them under the same header
              select(-se)                        


combo_fig_df <- bayes.flow_matched %>% 
  bind_rows(bayes.surv_matched, rec_matched)




ggplot(combo_fig_df, aes(trans_year, mean, colour = water_year_endo_for_legend, shape = water_year_endo_for_legend))+
  geom_point(size=3, stat="identity", position = position_dodge(width = 1))+
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=.5, size =1,  position = position_dodge(width = 1))+
  scale_colour_manual(name = "",values = c("#ffae19","#dd7509",
                                           "#ffae19",  "#dd7509",
                                           "#68BAE9", "#0338a8", 
                                           "#68BAE9", "#0338a8"))+
  scale_shape_manual("", values = c(17, 17,
                                    16,16,
                                    17,17,
                                    16,16))+
  #scale_x_discrete("Endophyte Status", labels = c("Irrigated _ 0" = "E-","Ambient _ 0" = "E-",
  #                                               "Irrigated _ 1" = "E+","Ambient _ 1" = "E+"))+
  labs(title = "Vital Rate Probabilities", y = "Vital Rate Probability", x = "Transition year")+
  # guides(color=guide_legend("Watering Regime and Year"))+
  theme(text = element_text(size=18))+
  theme_bw() + #+ theme(panel.border = element_blank()) +
  facet_grid(vital_rate~water, scales = "free")+
  theme(legend.position = "bottom")



### components of fertility

bayes.comp_fert_matched1 <- bayes.comp_fert_matched %>% 
                           select(-"2.5%", -"97.5%") %>% 
                           rename("2.5%" = CI_low,
                                  "97.5%" = CI_high) %>% 
                          mutate(fert_comp = case_when(type_mean == "infs" ~ "Number of inflorescences",
                                                       type_mean == "seeds" ~ "Seed mass (g)") )
 
bayes.percap_matched1 <- bayes.percap_matched %>% 
                         mutate(fert_comp = "Per capita seed mass (g)")

fert_combo_df<-bayes.comp_fert_matched1 %>% 
               bind_rows(bayes.percap_matched1)

ggplot(fert_combo_df, aes(trans_year, means_real, colour = water_year_endo_for_legend, shape = water_year_endo_for_legend))+
  geom_point(size=3, stat="identity", position = position_dodge(width = 1))+
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=.5, size =1,  position = position_dodge(width = 1))+
  scale_colour_manual(name = "",values = c("#ffae19","#dd7509",
                                           "#ffae19",  "#dd7509",
                                           "#68BAE9", "#0338a8", 
                                           "#68BAE9", "#0338a8"))+
  scale_shape_manual("", values = c(17, 17,
                                    16,16,
                                    17,17,
                                    16,16))+
  #scale_x_discrete("Endophyte Status", labels = c("Irrigated _ 0" = "E-","Ambient _ 0" = "E-",
  #                                               "Irrigated _ 1" = "E+","Ambient _ 1" = "E+"))+
  labs(title = "Components of the flowering vital rate", y = "Per capita inflorescence number or seed mass (g)", x = "Transition year")+
  # guides(color=guide_legend("Watering Regime and Year"))+
  theme(text = element_text(size=18))+
  theme_bw() + #+ theme(panel.border = element_blank()) +
  facet_grid(fert_comp~water, scales = "free")+
  theme(legend.position = "bottom") 



##### POSTERIOR PREDICTIVE CHECKS

prev.ppc.df <- as.data.frame(AGHY.endochange.20142016$BUGSoutput$sims.list$fit) %>% 
  dplyr::rename(actual = V1) %>%
           bind_cols(as.data.frame(AGHY.endochange.20142016$BUGSoutput$sims.list$fit.new)) %>% 
           dplyr::rename(fitted = V1) %>% 
           mutate(type = "Population Prevalence",
                  level = 1)

surv.ppc.df <- as.data.frame(AGHY.endochange.20142016$BUGSoutput$sims.list$fit.surv) %>% 
  dplyr::rename(actual = V1) %>%
                bind_cols(as.data.frame(AGHY.endochange.20142016$BUGSoutput$sims.list$fit.surv.new)) %>% 
                dplyr::rename(fitted = V1) %>% 
                mutate(type = "Survival", 
                       level = 2)

seed.ppc.df <- as.data.frame(AGHY.endochange.20142016$BUGSoutput$sims.list$fit.seed) %>% 
  dplyr::rename(actual = V1) %>%
  bind_cols(as.data.frame(AGHY.endochange.20142016$BUGSoutput$sims.list$fit.seed.new)) %>% 
  dplyr::rename(fitted = V1) %>% 
  mutate(type = "Seed",
         level = 3)

infs.ppc.df <- as.data.frame(AGHY.endochange.20142016$BUGSoutput$sims.list$fit.infs) %>% 
  dplyr::rename(actual = V1) %>%
  bind_cols(as.data.frame(AGHY.endochange.20142016$BUGSoutput$sims.list$fit.infs.new)) %>% 
  dplyr::rename(fitted = V1) %>% 
  mutate(type = "Known Inflorescences",
         level = 4)


infs.unk.ppc.df <- as.data.frame(AGHY.endochange.20142016$BUGSoutput$sims.list$fit.infs.unk) %>% 
  dplyr::rename(actual = V1) %>%
  bind_cols(as.data.frame(AGHY.endochange.20142016$BUGSoutput$sims.list$fit.infs.unk.new)) %>% 
  dplyr::rename(fitted = V1) %>% 
  mutate(type = "Unknown Inflorescences",
         level = 5)

linreg.ppc.df <- as.data.frame(AGHY.endochange.20142016$BUGSoutput$sims.list$fit.linreg) %>% 
  dplyr::rename(actual = V1) %>%
  bind_cols(as.data.frame(AGHY.endochange.20142016$BUGSoutput$sims.list$fit.linreg.new)) %>% 
  dplyr::rename(fitted = V1) %>% 
  mutate(type = "Seed mass to count",
         level = 6)

rec.ppc.df <- as.data.frame(AGHY.endochange.20142016$BUGSoutput$sims.list$fit.rec) %>% 
  dplyr::rename(actual = V1) %>%
  bind_cols(as.data.frame(AGHY.endochange.20142016$BUGSoutput$sims.list$fit.rec.new)) %>% 
  dplyr::rename(fitted = V1) %>% 
  mutate(type = "Recruitment",
         level = 7)

flow.ppc.df <- as.data.frame(AGHY.endochange.20142016$BUGSoutput$sims.list$fit.flow) %>% 
  dplyr::rename(actual = V1) %>%
  bind_cols(as.data.frame(AGHY.endochange.20142016$BUGSoutput$sims.list$fit.flow.new)) %>% 
  dplyr::rename(fitted = V1) %>% 
  mutate(type = "Flowering",
         level = 8)

vtrans.ppc.df <- as.data.frame(AGHY.endochange.20142016$BUGSoutput$sims.list$fit.vtrans) %>% 
  dplyr::rename(actual = V1) %>%
  bind_cols(as.data.frame(AGHY.endochange.20142016$BUGSoutput$sims.list$fit.vtrans.new))%>% 
  dplyr::rename(fitted = V1) %>% 
  mutate(type = "Vertical Transmission",
         level = 9)
  
posterior_pred_check_df <- prev.ppc.df %>% 
                          bind_rows(surv.ppc.df, seed.ppc.df, infs.ppc.df,
                                    linreg.ppc.df, rec.ppc.df, flow.ppc.df, vtrans.ppc.df)

ppc_fig <- ggplot(posterior_pred_check_df, aes(actual, fitted))+
  geom_point() +
  geom_abline(color ="gray", linetype ="dashed", size = 1) +
  facet_wrap(type~., scales = "free")+
  theme_bw()

#ggsave("../AGHY_Final/vital_rate_scripts/Figures/Figures_full_model/ppc_fig.png", ppc_fig)  


bayes.fecundity_em <- bayes.data.summary.surv %>% 
  filter(str_detect(info, "^em_seeds_")) %>% 
  rownames_to_column(var = "rowname") 

bayes.fecundity_ep <- bayes.data.summary.surv %>% 
  filter(str_detect(info, "^ep_seeds_")) %>% 
  rownames_to_column(var = "rowname") 

bayes.fecundity <- bayes.fecundity_ep %>% 
  bind_rows(bayes.fecundity_em) %>% 
  separate(info, c("endo","type","water","year"), remove = F) %>% 
  unite(water_year_endo, water, year, endo, sep = "_", remove = F) %>% 
  mutate(year = case_when(year == "14" ~ "2014",
                          year == "15" ~ "2015"),
         water = case_when(water == "ambient" ~ "1",
                           water == "irrigated" ~ "2"),
         endo = case_when(endo == "em" ~ "1",
                          endo == "ep" ~ "2")) %>% 
  unite(water_year_endo, c(water, year, endo), remove = F) %>% 
  mutate(water_year_endo_for_legend = (case_when(water_year_endo == "1_2014_1" ~ "Ambient E- '14 - '15",
                                                 water_year_endo == "1_2014_2" ~ "Ambient E+ '14 - '15",
                                                 water_year_endo == "1_2015_1" ~ "Ambient E- '15 - '16",
                                                 water_year_endo == "1_2015_2" ~ "Ambient E+ '15 - '16",
                                                 water_year_endo == "2_2014_1" ~ "Irrigated E- '14 - '15",
                                                 water_year_endo == "2_2014_2" ~ "Irrigated E+ '14 - '15",
                                                 water_year_endo == "2_2015_1" ~ "Irrigated E- '15 - '16",
                                                 water_year_endo == "2_2015_2" ~ "Irrigated E+ '15 - '16")),
         water_endo_for_legend = case_when(water_year_endo == "1_2014_1" ~ "Ambient E-",
                                           water_year_endo == "1_2014_2" ~ "Ambient E+",
                                           water_year_endo == "1_2015_1" ~ "Ambient E-",
                                           water_year_endo == "1_2015_2" ~ "Ambient E+",
                                           water_year_endo == "2_2014_1" ~ "Irrigated E-",
                                           water_year_endo == "2_2014_2" ~ "Irrigated E+",
                                           water_year_endo == "2_2015_1" ~ "Irrigated E-",
                                           water_year_endo == "2_2015_2" ~ "Irrigated E+"),
         water = case_when(water == "1" ~ "Ambient",
                           water == "2" ~ "Irrigated"),
         trans_year = case_when(year == "2014" ~ "'14 - '15",
                                year == "2015" ~ "'15 - '16"))

ggplot(bayes.fecundity, aes(trans_year, mean, colour = water_year_endo_for_legend, shape = water_year_endo_for_legend))+
  geom_point(size=3, stat="identity", position = position_dodge(width = 1))+
  geom_errorbar(aes(ymin=`2.5%`, ymax=`97.5%`), width=.5, size =1,  position = position_dodge(width = 1))+
  scale_colour_manual(name = "",values = c("#ffae19","#dd7509",
                                           "#ffae19",  "#dd7509",
                                           "#68BAE9", "#0338a8", 
                                           "#68BAE9", "#0338a8"))+
  scale_shape_manual("", values = c(21, 21,
                                    16,16,
                                    21,21,
                                    16,16))+
  #scale_x_discrete("Endophyte Status", labels = c("Irrigated _ 0" = "E-","Ambient _ 0" = "E-",
  #                                               "Irrigated _ 1" = "E+","Ambient _ 1" = "E+"))+
  labs(title = "Fecundity", y = "Fecundity (number of seeds)", x = "Transition year")+
  # guides(color=guide_legend("Watering Regime and Year"))+
  theme(text = element_text(size=18))+
  theme_classic() + #+ theme(panel.border = element_blank()) +
  facet_grid(~water)+
  theme(legend.position = "bottom")



#################################NEW
############################ NEW
########################## NEW
library(tidyverse)

## need to fix this -- survival and flowering are on different "years" -- make it transition year to keep as one df
Filter_list_SF <- c("^survival.vital.rate", 
                    "^flowering.vital.rate")
endo_SF <- rep(c("Negative", "Positive"), each = 1200, times = 4)
water_SF <- rep(c("Elevated", "Ambient"), each = 2400, times = 2)
year_SF <- rep(c('14-15', '15-16'), each = 2400, times = 2)

e_w_y_SF <- data.frame(endo_SF, water_SF, year_SF)

names(e_w_y_SF) <- c("endo", "water", "trans_year")

### new p rec

p_rec_post <- as.data.frame(AGHY.endochange.20142016$BUGSoutput$sims.list) %>% 
  gather("variable", "posterior") %>%
  filter(str_detect(variable, "rec_new")) %>% 
  mutate(vital_rate = "recruitment",
  variable = as.factor(variable)) %>% 
  mutate(water = case_when(str_detect(variable, "ambient") ~ "Ambient",
                           str_detect(variable, "irrigated") ~ "Elevated"),
         endo = case_when(str_detect(variable, "em") ~ "Negative",
                          str_detect(variable, "ep") ~ "Positive"),
         year = case_when(str_detect(variable, "14") ~ "14-15",
                          str_detect(variable, "15") ~ "15-16")) %>% 
  mutate(year = as.factor(year),
         water = as.factor(water)) %>% 
  unite_("water_endo", c("water", "endo"), sep="_", remove=F) %>% 
  unite_("water_yr", c("water", "year"), sep= "_", remove = F) %>% 
  mutate(water_endo_factor = as.factor(water_endo))

posteriors <- as.data.frame(AGHY.endochange.20142016$BUGSoutput$sims.list) %>% 
  gather("variable", "posterior") %>% 
  filter(str_detect(variable, Filter_list_SF)) %>% 
  mutate(vital_rate = case_when(str_detect(variable, "^survival.vital.rate") ~ "survival",
                                str_detect(variable, "^flowering.vital.rate") ~ "flowering"),
         variable = as.factor(variable)) %>% 
  mutate(water = case_when(variable == "survival.vital.rate.1" ~ "Elevated",
                           variable == "survival.vital.rate.2" ~ "Elevated",
                           variable == "survival.vital.rate.3" ~ "Ambient",
                           variable == "survival.vital.rate.4" ~ "Ambient",
                           variable == "survival.vital.rate.5" ~ "Elevated",
                           variable == "survival.vital.rate.6" ~ "Elevated",
                           variable == "survival.vital.rate.7" ~ "Ambient",
                           variable == "survival.vital.rate.8" ~ "Ambient",
                           
                           variable == "flowering.vital.rate.1" ~ "Elevated",
                           variable == "flowering.vital.rate.2" ~ "Elevated",
                           variable == "flowering.vital.rate.3" ~ "Ambient",
                           variable == "flowering.vital.rate.4" ~ "Ambient",
                           variable == "flowering.vital.rate.5" ~ "Elevated",
                           variable == "flowering.vital.rate.6" ~ "Elevated",
                           variable == "flowering.vital.rate.7" ~ "Ambient",
                           variable == "flowering.vital.rate.8" ~ "Ambient"),
         
         endo = case_when(variable == "survival.vital.rate.1" ~ "Negative",
                          variable == "survival.vital.rate.2" ~ "Positive",
                          variable == "survival.vital.rate.3" ~ "Negative",
                          variable == "survival.vital.rate.4" ~ "Positive",
                          variable == "survival.vital.rate.5" ~ "Negative",
                          variable == "survival.vital.rate.6" ~ "Positive",
                          variable == "survival.vital.rate.7" ~ "Negative",
                          variable == "survival.vital.rate.8" ~ "Positive",
                          
                          variable == "flowering.vital.rate.1" ~ "Negative",
                          variable == "flowering.vital.rate.2" ~ "Positive",
                          variable == "flowering.vital.rate.3" ~ "Negative",
                          variable == "flowering.vital.rate.4" ~ "Positive",
                          variable == "flowering.vital.rate.5" ~ "Negative",
                          variable == "flowering.vital.rate.6" ~ "Positive",
                          variable == "flowering.vital.rate.7" ~ "Negative",
                          variable == "flowering.vital.rate.8" ~ "Positive"),
         
         year = case_when(variable == "survival.vital.rate.1" ~ "14-15",
                          variable == "survival.vital.rate.2" ~ "14-15",
                          variable == "survival.vital.rate.3" ~ "14-15",
                          variable == "survival.vital.rate.4" ~ "14-15",
                          variable == "survival.vital.rate.5" ~ "15-16",
                          variable == "survival.vital.rate.6" ~ "15-16",
                          variable == "survival.vital.rate.7" ~ "15-16",
                          variable == "survival.vital.rate.8" ~ "15-16",
                          
                          variable == "flowering.vital.rate.1" ~ "14-15",
                          variable == "flowering.vital.rate.2" ~ "14-15",
                          variable == "flowering.vital.rate.3" ~ "14-15",
                          variable == "flowering.vital.rate.4" ~ "14-15",
                          variable == "flowering.vital.rate.5" ~ "15-16",
                          variable == "flowering.vital.rate.6" ~ "15-16",
                          variable == "flowering.vital.rate.7" ~ "15-16",
                          variable == "flowering.vital.rate.8" ~ "15-16")) %>% 
  mutate(year = as.factor(year),
         water = as.factor(water)) %>% 
unite_("water_endo", c("water", "endo"), sep="_", remove=F) %>% 
  unite_("water_yr", c("water", "year"), sep= "_", remove = F) %>% 
  mutate(water_endo_factor = as.factor(water_endo)) %>% 
  bind_rows(p_rec_post)
          
          ggplot(posteriors)+
  geom_density(aes(x=posterior,fill=endo),alpha = .85)+
 # scale_fill_manual(values=c("#ffae19","#dd7509", "#68BAE9", "#0338a8"))+
#  scale_color_manual(values=c("#ffae19","#dd7509", "#68BAE9", "#0338a8"))+
  facet_grid(vital_rate~water_yr, scales = "free")+
  theme_classic()
          

negative <- posteriors %>% 
            filter(endo == "Negative") %>% 
            select(variable, posterior, vital_rate, water, year, endo) %>% 
  dplyr::rename(endo_neg = endo,
                   posterior_neg = posterior)
                   
positive <- posteriors %>% 
            filter(endo == "Positive") %>% 
            select(posterior, endo) %>% 
  dplyr::rename(endo_pos = endo,
                   posterior_pos = posterior) 
          
        
delta_post <- negative %>% 
              bind_cols(positive) %>% 
              mutate(posterior_delta = posterior_pos - posterior_neg) %>% 
              unite("water_yr", c("water", "year"), remove = F)



ggplot(delta_post)+
  geom_density(aes(x=posterior_delta, y=..scaled..,fill=water),alpha = .90)+
 # scale_fill_manual(values =c("gray39", "black"))+
  
  scale_fill_manual(values=c("#ffae19", "#68BAE9", "#0338a8"))+
  scale_color_manual(values=c("#ffae19", "#68BAE9", "#0338a8"))+
  
 # scale_fill_manual(values=c("#ffae19","#dd7509", "#68BAE9", "#0338a8"))+
#  scale_color_manual(values=c("#ffae19","#dd7509", "#68BAE9", "#0338a8"))+
  facet_grid(vital_rate~year, scales = "free_y")+
  theme_classic(base_size = 18)+
  geom_vline(xintercept = 0, linetype = "dashed")
  




 post_flow <- posteriors %>% 
  filter(vital_rate == "flowering",
         water == "Ambient",
         year == "14-15")
          
    

other_list <- c("prob_ep_rec_new.ambient.14",
"prob_em_rec_new.ambient.14",

"prob_ep_rec_new.irrigated.14",
"prob_em_rec_new.irrigated.14",

"prob_ep_rec_new.ambient.15",
"prob_em_rec_new.ambient.15",

"prob_ep_rec_new.irrigated.15",
"prob_em_rec_new.irrigated.15",

# "ep_seeds_ambient.14",
# "ep_seeds_ambient.15",
# 
# "em_seeds_ambient.14",
# "em_seeds_ambient.15",
# 
# "ep_seeds_irrigated.14",
# "ep_seeds_irrigated.15",
# 
# "em_seeds_irrigated.14",
# "em_seeds_irrigated.15",


"^vtrans")


posteriors_other <- as.data.frame(AGHY.endochange.20142016$BUGSoutput$sims.list) %>% 
  gather("variable", "posterior") %>% 
  filter(str_detect(variable, other_list))








plot_data <- AGHY_plot_dat %>% 
  dplyr::rename(plot = newplot) %>% 
  mutate(water = case_when(water == 1 ~ "Ambient",
                           water == 2 ~ "Irrigated"))

bayes.eps <-  bayes.data.summary.surv %>% 
  filter(str_detect(info, "^eps")) %>% 
  separate(info, into = c("eps", "type", "year", "plot")) %>% 
  select(eps, type, year, plot, mean) %>% 
  dplyr::rename(mean_eps = mean) %>% 
  mutate(year = case_when(year == "1" ~ 2014,
                          year == "2" ~ 2015),
         plot = as.numeric(plot)) %>% 
  left_join(plot_data) %>% 
  filter(year != "NA")

save.image("AGHY_model_Bayes_full_output_update3.RData")


bayes.flow_frequency <- bayes.flow %>% 
  select(mean, endo, water, year) %>% 
  mutate(year = as.numeric(year)) %>% 
  mutate(vital_rate = "Flowering",
         type = "flow")


bayes.eps.flow <- bayes.eps %>% 
  filter(type == "flow")

bayes.flow_freq <- bayes.eps.flow %>% 
  left_join(bayes.flow_frequency, by = c("year", "water", "type")) %>% 
  left_join(prev_dat, by = c("plot", "year", "water")) %>% 
  unite(water_endo, c("endo", "water"), sep = "_", remove = F) %>% 
  mutate(plot_means = mean_eps + mean)

freq_plot_flow <- ggplot(bayes.flow_freq, aes(mean.prev, plot_means, color = water, linetype = water))+
  geom_point(aes(color = water), alpha =.75)+
  geom_smooth(aes(),se=F, size = .6)+
 scale_color_manual(values = c("orange3", "steelblue4"))+
  theme_classic()+
  facet_grid(endo~year)+
  labs(y = "Plot-specific means",
       x = "Population prevalence")


bayes.surv_frequency <- bayes.surv %>% 
  select(mean, endo, water, year) %>% 
  mutate(year = as.numeric(year)) %>% 
  mutate(vital_rate = "Survival",
         type ="surv")

bayes.eps.surv <- bayes.eps %>% 
  filter(type == "surv") %>% 
  mutate(year = case_when(year == 2014 ~ 2015,
                          year == 2015 ~ 2016))

bayes.surv_freq <- bayes.eps.surv %>% 
  left_join(bayes.surv_frequency, by = c("year", "water", "type")) %>% 
  left_join(prev_dat, by = c("plot", "year", "water")) %>% 
  mutate(plot_means = mean_eps + mean)

freq_plot_surv <- ggplot(bayes.surv_freq, aes(mean.prev, plot_means, linetype = water))+
  geom_point(aes(color = water))+
  geom_smooth(aes(color = water),se=F, size = .75)+
  scale_color_manual(values = c("orange3", "steelblue4"))+
  scale_linetype_manual(values = c("dashed", "solid"))+
  theme_classic()+
  facet_grid(endo~year)+
  labs(y = "Plot-specific epsilons",
       x = "Population prevalence")

bayes.surv_freq_pos <- bayes.surv_freq %>% 
  filter(endo == "Positive")

freq_plot_surv <- ggplot(bayes.surv_freq_pos, aes(mean.prev, plot_means, linetype = water))+
  geom_point(aes(color = water))+
  geom_smooth(aes(color = water),se=F, size = .75)+
  scale_color_manual(values = c("orange3", "steelblue4"))+
  scale_linetype_manual(values = c("dashed", "solid"))+
  theme_classic()+
  facet_grid(endo~year)+
  labs(y = "Plot-specific epsilons",
       x = "Population prevalence")

bayes.surv_freq_neg <- bayes.surv_freq %>% 
  filter(endo == "Negative")

freq_plot_surv_neg <- ggplot(bayes.surv_freq_neg, aes((mean.prev), plot_means, linetype = water))+
  geom_point(aes(color = water))+
  geom_smooth(aes(color = water),se=F, size = .75)+
  scale_color_manual(values = c("orange3", "steelblue4"))+
  scale_linetype_manual(values = c("dashed", "solid"))+
  theme_classic()+
  facet_grid(endo~year)+
  labs(y = "Plot-specific epsilons",
       x = "Population prevalence")

bayes.rec_frequency <- bayes.rec %>% 
  separate(info, c("blah", "endo", "type", "bdfsadf", "water", "year")) %>% 
  select(endo, water, mean, water, year) %>% 
  mutate(year = as.numeric(year),
         #endo = "Positive",
         water = case_when(water == "ambient" ~ "Ambient",
                           water == "irrigated" ~ "Irrigated")) %>% 
  mutate(vital_rate = "Recruitment",
         type = "rec") %>% 
  mutate(year = case_when(year == 14 ~ 2015,
                          year == 15 ~ 2016),
endo = case_when(endo == "em" ~ "Negative",
                 endo == "ep" ~ "Positive"))
  
  bayes.eps.rec <- bayes.eps %>% 
  filter(type == "rec") %>% 
  mutate(year = case_when(year == 2014 ~ 2015,
                          year == 2015 ~ 2016))

bayes.rec_freq <- bayes.eps.rec %>% 
  left_join(bayes.rec_frequency, by = c("year", "water", "type")) %>% 
  left_join(prev_dat, by = c("plot", "year", "water")) %>% 
  mutate(vital_rate = "Recruitment") %>% 
  mutate(plot_means = mean_eps + mean)

freq_plot_rec <- ggplot(bayes.rec_freq, aes(mean.prev, plot_means))+
  geom_point(aes(color = water))+
  geom_smooth(aes(color = water),se=F, linetype = "dashed", size = .75)+
  scale_color_manual(values = c("orange3", "steelblue4"))+
  theme_classic()+
  facet_grid(endo~year)+
  labs(y = "Plot-specific epsilons",
       x = "Population prevalence")


bayes_rec_surv_freq <- bayes.rec_freq %>% 
  bind_rows(bayes.surv_freq) %>% 
  mutate(trans_year = case_when(year == 2015 ~ "'14-'15",
                                year == 2016 ~ "'15-'16"))

freq_rec_surv_plot <- ggplot(bayes_rec_surv_freq, aes(mean.prev, mean_eps))+
  geom_point(aes(color = water, shape = endo), alpha = .6)+
  #geom_smooth(aes(color = water, linetype = endo),se=F, size = .6)+
  geom_smooth(aes(color = water, linetype = endo),se=F, size = .6, method = "lm")+
  scale_color_manual(values = c("orange3", "steelblue4"))+
  theme_classic()+
  facet_grid(vital_rate~trans_year, scales = "free_y")+
  labs(y = "Plot-specific epsilons",
       x = "Population prevalence")+
  theme(legend.position = "bottom")


bayes_seed_freq <- bayes.flow.percap.mu %>% 
  mutate(type = type_mean) %>% 
  mutate(year = as.numeric(year))

bayes.eps.seed_infs_freq <- bayes.eps %>% 
  filter(type == "seed" |
           type == "infs") %>% 
  mutate(year = as.numeric(year),
         type = ifelse(type == "seed", "seeds", type))


bayes.seeds_infs_flow_freq <- bayes.eps.seed_infs_freq %>% 
  left_join(bayes_seed_freq, by = c("year", "water", "type")) %>% 
  left_join(prev_dat, by = c("plot", "year", "water")) %>% 
  bind_rows(bayes.flow_freq) %>% 
  mutate(vital_rate = case_when(type == "infs" ~ "Inflorescences",
                                type == "seed" ~ "Seeds",
                                type == "flow" ~ "Flowering"))

bayes_infs_freq <- bayes.seeds_infs_flow_freq %>% 
  filter(type == "infs" |
           type == "seeds") %>% 
  mutate(mean = means_real,
  plot_means = means_real + mean_eps)



bayes_seeds_freq  <- bayes.seeds_infs_flow_freq %>% 
  left_join(bayes_seed_freq) %>% 
  left_join(prev_dat) %>% 
  bind_rows(bayes.flow_freq) %>% 
  mutate(plot_means = mean + mean_eps) %>% 
  filter(type != "infs") %>% 
  filter(type != "seeds")

new_df <- bayes_seeds_freq %>% 
  bind_rows(bayes_infs_freq) %>% 
  mutate(vital_rate = case_when(type == "infs" ~ "Inflorescences",
                                type == "seeds" ~ "Seeds",
                                type == "flow" ~ "Flowering"),
         trans_year = case_when(year == 2014 ~ "'14-'15",
                              year == 2015 ~ "'15-'16")) 


freq_seeds_infs_flow_plot <- ggplot(new_df, aes(mean.prev, plot_means))+
  geom_jitter(aes(color = water, shape = endo), alpha = 1)+
 # geom_smooth(aes(color = water),se=F, linetype = "dashed", size = .6)+
  geom_smooth(aes(color = water, linetype = endo),se=F, size = .6, method = "lm")+
  scale_color_manual(values = c("orange3", "steelblue4"))+
  theme_classic()+
  facet_grid(vital_rate~trans_year, scales = "free_y")+
  labs(y = "Plot-specific epsilons",
       x = "Population prevalence")+
  theme(legend.position = "bottom")


freq_all_df <- new_df %>% 
  bind_rows(bayes_rec_surv_freq)


### THIS IS THE FINAL FIGURE FOR FREQUENCY DEPENDENCE 
freq_all_plot <- ggplot(freq_all_df, aes(mean.prev, plot_means))+
  geom_jitter(aes(color = water, shape = endo), alpha = 1)+
  # geom_smooth(aes(color = water),se=F, linetype = "dashed", size = .6)+
  geom_smooth(aes(color = water, linetype = endo),se=F, size = .6, method = "lm")+
  scale_color_manual(values = c("orange3", "steelblue4"))+
  theme_classic()+
  facet_grid(vital_rate~trans_year, scales = "free_y")+
  labs(y = "Plot-specific epsilons",
       x = "Population prevalence")+
  theme(legend.position = "bottom")



library(cowplot)
plot_grid(freq_seeds_infs_flow_plot, freq_rec_surv_plot)




new_df_pos <- bayes_seeds_freq %>% 
  bind_rows(bayes_infs_freq) %>% 
  mutate(vital_rate = case_when(type == "infs" ~ "Inflorescences",
                                type == "seeds" ~ "Seeds",
                                type == "flow" ~ "Flowering")) %>% 
  filter(endo == "Positive")


freq_seeds_infs_flow_plot_pos <- ggplot(new_df_pos, aes(mean.prev, plot_means, color = water, shape = endo))+
  geom_jitter(aes(color = water, shape = endo), alpha = .7)+
  geom_smooth(aes(color = water),se=F, linetype = "dashed", size = .6)+
  scale_color_manual(values = c("orange3", "steelblue4"))+
  theme_classic()+
  facet_grid(vital_rate~year, scales = "free_y")+
  labs(y = "Endophyte positive Plot-specific epsilons",
       x = "Population prevalence")+
  theme(legend.position = "bottom")


new_df_neg <- bayes_seeds_freq %>% 
  bind_rows(bayes_infs_freq) %>% 
  mutate(vital_rate = case_when(type == "infs" ~ "Inflorescences",
                                type == "seeds" ~ "Seeds",
                                type == "flow" ~ "Flowering")) %>% 
  filter(endo == "Negative")


freq_seeds_infs_flow_plot_neg <- ggplot(new_df_neg, aes(mean.prev, plot_means, color = water, shape = endo))+
  geom_jitter(aes(color = water, shape = endo), alpha = .7)+
  geom_smooth(aes(color = water),se=F, linetype = "dashed", size = .6)+
  scale_color_manual(values = c("orange3", "steelblue4"))+
  theme_classic()+
  facet_grid(vital_rate~year, scales = "free_y")+
  labs(y = "Endophyte negative Plot-specific epsilons",
       x = "Population prevalence")+
  theme(legend.position = "bottom")
  

plot_grid(freq_seeds_infs_flow_plot_pos, freq_seeds_infs_flow_plot_neg)

