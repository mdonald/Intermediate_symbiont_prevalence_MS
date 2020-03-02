## Title: AGHY Full Bayes model, prevalence and vital rates figures
## Purpose: Use the full distribution from pop prev estimates to estimate vital rates
## Notes: Adjusted pop prev works for survival, flowering, fertility and recruitment! 
## Date Started: 2015
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
               
               
              
             
             ##  rec
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

              
               "prob.ep.surv", 
               "prob.em.surv",
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

              "p.rec",
              "ep_rec_plt_t1",
              "ep_rec_plt_t",
              "ep_rec_flow_t",
              "ep_seeds_plt_t",
              
              "ep_seeds_rec",
              "ep_rec_plt",

             "mu_vtrans",
              "mean.vtrans",

              

               "new_prob_ep_rec",
             
               "prob_ep_rec_new.ambient.14",
               "prob_ep_rec_new.ambient.15",
              
             "prob_ep_rec_new.irrigated.14",
              "prob_ep_rec_new.irrigated.15",

              "prob_em_rec_new.ambient.14",
              "prob_em_rec_new.ambient.15",
              
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
