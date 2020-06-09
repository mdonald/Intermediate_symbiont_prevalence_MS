### Intermediate symbiont prevalence manuscript
### JAGS parameterization from R
## pulls the data that are bundled in the "Intermediate_Prevalence_Data_Bundle_for_JAGS_Script.R"


load("Data_for_JAGS_Model.RData")
library(R2jags);load.module("mix")
library(mcmcplots)


## list of data for JAGS model
jag.data<-list(N.trt=N.trt,
               N.years=N.years,
               N.plots=N.plots,
               N.obs=N.obs,
               water=water,
               y.pos=y.pos, ## population prevalence data
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
               endo.k = survival_known$endo,
               water.k = survival_known$water,
               year.k = survival_known$year,
               plot.k = survival_known$plot,
               
               N.obs.surv.unknown = nrow(survival_unk_sbplt),
               alive.unk.t1 = survival_unk_sbplt$alive.t1,
               alive.unk.t = survival_unk_sbplt$alive.t,
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
               
               ### FERTILITY DATA
               N.obs.seed = nrow(seedmass_dat),
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
               
               ## ADJUST POP Prev by infs of known endo status

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
               
               ### flowering data for known endo status
               N.obs.flow.known = nrow(flowering_individuals),
               endo.flow.k = flowering_individuals$endo,
               water.flow.k = flowering_individuals$water,
               year.flow.k = flowering_individuals$year,
               flowering.k = flowering_individuals$flow_t,
               
               ### flowering data for unknown endo status
               N.plots.known.endo.flow = nrow(flowering_plts),
               plot.adj.flow = flowering_plts$plot,
               year.adj.flow = flowering_plts$year,
               total.plants.3.sbplts.f = flowering_plts$total_plants_f,
               known_ep_f = flowering_plts$total_known_ep_f,
               total_unknown_flow = flowering_plts$total_unknown_f,
               
               ### flowering data for weighted mean approach
               
               N.obs.flow.unknown = nrow(flowering_sbplts),
               year.f = flowering_sbplts$year_t,
               plot.f = flowering_sbplts$newplot,
               water.f = flowering_sbplts$water,
               total.unk.plants.f = flowering_sbplts$total_unknown,
               flow.unk = flowering_sbplts$total_unknown_flowered,
               
               
               ## linear regression seed mass and seed count data
               N.obs.seed.lin.reg = nrow(linreg_dat),
               seed.count.linreg = (linreg_dat$seed_count),
               seed.mass.linreg = (linreg_dat$seed_mass),
               
               
               ### recruitment 
               N.sbplts.avg.rec = nrow(rec_final),
               rec.sbplt = rec_final$total_plants_in_sbplt,
               water.r = rec_final$water,
               plot.r = rec_final$plot,
               year.r = rec_final$year,
               
               plots.rec_t = rec_plot_final_yr$plot,
               water.rec_t = rec_plot_final_yr$water,
               year.rec_t = rec_plot_final_yr$year,
               
               surv.rec.plot = survival_dat2$newplot,
               N.plots.surv.rec = nrow(survival_dat2),
               water.surv.rec = survival_dat2$water,
               year.surv.rec = survival_dat2$year,
               plots.rec_t1 = rec_plot_final_yr1$plot,
               year.rec_t1 = rec_plot_final_yr1$year,

               
               ## vertical transmission
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
                       q=runif(1,0,10),
                       
                       ## survival inits
                       mu_surv = array(rnorm(jag.data$N.endo*jag.data$N.water*jag.data$N.year),
                                       dim=c(jag.data$N.endo,jag.data$N.water,jag.data$N.year)),
                       sigma.surv = rlnorm(1),
                       ## FLOWERING AND RECRUITMENT
                       
                       mu.seed = array(rnorm(jag.data$N.endo*jag.data$N.water*jag.data$N.year),
                                       dim=c(jag.data$N.endo,jag.data$N.water,jag.data$N.year)),
                       mu.infs = array(rnorm(jag.data$N.endo*jag.data$N.water*jag.data$N.year),
                                       dim=c(jag.data$N.endo,jag.data$N.water,jag.data$N.year)),
                       sigma.res = rlnorm(1),
                       sigma.infs.overdisp = rlnorm(1),
                       sigma.plot = rlnorm(1),
                       sigma.plot.infs = rlnorm(1),
                       mu_flow = array(rnorm(jag.data$N.endo*jag.data$N.water*jag.data$N.year),
                                       dim=c(jag.data$N.endo,jag.data$N.water,jag.data$N.year)),
                                             sigma.flow = rlnorm(1),
                       w=runif(1,0,10),
                       
                       beta0_seed_linreg = rnorm(1,0,1),
                       slope_lin_reg = rnorm(1,0,1),
                       
                       sigma.count = rlnorm(1),
                       
                       ## recruitment
                       mu_rec = array(rnorm(jag.data$N.water*jag.data$N.years),
                                      dim=c(jag.data$N.water,jag.data$N.years)),
                                           sigma.rec = rlnorm(1),
                       ## VERTICAL TRANSMISSION
                       mu_vtrans = array(rnorm(jag.data$N.water*jag.data$N.year),
                                         dim=c(jag.data$N.water,jag.data$N.year)),
                       sigma.vt = rlnorm(1)
)
}

## Params to estimate
parameters<-c(
              # population prevalence
              "Eplus.add.pred", ## elevated precipitation fitted relationship
              "Eplus.control.pred", ## ambient precipitation fitted relationship
              "prev", ## population prevalence estimate
              
              # vital rates
              "survival.vital.rate", ## probabilty of survival vital rate
              
              "flowering.vital.rate", ## probability of flowering vital rate
              
              "seedmass.per.cap", ## per capita reproduction vital rate
        
              "prob_ep_rec_new.ambient.14", ## derived probability of E+ recruitment - ambient 2014-2015
              "prob_ep_rec_new.ambient.15", ## derived probabilty of E+ recruitment - ambient 2015-2016
              
              "prob_ep_rec_new.irrigated.14", ## derived probability of E+ recruitment - elevated 2014-2015
              "prob_ep_rec_new.irrigated.15", ## derived probability of E+ recruitment - elevated 2015-2016
              
              "prob_em_rec_new.ambient.14", ## derived probability of E- recruitment - ambient 2014-2015
              "prob_em_rec_new.ambient.15", ## derived probability of E- recruitment - ambient 2015-2016
              
              "prob_em_rec_new.irrigated.14", ## derived probability of E- recruitment - elevated 2014-2015
              "prob_em_rec_new.irrigated.15", ## derived probability of E- recruitment - elevated 2015-2016
              
              "vtrans", ## probability of vertial transmission
              
              ### random effects to check fequency dependence for survival, flowering, recruitment, seed mass, inflorescence production, and vertical transmission
              "eps.surv",
              "eps.flow",
              "eps.rec",
              "eps.seed",
              "eps.infs",
              "eps.vt",
              
              ## posterior predictive checks
              
              # population prevalence
              "fit",
              "fit.new",
              
              # survival
              "fit.surv", 
              "fit.surv.new",
              
              # seed mass
              "fit.seed",
              "fit.seed.new",
              
              # inflorescence production 
              "fit.infs",
              "fit.infs.new",
              
              # flowering
              "fit.flow",
              "fit.flow.new",
              
              # seed mass (linear regression)
              "fit.linreg",
              "fit.linreg.new",
              
              # recruitment
              "fit.rec",
              "fit.rec.new",
              
              # vertical transmission
              "fit.vtrans",
              "fit.vtrans.new"
              
)


## MCMC settings
ni<-10000
nb<-2000
nt<-20
nc<-3

## run JAGS
results_bayes<-jags(data=jag.data,inits=inits,parameters.to.save=parameters,
           model.file="./vital_rate_scripts/Bayes models/Intermediate_Prevalence_JAGS_script.txt",
           n.thin=nt,n.chains=nc,n.burnin=nb,
           n.iter=ni,working.directory=getwd())

#save.image("Intermediate_Prev_JAGS_output_May30.RData") ## data object for figures
