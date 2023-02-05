#Install package:
#install.packages("EpiModel", dependencies = TRUE)

#INITIAL METRICS CA PH 
#.	R0: 2-3 (unmitigated) or 0.9-1.5 (mitigated)
#.	Infection fatality ratio: 0.5-1.5%
#.	Hospitalization rate: 5-20%
#.	Doubling time before social distancing: 6-7 days
#.	Current hospitalized patients: county to determine
#.	Social distancing - % reduction: 30%-70%
#.	ICU percent: 1-5%
#.	Ventilated percent: 1-5%
#.	Hospital (non-ICU) length of stay: 7 days
#.	ICU length of stay (in addition to non-ICU length of stay): 9 days
#.	Length of time on a ventilator if intubated: 10 days
#.	Hospital market share: county to determine
#.	Regional population: county to determine
#.	Currently known regional infections: county to determine
#Our metrics for the SBD model for homeless

#ASSUMPTIONS 
#.	5760 estimate population
#.	>211 days: We have no intervention for SBD county homeless.
#.	Ro=Rt : 2.20
#.	5% Case fatality rate for homeless increases by a factor of 5 because of their co-morbidities and stressful living situation. It could be up to 50% due to this publication: 
#  .	23.6% Hospitalization rates: 9/38 (23.6%);ICU rates 5/38 (13.2%); deaths 1/38 (2.6%)
#.	15 days: Average length of stay: No current data on SB, 10-20 day hospitalization average 17. 
#.	8 days: Duration patient is infectious, I brought this up to 8 from the 2.9 default value to reflect the homeless lifestyle. 



##################################################LIBRARY2###################################################################################################################
#install.packages("tidyverse")
library(dplyr)
library(tidyverse)
library(magrittr)
library(lubridate)
library(stringr)
library(tibble)
library(broom)
library(ggplot2)
library(gt)
#install.packages("gt")
library(knitr)
#install.packages("devtools")
library(devtools)
#install.packages("DiagrammeR")
library(DiagrammeR)
library(parallel)
library(foreach)
#install.packages("tictoc")
library(tictoc)
#install.packages("network")
library(network)
suppressMessages(library(EpiModel))
#install.packages("incidence")
library(incidence)
#install.packages("earlyR")
library(earlyR)





tic("Time to complete")

#GITHUBGISTS
source_files <- c("_icm.mod.init.seiqhrf.R", "_icm.mod.status.seiqhrf.R", 
                  "_icm.mod.vital.seiqhrf.R", "_icm.control.seiqhrf.R", "_icm.utils.seiqhrf.R", 
                  "_icm.saveout.seiqhrf.R", "_icm.icm.seiqhrf.R")

src_path <- paste0("./_posts/2020-03-18-modelling-the-effects-of-public-health-", 
                   "interventions-on-covid-19-transmission-part-2/")

gist_url <- "https://gist.github.com/timchurches/92073d0ea75cfbd387f91f7c6e624bd7"

local_source <- FALSE

for (source_file in source_files) {
  if (local_source) {
    source(paste(src_path, source_file, sep = ""))
  } else {
    source_gist(gist_url, filename = source_file)
  }
}


# function to set-up and run the baseline simulations
simulate <- function(# control.icm params),
  type = "SEIQHRF", 
  nsteps = 366, 
  nsims = 8,
  ncores = 4,
  prog.rand = FALSE,
  rec.rand = FALSE,
  fat.rand = TRUE,
  quar.rand = FALSE,
  hosp.rand = FALSE,
  disch.rand = TRUE,
  infection.FUN = infection.seiqhrf.icm,
  recovery.FUN = progress.seiqhrf.icm,
  departures.FUN = departures.seiqhrf.icm,
  arrivals.FUN = arrivals.icm,
  get_prev.FUN = get_prev.seiqhrf.icm,
  # init.icm params
  s.num = 1878,
  e.num=1,
  i.num = 6,
  q.num=1,
  h.num=1,
  r.num = 1,
  f.num = 1,
  # param.icm params
  inf.prob.e = 0.02,#1/2171603*100=0.0004 
  act.rate.e = 10,
  inf.prob.i = 0.05,#3/2171603*100= 0.00013
  act.rate.i = 10,
  inf.prob.q = 0.02,#1/2171603= 0.000046 
  act.rate.q = 2.5,                    
  quar.rate = 1/30, 
  hosp.rate = 5/100,#12/100,15/100,16/100
  disch.rate = 1/15,
  prog.rate = 1/10,
  prog.dist.scale = 5,
  prog.dist.shape = 1.5,
  rec.rate = 1/20,
  rec.dist.scale = 35,
  rec.dist.shape = 1.5,
  fat.rate.base = 1/50,
  hosp.cap = 12118,# 95
  fat.rate.overcap = 1/25,
  fat.tcoeff = 0.5,
  vital = TRUE,
  a.rate = (10.5/365)/1000, 
  a.prop.e = 0.01,
  a.prop.i = 0.001,
  a.prop.q = 0.01,
  ds.rate = (7/365)/1000, 
  de.rate = (7/365)/1000, 
  di.rate = (7/365)/1000,
  dq.rate = (7/365)/1000,
  dh.rate = (20/365)/1000,
  dr.rate = (7/365)/1000,
  out="mean"
) 
{
  control <- control.icm(type = type, 
                         nsteps = nsteps, 
                         nsims = nsims,
                         ncores = ncores,
                         prog.rand = prog.rand,
                         rec.rand = rec.rand,
                         infection.FUN = infection.FUN,
                         recovery.FUN = recovery.FUN,
                         arrivals.FUN = arrivals.FUN,
                         departures.FUN = departures.FUN,
                         get_prev.FUN = get_prev.FUN)
  
  init <- init.icm(s.num = s.num,
                   e.num = e.num,
                   i.num = i.num,
                   q.num = q.num,
                   h.num = h.num,
                   r.num = r.num,
                   f.num = f.num)
  
  param <-  param.icm(inf.prob.e = inf.prob.e, 
                      act.rate.e = act.rate.e,
                      inf.prob.i = inf.prob.i, 
                      act.rate.i = act.rate.i,
                      inf.prob.q = inf.prob.q, 
                      act.rate.q = act.rate.q,                    
                      quar.rate = quar.rate,
                      hosp.rate = hosp.rate,
                      disch.rate = disch.rate,
                      prog.rate = prog.rate,
                      prog.dist.scale = prog.dist.scale,
                      prog.dist.shape = prog.dist.shape,
                      rec.rate = rec.rate,
                      rec.dist.scale = rec.dist.scale,
                      rec.dist.shape = rec.dist.shape,
                      fat.rate.base = fat.rate.base,
                      hosp.cap = hosp.cap,
                      fat.rate.overcap = fat.rate.overcap,
                      fat.tcoeff = fat.tcoeff,
                      vital = vital,
                      a.rate = a.rate, 
                      a.prop.e = a.prop.e,
                      a.prop.i = a.prop.i,
                      a.prop.q = a.prop.q,
                      ds.rate = ds.rate, 
                      de.rate = de.rate, 
                      di.rate = di.rate,
                      dq.rate = dq.rate,
                      dh.rate = dh.rate,
                      dr.rate = dr.rate)
  
  sim <- icm.seiqhrf(param, init, control)
  sim_df <- as.data.frame(sim, out=out)
  
  return(list(sim=sim, df=sim_df))
}

icm.seiqhrf <- function(param, init, control) {
  
  crosscheck.icm(param, init, control)
  verbose.icm(control, type = "startup")
  nsims <- control$nsims
  ncores <- ifelse(control$nsims == 1, 1, min(future::availableCores(), control$ncores))
  control$ncores <- ncores
  
  if (ncores == 1) {
    
    # Simulation loop start
    for (s in 1:control$nsims) {
      
      ## Initialization module
      if (!is.null(control[["initialize.FUN"]])) {
        dat <- do.call(control[["initialize.FUN"]], list(param, init, control))
      }
      
      
      # Timestep loop
      for (at in 2:control$nsteps) {
        
        ## User Modules
        um <- control$user.mods
        if (length(um) > 0) {
          for (i in 1:length(um)) {
            dat <- do.call(control[[um[i]]], list(dat, at))
          }
        }
        
        ## Infection
        if (!is.null(control[["infection.FUN"]])) {
          dat <- do.call(control[["infection.FUN"]], list(dat, at))
        }
        
        
        ## Recovery
        if (!is.null(control[["recovery.FUN"]])) {
          dat <- do.call(control[["recovery.FUN"]], list(dat, at))
        }
        
        
        ## Departure Module
        if (!is.null(control[["departures.FUN"]])) {
          dat <- do.call(control[["departures.FUN"]], list(dat, at))
        }
        
        
        ## Arrival Module
        if (!is.null(control[["arrivals.FUN"]])) {
          dat <- do.call(control[["arrivals.FUN"]], list(dat, at))
        }
        
        
        ## Outputs
        if (!is.null(control[["get_prev.FUN"]])) {
          dat <- do.call(control[["get_prev.FUN"]], list(dat, at))
        }
        
        
        ## Track progress
        verbose.icm(dat, type = "progress", s, at)
      }
      
      # Set output
      if (s == 1) {
        out <- saveout.seiqhrf.icm(dat, s)
      } else {
        out <- saveout.seiqhrf.icm(dat, s, out)
      }
      
    } # Simulation loop end
    
    class(out) <- "icm"
    
  } # end of single core execution
  
  if (ncores > 1) {  
    doParallel::registerDoParallel(ncores)
    sout <- foreach(s = 1:nsims, .export = c('infection.seiqhrf.icm','arrivals.seiqhrf.icm',
                                             'control.icm','cum_discr_si','departures.seiqhrf.icm',
                                             'get_prev.seiqhrf.icm','init_status.icm','initialize.icm',
                                             'saveout.seiqhrf.icm','infection.seiqhrf.icm','progress.seiqhrf.icm'),
                    
                    .packages = 'EpiModel') %dopar% {    
    # sout <- foreach(s = 1:nsims) %dopar% {
      
      control$nsims <- 1
      control$currsim <- s
      
      ## Initialization module
      if (!is.null(control[["initialize.FUN"]])) {
        dat <- do.call(control[["initialize.FUN"]], list(param, init, control))
      }
      
      # Timestep loop
      for (at in 2:control$nsteps) {
        
        ## User Modules
        um <- control$user.mods
        if (length(um) > 0) {
          for (i in 1:length(um)) {
            dat <- do.call(control[[um[i]]], list(dat, at))
          }
        }
        
        ## Infection
        if (!is.null(control[["infection.FUN"]])) {
          dat <- do.call(control[["infection.FUN"]], list(dat, at))
        }
        
        
        ## Recovery
        if (!is.null(control[["recovery.FUN"]])) {
          dat <- do.call(control[["recovery.FUN"]], list(dat, at))
        }
        
        
        ## Departure Module
        if (!is.null(control[["departures.FUN"]])) {
          dat <- do.call(control[["departures.FUN"]], list(dat, at))
        }
        
        
        ## Arrival Module
        if (!is.null(control[["arrivals.FUN"]])) {
          dat <- do.call(control[["arrivals.FUN"]], list(dat, at))
        }
        
        
        ## Outputs
        if (!is.null(control[["get_prev.FUN"]])) {
          dat <- do.call(control[["get_prev.FUN"]], list(dat, at))
        }
        
        
        ## Track progress
        verbose.icm(dat, type = "progress", s, at)
      }
      
      # Set output
      out <- saveout.seiqhrf.icm(dat, s=1)
      class(out) <- "icm"
      return(out)
      
    }
    
    # aggregate results collected from each thread
    collected_times <- list()
    
    # collect the times from sout then delete them
    for (i in 1:length(sout)) {
      collected_times[[paste0("sim", i)]] <- sout[[i]]$times$sim1 
      sout[[i]]$times <- NULL
    }
    
    # merge $epi structures
    merged.out <- sout[[1]]
    for (i in 2:length(sout)) {
      merged.out <- merge.seiqhrf.icm(merged.out, sout[[i]], param.error = FALSE)
    }
    out <- merged.out
    
    # add the collected timing data
    out$times <- collected_times
    
    class(out) <- "icm"
  } # end of parallel execution
  
  return(out)
}


#BASELINE SIMULATION
baseline_sim <- simulate(ncores = 4)









#######################################define a function to extract timings and assemble a data#############################################################################
# frame
get_times <- function(simulate_results) {
  
  sim <- simulate_results$sim
  
  for (s in 1:sim$control$nsims) {
    if (s == 1) {
      times <- sim$times[[paste("sim", s, sep = "")]]
      times <- times %>% mutate(s = s)
    } else {
      times <- times %>% bind_rows(sim$times[[paste("sim", 
                                                    s, sep = "")]] %>% mutate(s = s))
    }
  }
  
  times <- times %>% mutate(infTime = ifelse(infTime < 0, -5, 
                                             infTime), expTime = ifelse(expTime < 0, -5, expTime)) %>% 
    mutate(incubation_period = infTime - expTime, illness_duration = recovTime - 
             expTime, illness_duration_hosp = dischTime - expTime, 
           hosp_los = dischTime - hospTime, quarantine_delay = quarTime - 
             infTime, survival_time = fatTime - infTime) %>% 
    select(s, incubation_period, quarantine_delay, illness_duration, 
           illness_duration_hosp, hosp_los, survival_time) %>% 
    pivot_longer(-s, names_to = "period_type", values_to = "duration") %>% 
    mutate(period_type = factor(period_type, levels = c("incubation_period", 
                                                        "quarantine_delay", "illness_duration", "illness_duration_hosp", 
                                                        "hosp_los", "survival_time"), labels = c("Incubation period", 
                                                                                                 "Delay entering isolation", "Illness duration", "Illness duration (hosp)", 
                                                                                                 "Hospital care required duration", "Survival time of case fatalities"), 
                                ordered = TRUE))
  return(times)
}
#now get timing data
times <- get_times(baseline_sim)
#visualize it 
times %>% filter(duration <= 30) %>% ggplot(aes(x = duration)) + 
  geom_bar() + facet_grid(period_type ~ ., scales = "free_y") + 
  labs(title = "Duration frequency distributions", subtitle = "Baseline simulation")


#You could argue that these are day-only admissions, but really, hospital care should be required for at least one day




##############################################################VISUAL PREVALENCE#############################################################################################

#BASELINE PLOT

baseline_plot_df <- baseline_sim$df %>% # use only the prevalence columns
  select(time, s.num, e.num, i.num, q.num, h.num, r.num, f.num) %>% 
  # examine only the first 100 days since it is all over by
  # then using the default parameters
  filter(time <= 100) %>% pivot_longer(-c(time), names_to = "compartment", 
                                       values_to = "count")


# define a standard set of colours to represent compartments######
compcols <- c(s.num = "yellow", e.num = "orange", i.num = "red", 
              q.num = "cyan", h.num = "magenta", r.num = "lightgreen", 
              f.num = "black")
complabels <- c(s.num = "Susceptible", e.num = "Infected/asymptomatic", 
                i.num = "Infected/infectious", q.num = "Self-isolated", h.num = "Requires hospitalisation", 
                r.num = "Recovered", f.num = "Case fatality")

baseline_plot_df %>% ggplot(aes(x = time, y = count, colour = compartment)) + 
  geom_line(size = 2, alpha = 0.7) + scale_colour_manual(values = compcols, 
                                                         labels = complabels) + theme_dark() + labs(title = "Baseline simulation", 
       
                                                                                                    
                                                                                                    
#Notice that the whole population being infected                                                                                                    
                                                                                                                                                                                             x = "Days since beginning of epidemic", y = "Prevalence (persons)")

#leave out compartments############
baseline_plot_df %>% filter(compartment %in% c("e.num", "i.num", 
                                               "q.num", "h.num", "f.num")) %>% ggplot(aes(x = time, y = count, 
                                                                                          colour = compartment)) + geom_line(size = 2, alpha = 0.7) + 
  scale_colour_manual(values = compcols, labels = complabels) + 
  theme_dark() + labs(title = "Baseline simulation", x = "Days since beginning of epidemic", 
                      y = "Prevalence (persons)")

view(baseline_plot_df) #look at compartments 

#What can we see in this plot of our simulation of our hypothetical world of number of people 


#Checking the basic reproduction number of Ro######
# get the S-> E compartment flow, which is our daily
# incidence rate
incidence_counts <- baseline_sim$df %>% select(time, se.flow)
# uncount them
incident_case_dates <- incidence_counts %>% uncount(se.flow) %>% 
  pull(time)
# convert to an incidence object
incidence_all <- incident_case_dates %>% incidence(.)

# plot the incidence curve
plot(incidence_all)



View(baseline_plot_df)#Views table of number of infected over time for each run (corresponding to 3 CFRs)



###########################################################find the peak of the epidemic curve###############################################################################
peak_of_epidemic_curve <- find_peak(incidence_all)

# repeat with just the growth part of the epidemic curve
incident_case_dates_growth_phase <- incidence_counts %>% filter(time <= 
                                                                  peak_of_epidemic_curve) %>% select(time, se.flow) %>% uncount(se.flow) %>% 
  pull(time)

incidence_growth_phase <- incident_case_dates_growth_phase %>% 
  incidence(., last_date = peak_of_epidemic_curve)
# specify serial interval mean and SD since the last blog
# post new studies have appeared suggesting 4.5 is a better
# mean for the SI
si_mean <- 4.5
si_sd <- 3.4

# get a set of MLE estimates for R0 and plot
res <- get_R(incidence_growth_phase, si_mean = si_mean, si_sd = si_sd)
plot(res, "R")



###########################################Running an intervention experiment################################################################################################
quar_rate_ramp <- function(t) {
  ifelse(t < 15, 0.0333, ifelse(t <= 30, 0.0333 + (t - 15) * 
                                  (0.5 - 0.0333)/15, 0.5))
}



ramp_quar_rate_sim <- simulate(quar.rate = quar_rate_ramp(1:366))

#Let's model the effect of increasing the rate at which people enter (self-)isolation
#ecause we can use time-variant parameters, we will! Thus we ramp up the isolation rate (quar.rate) 
#from it's low level of 0.033 (as in the baseline simulation)
#starting at day 15, up to 0.5 at day 30. We'll write a little function to do that.



ramp_quar_rate_sim_plot_df <- ramp_quar_rate_sim$df %>% # DEFINING NEW DF with SIMS
  select(time, s.num, e.num, i.num, q.num, h.num, r.num, f.num) %>% 
  # examine only the first 100 days since it is all over by
  # then using the default parameters
  filter(time <= 100) %>% pivot_longer(-c(time), names_to = "compartment", 
                                       values_to = "count")

# define a standard set of colours to represent compartments
compcols <- c(s.num = "yellow", e.num = "orange", i.num = "red", 
              q.num = "cyan", h.num = "magenta", r.num = "lightgreen", 
              f.num = "black")
complabels <- c(s.num = "Susceptible", e.num = "Infected/asymptomatic", 
                i.num = "Infected/infectious", q.num = "Self-isolated", h.num = "Requires hospitalisation", 
                r.num = "Recovered", f.num = "Case fatality")

baseline_plot_df %>% ggplot(aes(x = time, y = count, colour = compartment)) + 
  geom_line(size = 2, alpha = 0.7) + scale_colour_manual(values = compcols, 
                                                         labels = complabels) + theme_dark() + labs(title = "Baseline simulation", 
                                                                                                    x = "Days since beginning of epidemic", y = "Prevalence (persons)")





baseline_plot_df %>% mutate(experiment = "Baseline") %>% bind_rows(ramp_quar_rate_sim_plot_df %>% 
                                                                     mutate(experiment = "Ramp up isolation")) %>% filter(compartment %in% 
                                                                                                                            c("e.num", "i.num", "q.num", "h.num", "f.num")) %>% ggplot(aes(x = time, 
                                                                                                                                                                                           y = count, colour = compartment)) + geom_line(size = 2, alpha = 0.7) + 
  facet_grid(experiment ~ .) + scale_colour_manual(values = compcols, 
                                                   labels = complabels) + theme_dark() + labs(title = "Baseline vs ramping up isolation simulations", 
                                                                                              x = "Days since beginning of epidemic", y = "Prevalence (persons)")


#####################################################experiment 2 hospital beds#############################################################################################
#Over a four week period, let's ramp up hospital capacity to triple the baseline level, 
#starting at day 15. Hey, China built a 1000+ bed COVID-19 hospital in Wuhan in just 10 days.

hosp_cap_ramp <- function(t) {
  ifelse(t < 15, 40, ifelse(t <= 45, 40 + (t - 15) * (120 - 
                                                        40)/30, 120))
}

raise_hosp_cap_sim <- simulate(hosp.cap = hosp_cap_ramp(1:366))

raise_hosp_cap_df <- raise_hosp_cap_sim$df %>% # DEFINING NEW DF with SIMS
  select(time, s.num, e.num, i.num, q.num, h.num, r.num, f.num) %>% 
  # examine only the first 100 days since it is all over by
  # then using the default parameters
  filter(time <= 100) %>% pivot_longer(-c(time), names_to = "compartment", 
                                       values_to = "count")

# define a standard set of colours to represent compartments
compcols <- c(s.num = "yellow", e.num = "orange", i.num = "red", 
              q.num = "cyan", h.num = "magenta", r.num = "lightgreen", 
              f.num = "black")
complabels <- c(s.num = "Susceptible", e.num = "Infected/asymptomatic", 
                i.num = "Infected/infectious", q.num = "Self-isolated", h.num = "Requires hospitalisation", 
                r.num = "Recovered", f.num = "Case fatality")

baseline_plot_df %>% ggplot(aes(x = time, y = count, colour = compartment)) + 
  geom_line(size = 2, alpha = 0.7) + scale_colour_manual(values = compcols, 
                                                         labels = complabels) + theme_dark() + labs(title = "Baseline simulation", 
                                                                                                    x = "Days since beginning of epidemic", y = "Prevalence (persons)")

baseline_plot_df %>% ggplot(aes(x = time, y = count, colour = compartment)) + 
  geom_line(size = 2, alpha = 0.7) + scale_colour_manual(values = compcols, 
                                                         labels = complabels) + theme_dark() + labs(title = "Baseline simulation", 
                                                                                                    x = "Days since beginning of epidemic", y = "Prevalence (persons)")





baseline_plot_df %>% mutate(experiment = "Baseline") %>% bind_rows(ramp_quar_rate_sim_plot_df %>% 
                                                                     mutate(experiment = "More Hospital Beds")) %>% filter(compartment %in% 
                                                                                                                            c("e.num", "i.num", "q.num", "h.num", "f.num")) %>% ggplot(aes(x = time, 
                                                                                                                                                                                           y = count, colour = compartment)) + geom_line(size = 2, alpha = 0.7) + 
  facet_grid(experiment ~ .) + scale_colour_manual(values = compcols, 
                                                   labels = complabels) + theme_dark() + labs(title = "Baseline vs More Hospital Beds", 
                                                                                              x = "Days since beginning of epidemic", y = "Prevalence (persons)")



######################################################################EXPERIMENT 4 - More social distancing at day 15#######################################################
social_distancing_day15_ramp <- function(t) {
  ifelse(t < 15, 10, ifelse(t <= 30, 10 - (t - 15) * (10 - 
                                                        5)/15, 5))
}

t15_social_distancing_sim <- simulate(act.rate.i = social_distancing_day15_ramp(1:366), 
                                      act.rate.e = social_distancing_day15_ramp(1:366))


t15_social_distancing_df <- t15_social_distancing_sim$df %>% # DEFINING NEW DF with SIMS
  select(time, s.num, e.num, i.num, q.num, h.num, r.num, f.num) %>% 
  # examine only the first 100 days since it is all over by
  # then using the default parameters
  filter(time <= 100) %>% pivot_longer(-c(time), names_to = "compartment", 
                                       values_to = "count")

# define a standard set of colours to represent compartments
compcols <- c(s.num = "yellow", e.num = "orange", i.num = "red", 
              q.num = "cyan", h.num = "magenta", r.num = "lightgreen", 
              f.num = "black")
complabels <- c(s.num = "Susceptible", e.num = "Infected/asymptomatic", 
                i.num = "Infected/infectious", q.num = "Self-isolated", h.num = "Requires hospitalisation", 
                r.num = "Recovered", f.num = "Case fatality")

baseline_plot_df %>% ggplot(aes(x = time, y = count, colour = compartment)) + 
  geom_line(size = 2, alpha = 0.7) + scale_colour_manual(values = compcols, 
                                                         labels = complabels) + theme_dark() + labs(title = "Baseline simulation", 
                                                                                                    x = "Days since beginning of epidemic", y = "Prevalence (persons)")

baseline_plot_df %>% ggplot(aes(x = time, y = count, colour = compartment)) + 
  geom_line(size = 2, alpha = 0.7) + scale_colour_manual(values = compcols, 
                                                         labels = complabels) + theme_dark() + labs(title = "Baseline simulation", 
                                                                                                    x = "Days since beginning of epidemic", y = "Prevalence (persons)")





baseline_plot_df %>% mutate(experiment = "Baseline") %>% bind_rows(t15_social_distancing_df %>% 
                                                                     mutate(experiment = "Social Distancing Day 15")) %>% filter(compartment %in% 
                                                                                                                             c("e.num", "i.num", "q.num", "h.num", "f.num")) %>% ggplot(aes(x = time, 
                                                                                                                                                                                            y = count, colour = compartment)) + geom_line(size = 2, alpha = 0.7) + 
  facet_grid(experiment ~ .) + scale_colour_manual(values = compcols, 
                                                   labels = complabels) + theme_dark() + labs(title = "Baseline vs Social Distancing at Day 15", 
                                                                                              x = "Days since beginning of epidemic", y = "Prevalence (persons)")


#############################################################EXPERIMENT 4 - More social distancing at day 30################################################################
social_distancing_day30_ramp <- function(t) {
  ifelse(t < 30, 10, ifelse(t <= 45, 10 - (t - 30) * (10 - 
                                                        5)/15, 5))
}

t30_social_distancing_sim <- simulate(act.rate.i = social_distancing_day30_ramp(1:366), 
                                      act.rate.e = social_distancing_day30_ramp(1:366))



t30_social_distancing_df <- t30_social_distancing_sim$df %>% # DEFINING NEW DF with SIMS
  select(time, s.num, e.num, i.num, q.num, h.num, r.num, f.num) %>% 
  # examine only the first 100 days since it is all over by
  # then using the default parameters
  filter(time <= 100) %>% pivot_longer(-c(time), names_to = "compartment", 
                                       values_to = "count")

# define a standard set of colours to represent compartments
compcols <- c(s.num = "yellow", e.num = "orange", i.num = "red", 
              q.num = "cyan", h.num = "magenta", r.num = "lightgreen", 
              f.num = "black")
complabels <- c(s.num = "Susceptible", e.num = "Infected/asymptomatic", 
                i.num = "Infected/infectious", q.num = "Self-isolated", h.num = "Requires hospitalisation", 
                r.num = "Recovered", f.num = "Case fatality")

baseline_plot_df %>% ggplot(aes(x = time, y = count, colour = compartment)) + 
  geom_line(size = 2, alpha = 0.7) + scale_colour_manual(values = compcols, 
                                                         labels = complabels) + theme_dark() + labs(title = "Baseline simulation", 
                                                                                                    x = "Days since beginning of epidemic", y = "Prevalence (persons)")

baseline_plot_df %>% ggplot(aes(x = time, y = count, colour = compartment)) + 
  geom_line(size = 2, alpha = 0.7) + scale_colour_manual(values = compcols, 
                                                         labels = complabels) + theme_dark() + labs(title = "Baseline simulation", 
                                                                                                    x = "Days since beginning of epidemic", y = "Prevalence (persons)")





baseline_plot_df %>% mutate(experiment = "Baseline") %>% bind_rows(t30_social_distancing_df %>% 
                                                                     mutate(experiment = "Social Distancing Day 30")) %>% filter(compartment %in% 
                                                                                                                                   c("e.num", "i.num", "q.num", "h.num", "f.num")) %>% ggplot(aes(x = time, 
                                                                                                                                                                                                  y = count, colour = compartment)) + geom_line(size = 2, alpha = 0.7) + 
  facet_grid(experiment ~ .) + scale_colour_manual(values = compcols, 
                                                   labels = complabels) + theme_dark() + labs(title = "Baseline vs Social Distancing at Day 30", 
                                                                                              x = "Days since beginning of epidemic", y = "Prevalence (persons)")



##################################Experiment 5 - increase both social distancing and increased self-isolation rates starting day 15#########################################
quar_rate_ramp <- function(t) {
  ifelse(t < 15, 0.0333, ifelse(t <= 30, 0.0333 + (t - 15) * 
                                  (0.5 - 0.0333)/15, 0.5))
}

ramp_quar_rate_sim <- simulate(quar.rate = quar_rate_ramp(1:366))
t15_soc_dist_quar_sim <- simulate(act.rate.i = social_distancing_day15_ramp(1:366), 
                                  act.rate.e = social_distancing_day15_ramp(1:366), quar.rate = quar_rate_ramp(1:366))


t15_soc_dist_quar_df <- t15_soc_dist_quar_sim$df %>% # DEFINING NEW DF with SIMS
  select(time, s.num, e.num, i.num, q.num, h.num, r.num, f.num) %>% 
  # examine only the first 100 days since it is all over by
  # then using the default parameters
  filter(time <= 100) %>% pivot_longer(-c(time), names_to = "compartment", 
                                       values_to = "count")

# define a standard set of colours to represent compartments
compcols <- c(s.num = "yellow", e.num = "orange", i.num = "red", 
              q.num = "cyan", h.num = "magenta", r.num = "lightgreen", 
              f.num = "black")
complabels <- c(s.num = "Susceptible", e.num = "Infected/asymptomatic", 
                i.num = "Infected/infectious", q.num = "Self-isolated", h.num = "Requires hospitalisation", 
                r.num = "Recovered", f.num = "Case fatality")

baseline_plot_df %>% ggplot(aes(x = time, y = count, colour = compartment)) + 
  geom_line(size = 2, alpha = 0.7) + scale_colour_manual(values = compcols, 
                                                         labels = complabels) + theme_dark() + labs(title = "Baseline simulation", 
                                                                                                    x = "Days since beginning of epidemic", y = "Prevalence (persons)")

baseline_plot_df %>% ggplot(aes(x = time, y = count, colour = compartment)) + 
  geom_line(size = 2, alpha = 0.7) + scale_colour_manual(values = compcols, 
                                                         labels = complabels) + theme_dark() + labs(title = "Baseline simulation", 
                                                                                                    x = "Days since beginning of epidemic", y = "Prevalence (persons)")





baseline_plot_df %>% mutate(experiment = "Baseline") %>% bind_rows(t15_soc_dist_quar_df %>% 
                                                                     mutate(experiment = "Social Distancing Day 15 amd Increased Self-Isolation")) %>% filter(compartment %in% 
                                                                                                                                   c("e.num", "i.num", "q.num", "h.num", "f.num")) %>% ggplot(aes(x = time, 
                                                                                                                                                                                                  y = count, colour = compartment)) + geom_line(size = 2, alpha = 0.7) + 
  facet_grid(experiment ~ .) + scale_colour_manual(values = compcols, 
                                                   labels = complabels) + theme_dark() + labs(title = "Baseline vs Social Distancing and Self-Quarantine", 
                                                                                              x = "Days since beginning of epidemic", y = "Prevalence (persons)")


