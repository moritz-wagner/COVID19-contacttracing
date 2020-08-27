require(patchwork)
require(ringbp)
require(tidyverse)
require(data.table)

source('R/contact_functions.R')
source('R/outbreak_functions.R')

## Create synthetic population based on contact data-----------
freq = c("Daily","Often","Regularly","Rarely","Never") #Can adjust this depending on which types of contacts to include
p.urban = .3
t <- 7
pop <- 230*1000

# Age-specific relative susceptibility estimates from https://www.nature.com/articles/s41591-020-0962-9
rel.susc.age = c(.4,.4,(.4+.38)/2,.38,(.79+.86+.8)/3,(.82+.88+.74)/3) 
rel.susc.age <- rel.susc.age/max(rel.susc.age)

CD <- count.contacts4(contact_data,t=t,p.urban = p.urban,freq=freq,pop=pop,rel.susc.age=rel.susc.age)
CD[[1]] %>% nrow #Population size
study_participants <- CD[[1]]
two_week_contacts_subset <- CD[[3]]
setnames(CD[[1]],old = "caseid",new = "csid")
CD[[4]] <- two_week_contacts_subset[,.(nhh.contact_ids=list(grep("nhh",unlist(contact_id),value=TRUE)),
                                       nhh_weight=list(unlist(weight)[grep("nhh",unlist(contact_id))])),by=age]
CD[[4]] <- CD[[4]][,`:=`(N=length(unlist(nhh_weight))),by=age]
CD[[5]] <- study_participants[,.(csid=list(csid)),by=age]


# Setup infection and intervention paramaters---------------
sims <- 100

# Setup - When do interventions start? (not used in paper)
cap_cases <- Inf # Sim stops with this many cumulative cases
cap_weeks <- 8 # Sim stops after this many weeks
cap_max_days <- (cap_weeks-1)*7
cap_contacts_per_day <- Inf # Stop intervention if daily contacts (HH or non-HH) to trace go above this
cap_max_hhcon <- cap_contacts_per_day*cap_max_days
cap_max_nhhcon <- cap_contacts_per_day*cap_max_days
intervention.trigger <- 1 #number of cumulative cases when intervention should be triggered

# Setup infection parameters
num.initial.cases <- 5
r0.lev <- "medium" #R0 - low, medium, high
disp <- "medium" #Overdispersion - low, medium, high
iso_delay <- "medium" #Delay to isolation - short, medium, long
pre_symp <- "medium" #Delay to isolation - low, medium, high
rel.inf.asym <- 1 #Relative infectiousness of asymptomatics (compared to sympatomatics)
# scale.symp <- FALSE

r0community <- ifelse(r0.lev=="high",3,ifelse(r0.lev=="medium",2.5,2))
disp.com <- ifelse(disp=="high",0.16,ifelse(disp=="medium",0.58,2))
r0isolated = 0;disp.iso = 1
k <- ifelse(pre_symp=="high",0.2,ifelse(pre_symp=="medium",0.7,1.9))
if (iso_delay=="medium") {
  delay_shape = 1.651524;delay_scale = 4.287786
} else if (iso_delay=="short") {
  delay_shape = 1.5;delay_scale = 1.5
} else {
  delay_shape = 2.5;delay_scale = 5 
}
# Age-specific clinical fraction estimates from https://www.nature.com/articles/s41591-020-0962-9
prop.sym.age = c(0.29,0.29,(.21+.28)/2,.21,(.27+.33+.4)/3,(.49+.63+.69)/3)
prop.asym <- 1-prop.sym.age

##roughly scale up r0 if asymptomatics are less infectious
if (rel.inf.asym!=1) {
  # if (scale.symp) {
  #   r0community <- r0community/.75
  # } else {
    r0.scale <- sum(CD[[1]][,.(con=sum(hh.con+nhh.con)),by=age][,`:=`(con=con/sum(con))][c(4,2,3,1,5,6),2]*prop.asym)
    r0community <- r0community/r0.scale
  # }
}


## Run sims for different interventions-------------
# No Isolation (assumes prop asymtomatic = 1)-----
quarantine <- FALSE
HH_trace <- FALSE
prop.ascertain <- 0
age.limits <- levels(contact_data$age_class_part)
intervention.age <- age.limits
intervention.age <- NULL
prop.non.HH <- 1
lim.non.HH <- Inf

time.ext.par<-system.time(res.ext <- scenario_sim_ext_parallel(n.sim = sims,num.initial.cases = num.initial.cases, 
                                                               prop.ascertain = prop.ascertain, cap_max_days = cap_max_days, 
                                                               cap_cases = cap_cases, cap_max_hhcon=cap_max_hhcon, cap_max_nhhcon=cap_max_nhhcon, 
                                                               r0isolated = r0isolated, r0community = r0community, 
                                                               rel.inf.asym=rel.inf.asym,
                                                               disp.iso = disp.iso, disp.com = disp.com, 
                                                               delay_shape = Inf, 
                                                               delay_scale = delay_scale, 
                                                               k = k, prop.asym = prop.asym, 
                                                               quarantine = quarantine,HH_trace = HH_trace,
                                                               contact_data=CD,p.urban = p.urban,
                                                               t=t,freq = freq,
                                                               pop=pop,
                                                               intervention.trigger=intervention.trigger,
                                                               intervention.age=intervention.age,
                                                               prop.non.HH=prop.non.HH,lim.non.HH=lim.non.HH))

time.ext.par/60

res.no.int <- res.ext$res.week
res.no.int$intervention <- "No Interventions"

res.no.int.gen <- res.ext$res.gen
res.no.int.gen$intervention <- "No Interventions"

# Isolation only symptomatic (i.e. 20% symp)-----
quarantine <- FALSE
HH_trace <- FALSE
prop.ascertain <- 0
age.limits <- levels(contact_data$age_class_part)
intervention.age <- age.limits
intervention.age <- NULL
prop.non.HH <- 1
lim.non.HH <- Inf

time.ext.par<-system.time(res.ext <- scenario_sim_ext_parallel(n.sim = sims,num.initial.cases = num.initial.cases, 
                                                               prop.ascertain = prop.ascertain, cap_max_days = cap_max_days, 
                                                               cap_cases = cap_cases, cap_max_hhcon=cap_max_hhcon, cap_max_nhhcon=cap_max_nhhcon, 
                                                               r0isolated = r0isolated, r0community = r0community, 
                                                               rel.inf.asym=rel.inf.asym,
                                                               disp.iso = disp.iso, disp.com = disp.com, delay_shape = delay_shape, 
                                                               delay_scale = delay_scale, k = k, prop.asym = prop.asym, 
                                                               quarantine = quarantine,HH_trace = HH_trace,
                                                               contact_data=CD,p.urban = p.urban,
                                                               t=t,freq = freq,
                                                               pop=pop,
                                                               intervention.trigger=intervention.trigger,
                                                               intervention.age=intervention.age,
                                                               prop.non.HH=prop.non.HH,lim.non.HH=lim.non.HH))

time.ext.par/60

res.no.tracing2 <- res.ext$res.week
res.no.tracing2$intervention <- "Isolation (all symptomatic)"

res.no.tracing.gen2 <- res.ext$res.gen
res.no.tracing.gen2$intervention <- "Isolation (all symptomatic)"

# Isolation only symptomatic and no SSEs-----
quarantine <- FALSE
HH_trace <- FALSE
prop.ascertain <- 0
age.limits <- levels(contact_data$age_class_part)
intervention.age <- age.limits
prop.non.HH <- 1
lim.non.HH <- Inf

time.ext.par<-system.time(res.ext <- scenario_sim_ext_parallel(n.sim = sims,num.initial.cases = num.initial.cases, 
                                                               prop.ascertain = prop.ascertain, cap_max_days = cap_max_days, 
                                                               cap_cases = cap_cases, cap_max_hhcon=cap_max_hhcon, cap_max_nhhcon=cap_max_nhhcon, 
                                                               r0isolated = r0isolated, r0community = r0community, 
                                                               rel.inf.asym=rel.inf.asym,
                                                               disp.iso = disp.iso, disp.com = disp.com, delay_shape = delay_shape, 
                                                               delay_scale = delay_scale, k = k, prop.asym = prop.asym, 
                                                               quarantine = quarantine,HH_trace = HH_trace,
                                                               contact_data=CD,p.urban = p.urban,
                                                               t=t,freq = freq,
                                                               pop=pop,
                                                               intervention.trigger=intervention.trigger,
                                                               intervention.age=intervention.age,
                                                               prop.non.HH=prop.non.HH,lim.non.HH=lim.non.HH))

time.ext.par/60

res.iso.SSE <- res.ext$res.week
res.iso.SSE$intervention <- "Isolation (all symptomatic) + no SSEs"

res.iso.SSE.gen2 <- res.ext$res.gen
res.iso.SSE.gen2$intervention <- "Isolation (all symptomatic) + no SSEs"

# Tracing HH + 25% non-HH-----
quarantine <- FALSE
HH_trace <- TRUE
prop.ascertain <- 0.25
age.limits <- levels(contact_data$age_class_part)
intervention.age <- age.limits
intervention.age <- NULL
prop.non.HH <- 1
lim.non.HH <- Inf

time.ext.par<-system.time(res.ext <- scenario_sim_ext_parallel(n.sim = sims,num.initial.cases = num.initial.cases, 
                                                               prop.ascertain = prop.ascertain, cap_max_days = cap_max_days, 
                                                               cap_cases = cap_cases, cap_max_hhcon=cap_max_hhcon, cap_max_nhhcon=cap_max_nhhcon, 
                                                               r0isolated = r0isolated, r0community = r0community, 
                                                               rel.inf.asym=rel.inf.asym,
                                                               disp.iso = disp.iso, disp.com = disp.com, delay_shape = delay_shape, 
                                                               delay_scale = delay_scale, k = k, prop.asym = prop.asym, 
                                                               quarantine = quarantine,HH_trace = HH_trace,
                                                               contact_data=CD,p.urban = p.urban,
                                                               t=t,freq = freq,
                                                               pop=pop,
                                                               intervention.trigger=intervention.trigger,
                                                               intervention.age=intervention.age,
                                                               prop.non.HH=prop.non.HH,lim.non.HH=lim.non.HH))

time.ext.par/60

res.tracing.HH2 <- res.ext$res.week
res.tracing.HH2$intervention <- "Tracing HH + 25% non-HH"

res.tracing.HH2.gen <- res.ext$res.gen
res.tracing.HH2.gen$intervention <- "Tracing HH + 25% non-HH"

# Quaranting HH + 25% non-HH-----
quarantine <- TRUE
HH_trace <- TRUE
prop.ascertain <- 0.25
age.limits <- levels(contact_data$age_class_part)
intervention.age <- age.limits
intervention.age <- NULL
prop.non.HH <- 1
lim.non.HH <- Inf

time.ext.par<-system.time(res.ext <- scenario_sim_ext_parallel(n.sim = sims,num.initial.cases = num.initial.cases, 
                                                               prop.ascertain = prop.ascertain, cap_max_days = cap_max_days, 
                                                               cap_cases = cap_cases, cap_max_hhcon=cap_max_hhcon, cap_max_nhhcon=cap_max_nhhcon, 
                                                               r0isolated = r0isolated, r0community = r0community, 
                                                               rel.inf.asym=rel.inf.asym,
                                                               disp.iso = disp.iso, disp.com = disp.com, delay_shape = delay_shape, 
                                                               delay_scale = delay_scale, k = k, prop.asym = prop.asym, 
                                                               quarantine = quarantine,HH_trace = HH_trace,
                                                               contact_data=CD,p.urban = p.urban,
                                                               t=t,freq = freq,
                                                               pop=pop,
                                                               intervention.trigger=intervention.trigger,
                                                               intervention.age=intervention.age,
                                                               prop.non.HH=prop.non.HH,lim.non.HH=lim.non.HH))

time.ext.par/60

res.quarantine2 <- res.ext$res.week
res.quarantine2$intervention <- "Quarantine HH + 25% non-HH"

res.quarantine2.gen <- res.ext$res.gen
res.quarantine2.gen$intervention <-  "Quarantine HH + 25% non-HH"


# General SD, reduction in non-HH contacts to max 5 (distinct in 7 days)-----
quarantine <- FALSE
HH_trace <- FALSE
age.limits <- levels(contact_data$age_class_part)
intervention.age <- age.limits
prop.non.HH <- 1
lim.non.HH <- 5
prop.ascertain <- 0

time.ext.par<-system.time(res.ext <- scenario_sim_ext_parallel(n.sim = sims,num.initial.cases = num.initial.cases, 
                                                               prop.ascertain = prop.ascertain, cap_max_days = cap_max_days, 
                                                               cap_cases = cap_cases, cap_max_hhcon=cap_max_hhcon, cap_max_nhhcon=cap_max_nhhcon, 
                                                               r0isolated = r0isolated, r0community = r0community, 
                                                               rel.inf.asym=rel.inf.asym,
                                                               disp.iso = disp.iso, disp.com = disp.com, delay_shape = delay_shape, 
                                                               delay_scale = delay_scale, k = k, prop.asym = prop.asym,
                                                               quarantine = quarantine,HH_trace = HH_trace,
                                                               contact_data=CD,p.urban = p.urban,
                                                               t=t,freq = freq,
                                                               pop=pop,
                                                               intervention.trigger=intervention.trigger,
                                                               intervention.age=intervention.age,
                                                               prop.non.HH=prop.non.HH,lim.non.HH=lim.non.HH))

time.ext.par/60

res.sd2 <- res.ext$res.week
res.sd2$intervention <- "Social Distancing (Max 5 non-HH contacts)"

res.sd2.gen <- res.ext$res.gen
res.sd2.gen$intervention <-  "Social Distancing (Max 5 non-HH contacts)"

# Shielding of elderly, reduction in non-HH contacts to none-----
quarantine <- FALSE
HH_trace <- FALSE
age.limits <- levels(contact_data$age_class_part)
intervention.age <- age.limits[6]
prop.non.HH <- 1
lim.non.HH <- 0
prop.ascertain <- 0

time.ext.par<-system.time(res.ext <- scenario_sim_ext_parallel(n.sim = sims,num.initial.cases = num.initial.cases, 
                                                               prop.ascertain = prop.ascertain, cap_max_days = cap_max_days, 
                                                               cap_cases = cap_cases, cap_max_hhcon=cap_max_hhcon, cap_max_nhhcon=cap_max_nhhcon, 
                                                               r0isolated = r0isolated, r0community = r0community, 
                                                               rel.inf.asym=rel.inf.asym,
                                                               disp.iso = disp.iso, disp.com = disp.com, delay_shape = delay_shape, 
                                                               delay_scale = delay_scale, k = k, prop.asym = prop.asym,
                                                               quarantine = quarantine,HH_trace = HH_trace,
                                                               contact_data=CD,p.urban = p.urban,
                                                               t=t,freq = freq,
                                                               pop=pop,
                                                               intervention.trigger=intervention.trigger,
                                                               intervention.age=intervention.age,
                                                               prop.non.HH=prop.non.HH,lim.non.HH=lim.non.HH))

time.ext.par/60

res.shielding2 <- res.ext$res.week
res.shielding2$intervention <- "Shielding elderly (0 non-HH contacts for over 50 year olds)"

res.shielding2.gen <- res.ext$res.gen
res.shielding2.gen$intervention <-  "Shielding elderly (0 non-HH contacts for over 50 year olds)"

# School closures, reduction in non-HH contacts of 6-19 year olds to none-----
quarantine <- FALSE
HH_trace <- FALSE
age.limits <- levels(contact_data$age_class_part)
intervention.age <- age.limits[c(3,4)]
prop.non.HH <- 1
lim.non.HH <- 0
prop.ascertain <- 0

time.ext.par<-system.time(res.ext <- scenario_sim_ext_parallel(n.sim = sims,num.initial.cases = num.initial.cases, 
                                                               prop.ascertain = prop.ascertain, cap_max_days = cap_max_days, 
                                                               cap_cases = cap_cases, cap_max_hhcon=cap_max_hhcon, cap_max_nhhcon=cap_max_nhhcon, 
                                                               r0isolated = r0isolated, r0community = r0community, 
                                                               rel.inf.asym=rel.inf.asym,
                                                               disp.iso = disp.iso, disp.com = disp.com, delay_shape = delay_shape, 
                                                               delay_scale = delay_scale, k = k, prop.asym = prop.asym,
                                                               quarantine = quarantine,HH_trace = HH_trace,
                                                               contact_data=CD,p.urban = p.urban,
                                                               t=t,freq = freq,
                                                               pop=pop,
                                                               intervention.trigger=intervention.trigger,
                                                               intervention.age=intervention.age,
                                                               prop.non.HH=prop.non.HH,lim.non.HH=lim.non.HH))

time.ext.par/60

res.shool <- res.ext$res.week
res.shool$intervention <- "School closures (0 non-HH contacts for 6-19 year olds)"

res.shool.gen <- res.ext$res.gen
res.shool.gen$intervention <-  "School closures (0 non-HH contacts for 6-19 year olds)"

##Put them together----------

res <- bind_rows(list(res.no.int,res.no.tracing2,res.iso.SSE,
                      res.shielding2,res.shool,res.sd2,
                      res.tracing.HH2,
                      res.quarantine2))
res.gen <- bind_rows(list(res.no.int.gen,res.no.tracing.gen2,res.iso.SSE.gen2,
                          res.shielding2.gen,res.shool.gen,res.sd2.gen,
                          res.tracing.HH2.gen,
                          res.quarantine2.gen))

res %<>% mutate(week=week+1)

res$Intervention <- factor(res$intervention,levels = unique(res$intervention))
res.gen$Intervention <- factor(res.gen$intervention,levels = unique(res.gen$intervention))

save(res,res.gen,CD,file = paste0('R/runs/agesymp_r0',r0.lev,'_disp',disp,'_relasym',rel.inf.asym*100,'_delay_',iso_delay,'_pre_',pre_symp,'_init5_',sims,'sims_',pop/1000,'pop.RData'))

# # SD + Quaranting HH + # of non-HH-----
# # sims <- 100
quarantine <- TRUE
HH_trace <- FALSE
prop.ascertain <- seq(0,1,length.out = 5)
age.limits <- levels(contact_data$age_class_part)
intervention.age <- age.limits
lim.non.HH <- c(0,1,5,10,Inf)
prop.non.HH <- 1
res.ext <- list()
tic <- Sys.time()
for (j in seq_along(lim.non.HH)) {
  for (i in seq_along(prop.ascertain)) {
    ind <- i+(j-1)*length(lim.non.HH)
    time.ext.par<-system.time(res.ext[[ind]] <- scenario_sim_ext_parallel(n.sim = sims,num.initial.cases = num.initial.cases,
                                                                          prop.ascertain = prop.ascertain[i], cap_max_days = cap_max_days,
                                                                          cap_cases = cap_cases, cap_max_hhcon=cap_max_hhcon, cap_max_nhhcon=cap_max_nhhcon,
                                                                          r0isolated = r0isolated, r0community = r0community,
                                                                          rel.inf.asym=rel.inf.asym,
                                                                          disp.iso = disp.iso, disp.com = disp.com, delay_shape = delay_shape,
                                                                          delay_scale = delay_scale, k = k, prop.asym = prop.asym,
                                                                          quarantine = quarantine,HH_trace = HH_trace,
                                                                          contact_data=CD,p.urban = p.urban,
                                                                          t=t,freq = freq,
                                                                          pop=pop,
                                                                          intervention.trigger=intervention.trigger,
                                                                          intervention.age=intervention.age,
                                                                          prop.non.HH=prop.non.HH,lim.non.HH=lim.non.HH[j]))
    
    res.ext[[ind]]$res.week$prop.ascertain <- prop.ascertain[i]
    res.ext[[ind]]$res.week$prop.non.HH <- prop.non.HH
    res.ext[[ind]]$res.week$lim.non.HH <- lim.non.HH[j]
    
    res.ext[[ind]]$res.gen$prop.ascertain <- prop.ascertain[i]
    res.ext[[ind]]$res.gen$prop.non.HH <- prop.non.HH
    res.ext[[ind]]$res.gen$lim.non.HH <- lim.non.HH[j]
    
    print(ind)
    print(prop.ascertain[i])
    print(lim.non.HH[j])
    print(time.ext.par)
  }
}
toc <- Sys.time()
toc-tic

save(res.ext,file = paste0('R/runs/agesymp_quarantineSD2_r0',r0.lev,'_disp',disp,'_relasym',rel.inf.asym*100,'_delay_',iso_delay,'_pre_',pre_symp,'_init5_',sims,'sims_',pop/1000,'.RData'))

