# First working attempt at extending the LSHTM model to allow for depleting susceptibles
# as informed by contact data
# The popolation being modelled is closed and has household structure
# ==============================================================================
# --- Unique contacts over a 2-week period
source("R/new_contacts.R")
t <- 14 
participants=unique(cd$csid)
two_week_contacts=count.contacts(participants[1],cd,t) 
for(p in participants[2:length(participants)]) two_week_contacts=rbind(two_week_contacts,count.contacts(p,cd,t))
two_week_contacts=two_week_contacts[,c(9,1:8)]
# --- some differences between daily data and two_week data
sum(two_week_contacts$how_often=='Daily' | two_week_contacts$how_often=='Never')
sum(cd$how_often=='Daily' | cd$how_often=='Never')

# --- Restricting contacts of interest to 'Daily' and 'Never' contacts for now 
# so as to work with a smaller matrix
two_week_contacts_subset=two_week_contacts[(two_week_contacts$how_often=='Daily' | two_week_contacts$how_often=='Never'),]
# Weight for the types of contact
c.type=c(1, 0.000001, 1.5/7, 0.5/30, 1.5/30, 2/365)# Daily, Missing, Often, Rarely, Regularly, Never
two_week_contacts_subset$weight=c.type[as.numeric(two_week_contacts_subset$how_often)]
# ID's for study participants
study_participants=unique(two_week_contacts_subset$csid)
# ID's for their contacts
two_week_contacts_subset$contact=paste(two_week_contacts_subset$csid,two_week_contacts_subset$contact_id,sep='.')

# --- Study participant, age, no. of contacts, no. of nHH contacts
study_participants_demo<-data.frame(id=study_participants,
                                    age=two_week_contacts_subset$age_class_part[match(study_participants,two_week_contacts_subset$csid)],
                                    no.con= sapply(study_participants, function(p)sum(two_week_contacts_subset$csid==p)),
                                    no.nHH.con=sapply(study_participants, function(p)sum(two_week_contacts_subset$csid==p & two_week_contacts_subset$share_hh=="non-HH"))) 

# ==============================================================================
# --- Set up the outbreak
# This function is similar to the 'outbreak_setup' in the 'ringbp' library except
# that the case id is sampled from the participants in the contact study
outbreak_setup_ext<-function (num.initial.cases,case_ids, incfn, delayfn, k, prop.asym) 
{
  inc_samples <- incfn(num.initial.cases) # determine incubation period for each case
  case_data <- data.table(exposure = rep(0, num.initial.cases), # assign exposure time of 0 for all index cases
                          asym = purrr::rbernoulli(num.initial.cases, prop.asym), # determine which index cases were asymptomatic
                          caseid = sample(case_ids,num.initial.cases,replace = F), # assign case id
                          infector = 0, # assign infection source as 0(unknown)
                          missed = TRUE, 
                          onset = inc_samples, # assign symptom onset day based on incubation periods
                          new_cases = NA) # no transmission events as yet
  case_data <- case_data[, `:=`(isolated_time, onset + delayfn(1))][, 
                                                                    `:=`(isolated, FALSE)] # determine isolation time
  case_data$isolated_time[case_data$asym] <- Inf # if case(s) were asymptomatic, then were never isolated
  return(case_data)
}
incfn <- dist_setup(dist_shape = 2.322737, dist_scale = 6.492272)
delay_shape = 2.5;delay_scale = 5
delayfn <- dist_setup(delay_shape, delay_scale)

case_data_setup <- outbreak_setup_ext(num.initial.cases = 3,study_participants,
                                  incfn,delayfn,k=1.95,prop.asym=0)
# --- Create the next generation of cases
# This function is an extension of the 'outbreak_step' in the 'ringbp' library
# the tries to account for clustering of contacts and hence depletion of local 
# pool of susceptibles
outbreak_step_ext<-function (case_data = NULL, disp.iso = NULL, disp.com = NULL, 
                             r0isolated = NULL, r0community = NULL, prop.asym = NULL, 
                             incfn = NULL, delayfn = NULL, prop.ascertain = NULL, k = NULL, 
                             quarantine = NULL) 
{
  
  vect_isTRUE <- function(x) {
    purrr::map_lgl(x, isTRUE)
  }
  vect_max <- function(x, y) {
    purrr::map2_dbl(x, y, max)
  }
  vect_min <- function(x, y) {
    purrr::map2_dbl(x, y, min)
  }
  
  # Assigning the number of new cases based on Ro and dispersion parameter
  case_data[, `:=`(new_cases, purrr::map2_dbl(ifelse(vect_isTRUE(isolated), 
                                                     disp.iso, disp.com), # choice of dispersion parameter if isolated
                                              ifelse(vect_isTRUE(isolated), r0isolated,
                                                     r0community),        # choice of mean Ro if isolated
                                              ~rnbinom(1, size = .x, mu = .y)))] # sample Ro for the case(s) 
  
  new_case_data <- case_data[new_cases > 0]
  total_new_cases <- case_data[, sum(new_cases), ]
  # If there are no new cases
  if (total_new_cases == 0) {
    case_data$isolated <- TRUE
    effective_r0 <- 0
    cases_in_gen <- 0
    out <- list(case_data, effective_r0, cases_in_gen)
    names(out) <- c("cases", "effective_r0", "cases_in_gen")
    return(out)
  }
  
  # Assign incubation periods for each new case
  inc_samples <- incfn(total_new_cases) 
  # Assign variables for the potential secondary cases
  prob_samples <- data.table(exposure = unlist(purrr::map2(new_case_data$new_cases, 
                                                           new_case_data$onset, function(x, y) {
                                                             inf_fn(rep(y, x), k)
                                                           })),
                             infector = unlist(purrr::map2(new_case_data$caseid,
                                                           new_case_data$new_cases, function(x, y) {
                                                             rep(x, y)
                                                           })),
                             infector_iso_time = unlist(purrr::map2(new_case_data$isolated_time,
                                                                    new_case_data$new_cases, function(x, y) {
                                                                      rep(x, as.integer(y))
                                                                    })),
                             infector_asym = unlist(purrr::map2(new_case_data$asym,
                                                                new_case_data$new_cases, function(x, y) {
                                                                  rep(x, y)
                                                                })),
                             asym = purrr::rbernoulli(n = total_new_cases, p = prop.asym), 
                             missed = purrr::rbernoulli(n = total_new_cases, p = 1 - 
                                                          prop.ascertain),
                             incubfn_sample = inc_samples,
                             isolated = FALSE, 
                             new_cases = NA)
  # If infection time was before source isolation time, then keep the new cases and assign an onset time
  prob_samples <- prob_samples[exposure < infector_iso_time][, 
                                                             `:=`(onset = exposure + incubfn_sample)]
  # If no secondary case exposures occured pre-isolation
  if (nrow(prob_samples)==0) {
    case_data$isolated <- TRUE
    effective_r0 <- 0
    cases_in_gen <- 0
    out <- list(case_data, effective_r0, cases_in_gen)
    names(out) <- c("cases", "effective_r0", "cases_in_gen")
    return(out)
  }
  # Assign case id according to local contact network
  #-----------------------------------------------------------------------------
  case.id.assignement<-function(infector.i){
    # --- Contacts
    # if case is a study participant
    if(any(study_participants==infector.i)){
      # contacts
      con=two_week_contacts_subset$contact[which(two_week_contacts_subset$csid==infector.i)]
      # Weight of contact
      weight=two_week_contacts_subset$weight[which(two_week_contacts_subset$csid==infector.i)]
      
    }else{
      # --- contacts
      # study participant linked to contact
      part.link=floor(as.numeric(infector.i))
      # - HH contacts
      hh.con=c(part.link,two_week_contacts_subset$contact[which(two_week_contacts_subset$csid==part.link)])
      # - Non-HH contacts (sample from pool of entire population)
      # Age of infector
      age.infector=two_week_contacts_subset$age_class_cont[which(two_week_contacts_subset$contact==infector.i)]
      # Study participants of similar age
      age.similar=which(study_participants_demo$age==age.infector)
      # Number of expected contacts based on age group
      nHH.con.infector=sample(study_participants_demo$no.nHH.con[age.similar],size = 1)
      # id of contacts
      nhh.con=sample(c(setdiff(study_participants,part.link),unique(two_week_contacts_subset$contact)),
                     size=nHH.con.infector,replace = F)
      # All contacts
      con=c(hh.con,nhh.con)
      # Weight of contact(randomly assign)
      weight=sample(two_week_contacts_subset$weight,size=length(con))
    }  
    
    # --- Susceptibility of contacts
    if(any(con %in% case_data$caseid)){
      weight=weight[-which(con %in% case_data$caseid)]
      con=con[-which(con %in% case_data$caseid)]
    }
    # --- Required number of secondary cases
    sec.cases=sum(prob_samples$infector==infector.i)
    
    # --- Weighted sampling of which of these contacts gets the infection
    # if there are susceptible contacts
    if(length(con)>0){
      # if pool of susceptible is greater than required secondary contacts
      if(sec.cases<length(con)){
        sec.case.actual=sample(x = con ,size = sec.cases,replace = F,
                               prob = weight)
        #prob_samples$caseid[prob_samples$infector==infector.i]<-sec.case.actual
        sec.case.id<-sec.case.actual
      }else{
        x=which(prob_samples$infector==infector.i)
        #keep=x[1:length(con)];prob_samples$caseid[keep]<-con;prob_samples$caseid[x!=keep]<-NaN
        sec.case.id<-rep(NaN,sec.cases);sec.case.id[1:length(con)]<-con
      }
    }else{
      #prob_samples$caseid[prob_samples$infector==infector.i]<- NaN
      sec.case.id<-rep(NaN,sec.cases)
    }
    return(sec.case.id) 
  }
  
  prob_samples$caseid<-c(unlist(sapply(unique(prob_samples$infector),case.id.assignement)))
  
  # If after considering local network of susceptibles there are no secondary cases
  if (all(is.na(prob_samples$caseid))) {
    case_data$isolated <- TRUE
    effective_r0 <- 0
    cases_in_gen <- 0
    out <- list(case_data, effective_r0, cases_in_gen)
    names(out) <- c("cases", "effective_r0", "cases_in_gen")
    return(out)
  }
  # If some sec cases did not occur due to lack of susceptible in the local network
  prob_samples=prob_samples[!is.na(as.numeric(prob_samples$caseid)),]
  # If caseIDs were repeated, i.e. allocated to more than 1 source, remove the 
  # duplicate sec case.
  prob_samples<-prob_samples[!duplicated(prob_samples$caseid),]
  #-----------------------------------------------------------------------------
  # If the infector was asymptomatic, then all the secondary cases were missed (i.e. not traced)
  prob_samples$missed[vect_isTRUE(prob_samples$infector_asym)] <- TRUE
  # Assign an isolation time for each new case
  prob_samples[, `:=`(isolated_time, ifelse(vect_isTRUE(asym), 
                                            Inf, ifelse(vect_isTRUE(missed), onset + delayfn(1), 
                                                        ifelse(!vect_isTRUE(rep(quarantine, total_new_cases)), 
                                                               vect_min(onset + delayfn(1), vect_max(onset, 
                                                                                                     infector_iso_time)), infector_iso_time))))]
  # remove unwanted variables
  prob_samples[, `:=`(c("incubfn_sample", "infector_iso_time", 
                        "infector_asym"), NULL)]
  
  # Other variables
  cases_in_gen <- nrow(prob_samples)
  effective_r0 <- nrow(prob_samples)/nrow(case_data[!vect_isTRUE(case_data$isolated)])
  case_data$isolated <- TRUE
  case_data <- data.table::rbindlist(list(as.data.table(case_data), as.data.table(prob_samples)), 
                                     use.names = TRUE)
  out <- list(case_data, effective_r0, cases_in_gen)
  names(out) <- c("cases", "effective_r0", "cases_in_gen")
  return(out)
  
}

# Generation 1
case_data <- outbreak_step_ext(case_data =case_data_setup,disp.iso =1,disp.com =0.16,r0isolated =0,
                           r0community =2.5,prop.asym = 0,incfn =incfn,delayfn =delayfn,
                           prop.ascertain =0,k=1.95,quarantine =FALSE)
case_data
# Generation 2
case_data <- outbreak_step_ext(case_data =case_data$cases,disp.iso =1,disp.com =0.16,r0isolated =0,
                           r0community =2.5,prop.asym = 0,incfn =incfn,delayfn =delayfn,
                           prop.ascertain =0,k=1.95,quarantine =FALSE)
case_data
# Generation 3
case_data <- outbreak_step_ext(case_data =case_data$cases,disp.iso =1,disp.com =0.16,r0isolated =0,
                           r0community =2.5,prop.asym = 0,incfn =incfn,delayfn =delayfn,
                           prop.ascertain =0,k=1.95,quarantine =FALSE)
# Generation 4
case_data <- outbreak_step_ext(case_data =case_data$cases,disp.iso =1,disp.com =0.16,r0isolated =0,
                           r0community =2.5,prop.asym = 0,incfn =incfn,delayfn =delayfn,
                           prop.ascertain =0,k=1.95,quarantine =FALSE)
case_data
# --- Run a single instance of the branching process model
# This function is more of less the same function in the 'ringbp' library
outbreak_model_ext<-function (num.initial.cases = NULL, prop.ascertain = NULL, cap_max_days = NULL, 
                          cap_cases = NULL, r0isolated = NULL, r0community = NULL, 
                          disp.iso = NULL, disp.com = NULL, k = NULL, delay_shape = NULL, 
                          delay_scale = NULL, prop.asym = NULL, quarantine = NULL) 
{
  incfn <- dist_setup(dist_shape = 2.322737, dist_scale = 6.492272)
  delayfn <- dist_setup(delay_shape, delay_scale)
  total.cases <- num.initial.cases
  latest.onset <- 0
  extinct <- FALSE
  case_data <- outbreak_setup_ext(num.initial.cases = num.initial.cases,case_ids = study_participants, 
                              incfn = incfn, prop.asym = prop.asym, delayfn = delayfn, 
                              k = k)
  effective_r0_vect <- c()
  cases_in_gen_vect <- c()
  while (latest.onset < cap_max_days & total.cases < cap_cases & 
         !extinct) {
    out <- outbreak_step_ext(case_data = case_data, disp.iso = disp.iso, 
                             disp.com = disp.com, r0isolated = r0isolated, r0community = r0community, 
                             incfn = incfn, delayfn = delayfn, prop.ascertain = prop.ascertain, 
                             k = k, quarantine = quarantine, prop.asym = prop.asym)
    case_data <- out[[1]]
    effective_r0_vect <- c(effective_r0_vect, out[[2]])
    cases_in_gen_vect <- c(cases_in_gen_vect, out[[3]])
    total.cases <- nrow(case_data)
    latest.onset <- max(case_data$onset)
    extinct <- all(case_data$isolated)
  }
  weekly_cases <- case_data[, `:=`(week, floor(onset/7))][, 
                                                          .(weekly_cases = .N), by = week]
  max_week <- floor(cap_max_days/7)
  missing_weeks <- (0:max_week)[!(0:max_week %in% weekly_cases$week)]
  if (length(missing_weeks > 0)) {
    weekly_cases <- data.table::rbindlist(list(weekly_cases, 
                                               data.table(week = missing_weeks, weekly_cases = 0)))
  }
  weekly_cases <- weekly_cases[order(week)][, `:=`(cumulative, 
                                                   cumsum(weekly_cases))]
  weekly_cases <- weekly_cases[week <= max_week]
  weekly_cases <- weekly_cases[, `:=`(mean_effective_r0 = mean(effective_r0_vect, 
                                                               na.rm = TRUE),
                                      effective_r0_per_gen=list(effective_r0_vect),
                                      cases_per_gen = list(cases_in_gen_vect))]
  return(weekly_cases)
}
outbreak_data<-outbreak_model_ext(num.initial.cases = 5,prop.ascertain = 0.5,
                              cap_max_days = 365,cap_cases = 2000,r0isolated = 0.5,
                              r0community = 2.5,disp.iso = 1,disp.com = 0.16,
                              k = 0.7,delay_shape = 2.5,delay_scale = 5,
                              prop.asym = 0.3,quarantine = TRUE)
# --- Running multiple simulations
# This function is the same as the LSHTM one, just changed the name of the model
# function 
scenario_sim_ext<-function (n.sim = NULL, prop.ascertain = NULL, cap_max_days = NULL, 
                            cap_cases = NULL, r0isolated = NULL, r0community = NULL, 
                            disp.iso = NULL, disp.com = NULL, k = NULL, delay_shape = NULL, 
                            delay_scale = NULL, num.initial.cases = NULL, prop.asym = NULL, 
                            quarantine = NULL) 
{
  res <- purrr::map(.x = 1:n.sim, ~outbreak_model_ext(num.initial.cases = num.initial.cases, 
                                                  prop.ascertain = prop.ascertain, cap_max_days = cap_max_days, 
                                                  cap_cases = cap_cases, r0isolated = r0isolated, r0community = r0community, 
                                                  disp.iso = disp.iso, disp.com = disp.com, delay_shape = delay_shape, 
                                                  delay_scale = delay_scale, k = k, prop.asym = prop.asym, 
                                                  quarantine = quarantine))
  res <- data.table::rbindlist(res)
  res[, `:=`(sim, rep(1:n.sim, rep(floor(cap_max_days/7) + 
                                     1, n.sim))), ]
  return(res)
}

library(ringbp)
# --- For the same parameters, how do the range of results from our extension
# compare to the original model

# Our extension
time.ext<-system.time(res.ext <- scenario_sim_ext(n.sim = 1000,num.initial.cases = 5,cap_max_days = 365,
                                                  cap_cases = 2000,r0isolated = 0,r0community = 2.5,
                                                  disp.iso = 1,disp.com = 0.16,k = 0.7,delay_shape = 2.5,
                                                  delay_scale = 5,prop.asym = 0,prop.ascertain = 0))
# LSHTM original
time.lshtm<-system.time(res.lshtm <- scenario_sim(n.sim = 1000,num.initial.cases = 5,cap_max_days = 365,
                                                  cap_cases = 2000,r0isolated = 0,r0community = 2.5,
                                                  disp.iso = 1,disp.com = 0.16,k = 0.7,delay_shape = 2.5,
                                                  delay_scale = 5,prop.asym = 0,prop.ascertain = 0))




weekly.cases.ext=matrix(res.ext$weekly_cases,byrow = F,ncol=1000)
weekly.cases.lshtm=matrix(res.lshtm$weekly_cases,byrow = F,ncol=1000)
# Weekly cases
matplot(res.ext$weekweekly.cases.lshtm,type='l',col='grey',lty=1,lwd=2)
matlines(weekly.cases.ext,col='pink',lwd=2)
legend('topright',col=c('grey','pink'),lty=1,lwd=2,legend = c('lshtm','extension'),
       bty='n') 
# cases per generation
ext.gen=res.ext$cases_per_gen[[1]];lshtm.gen=res.lshtm$cases_per_gen[[1]]
plot(ext.gen,type='l',lwd=2,col='pink',ylim=c(0,max(ext.gen,lshtm.gen)))
lines(lshtm.gen,lwd=2)
