# Outbreak functions
# as informed by contact data
# The popolation being modelled is closed and includes household structure

#ringbp package includes the functions LSHTM uses as well as auxiliary functions that we make use of
require(ringbp)


# --- Set up the outbreak
# This function is similar to the 'outbreak_setup' in the 'ringbp' library except
# that the case id is sampled from the participants in the contact study
outbreak_setup_ext<-function (num.initial.cases,case_ids, incfn, delayfn, k, prop.asym) 
{
  inc_samples <- incfn(num.initial.cases) # determine incubation period for each case
  case_data <- data.table(exposure = rep(0, num.initial.cases), # assign exposure time of 0 for all index cases
                          asym = purrr::rbernoulli(num.initial.cases, prop.asym), # determine which index cases were asymptomatic
                          caseid = case_ids, # assign case id
                          infector = 0, # assign infection source as 0(unknown)
                          missed = TRUE, 
                          onset = inc_samples, # assign symptom onset day based on incubation periods
                          new_cases = NA) # no transmission events as yet
  case_data <- case_data[, `:=`(isolated_time, onset + delayfn(1))][, 
                                                                    `:=`(isolated, FALSE)] # determine isolation time
  case_data$isolated_time[case_data$asym] <- Inf # if case(s) were asymptomatic, then were never isolated
  return(case_data)
}


# --- Create the next generation of cases
# This function is an extension of the 'outbreak_step' in the 'ringbp' library
# to include contact structure
# Keep track of depleteing susc pool
# Allow two levels of isolation, HH contacts and all contacts
outbreak_step_ext<-function (case_data = NULL, disp.iso = NULL, disp.com = NULL, 
                             r0isolated = NULL, r0community = NULL, prop.asym = NULL, 
                             incfn = NULL, delayfn = NULL, prop.ascertain = NULL, k = NULL, 
                             quarantine = NULL,
                             HH_trace=NULL,
                             contact_data=NULL) 
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
  
  # Get contact data
  study_participants <- contact_data[[1]]
  nhh.study_participants <- contact_data[[2]]
  two_week_contacts_subset <- contact_data[[3]]
  
  #Remove any study participants already infected
  study_participants <- study_participants[!csid%in%case_data[is.na(new_cases),caseid]]
  
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
    # cases_in_gen <- 0
    out <- list(case_data, effective_r0,0,0,0,0,0)
    names(out) <- c("cases", "effective_r0","dep.iso","dep.con","dep.dup","dep.HH","dep.nHH")
    return(out)
  }
  
  # Assign incubation periods for each new case
  inc_samples <- incfn(total_new_cases) 
  # Assign variables for the potential secondary cases
  prob_samples <- data.table(exposure = unlist(purrr::pmap(list(new_case_data$new_cases,
                                                                new_case_data$onset,
                                                                new_case_data$exposure),
                                                           function(x, y, z) {
                                                             inf_fn(rep(y-z, x), k)+z
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
  
  dep.iso <- total_new_cases-nrow(prob_samples) ## keep track of  depletion by infector isolation
  
  # If no secondary case exposures occured pre-isolation
  if (nrow(prob_samples)==0) {
    case_data$isolated <- TRUE
    effective_r0 <- 0
    # cases_in_gen <- 0
    out <- list(case_data, effective_r0,0,0,0,0,0)
    names(out) <- c("cases", "effective_r0","dep.iso","dep.con","dep.dup","dep.HH","dep.nHH")
    return(out)
  }
  
  ##Sample contacts from contact data
  #Get number of cases by infector
  new_case_data2 <- prob_samples[,.N,by=infector]
  #Add their contacts and weights
  new_case_data2 <- two_week_contacts_subset[,`:=`(infector=csid)][new_case_data2,on="infector"]
  #Get rid of any contacts that are already infected
  new_case_data2[,`:=`(inf_inds=list(!unlist(contact_id)%in%case_data$caseid )),by=infector][
    ,`:=`(contact_id=list(unlist(contact_id)[unlist(inf_inds)]),
          weight=list(unlist(weight)[unlist(inf_inds)])),by=infector]
  #Count the number of contacts available
  new_case_data2[,`:=`(con=length(unlist(contact_id))),by=infector]
  
  #Keep track of contact depletion if no contacts are available and get rid of any infectors that have no contacts
  dep.con <- new_case_data2[con==0,N] %>% sum
  prob_samples <- prob_samples[!infector%in%new_case_data2[con==0,csid],]
  new_case_data2 <- new_case_data2[!con==0,]
  
  #Sample infecteds from the contacts available by weight
  z <- new_case_data2[,`:=`(inf_id=list(sample(x = unlist(contact_id),
                                               size = min(N,con),
                                               prob = unlist(weight)))),by=infector][
                                                 ,.(inf_id=unlist(inf_id)),by=infector]
  
  ##Keep track of contact depletion and get rid of infections that didn't occur due to lack of contacts
  con.dep <- new_case_data2[,.(diff=max(N-con,0)),by=.(infector,N)]
  dep.con <- dep.con+sum(con.dep$diff)
  con.dep[,`:=`(ind1=1+cumsum(N)-N)]
  con.dep[,`:=`(ind2=ind1+diff-1)]
  dep.inds <- unlist(con.dep[,`:=`(inds=ifelse(ind1>ind2,list(NA),list(ind1:ind2))),by=infector][,.(unlist(inds[!is.na(inds)]))])
  if (!is.null(dep.inds)) {
    prob_samples <- prob_samples[!dep.inds,]
  }
  
  ##Remove and count duplicates
  dep.dup <- 0
  if (any(duplicated(z$inf_id))) {
    dep.dup <- duplicated(z$inf_id) %>% sum
    prob_samples <- prob_samples[!duplicated(z$inf_id)]
    z <- z[!duplicated(inf_id)]
  }
  
  ##Remove any duplicates 
  dep.nHH <- 0
  if (any(duplicated(z[grepl("nhh",inf_id),inf_id]))) {
    dep.nHH <- duplicated(z[grepl("nhh",inf_id),inf_id]) %>% sum
    z[grepl("nhh",inf_id)] <- z[grepl("nhh",inf_id)][!duplicated(z[grepl("nhh",inf_id),inf_id])]
  }
  
  #Now need to match the non-HH contacts to other study participants
  if (length(z[grepl("nhh",inf_id),inf_id])>0) {
    #Study particiapnts by age to sample from
    study_participants.age <- study_participants[,.(csid=list(csid)),by=age]
    ##Non HH contacts ages
    nhh.inf.ages <- nhh.study_participants[csid%in%z[grepl("nhh",inf_id),inf_id],]
    #Sample
    new.nhh.ids <- nhh.inf.ages[study_participants.age,nomatch=0L,on="age"][,.(new.csid=sample(unlist(i.csid),1)),by=csid]
    z[which(inf_id%in%new.nhh.ids$csid),`:=`(inf_id=new.nhh.ids[,new.csid])]
  }
  
  ##Assign all the new caseids
  prob_samples$caseid <- z$inf_id
  
  #  Add the type of contact (HH or non-HH)
  prob_samples[,`:=`(case.link=floor(as.numeric(caseid)),
                     infector.link=floor(as.numeric(infector)))][
                       ,`:=`(type,
                     ifelse(case.link==infector.link,
                            "HH","non-HH")),
               by=.(caseid)]
  
  
  # If after considering local network of susceptibles there are no secondary cases
  if (all(is.na(prob_samples$caseid))) {
    case_data$isolated <- TRUE
    effective_r0 <- 0
    # cases_in_gen <- 0
    out <- list(case_data, effective_r0,0,0,0,0,0)
    names(out) <- c("cases", "effective_r0","dep.iso","dep.con","dep.dup","dep.HH","dep.nHH")
    # out <- list(case_data, effective_r0)
    # names(out) <- c("cases", "effective_r0")
    return(out)
  }
  
  ##Different depletion scenarios
  ##Infector isolated
  dep.iso = dep.iso/total_new_cases
  #Contacts
  dep.con = dep.con/total_new_cases
  #Dupplicate assignment
  dep.dup = (dep.dup+nrow(prob_samples[duplicated(caseid),]))/total_new_cases
  # HH already infected
  dep.HH = nrow(prob_samples[type=="HH" & caseid%in%case_data$caseid])/total_new_cases
  # nHH already infected + those that were chosen twice
  dep.nHH = (dep.nHH+nrow(prob_samples[type=="non-HH" & caseid%in%case_data$caseid]))/total_new_cases
  
  ##Now get rid of duplicates
  # If caseIDs were repeated, i.e. allocated to more than 1 source, remove the 
  # duplicate sec case.
  prob_samples <- prob_samples[!duplicated(caseid),]
  #Remove assigned ids that are already infected (HH and non-HH)
  prob_samples <- prob_samples[!caseid%in%case_data$caseid]
  
  if (all(is.na(prob_samples$caseid))) {
    case_data$isolated <- TRUE
    effective_r0 <- 0
    # cases_in_gen <- 0
    out <- list(case_data, effective_r0,0,0,0,0,0)
    names(out) <- c("cases", "effective_r0","dep.iso","dep.con","dep.dup","dep.HH","dep.nHH")
    # out <- list(case_data, effective_r0)
    # names(out) <- c("cases", "effective_r0")
    return(out)
  }
  
  # --- If HH_tracing, change missed status for HH contacts to FALSE. Assuming all HH contacts
  # can be traced
  if (HH_trace) prob_samples$missed[prob_samples$type=='HH']<-FALSE
  # ## If quarantine, we don't miss
  # if (quarantine)
  # If the infector was asymptomatic, then all the secondary cases were missed (i.e. not traced)
  prob_samples$missed[vect_isTRUE(prob_samples$infector_asym)] <- TRUE
  # Assign an isolation time for each new case
  prob_samples[, `:=`(isolated_time, ifelse(vect_isTRUE(asym & !rep(quarantine, nrow(prob_samples))), 
                                            Inf, ifelse(vect_isTRUE(missed), 
                                                        onset + delayfn(1), 
                                                        ifelse(!vect_isTRUE(rep(quarantine, nrow(prob_samples))), 
                                                               vect_min(onset + delayfn(1), 
                                                                        vect_max(onset,
                                                                                 infector_iso_time)), 
                                                               infector_iso_time))))]
  # remove unwanted variables
  prob_samples[, `:=`(c("incubfn_sample", "infector_iso_time", 
                        "infector_asym","case.link","infector.link"), NULL)]
  
  # Other variables
  # cases_in_gen <- nrow(prob_samples)
  effective_r0 <- nrow(prob_samples)/nrow(case_data[!vect_isTRUE(case_data$isolated)])
  case_data$isolated <- TRUE
  case_data <- data.table::rbindlist(list(as.data.table(case_data), as.data.table(prob_samples)), 
                                     use.names = TRUE,fill=TRUE)
  
  case_data[is.na(new_cases)]$caseid
  
  out <- list(case_data, effective_r0,dep.iso,dep.con,dep.dup,dep.HH,dep.nHH)
  names(out) <- c("cases", "effective_r0","dep.iso","dep.con","dep.dup","dep.HH","dep.nHH")
  return(out)
  
}

# --- Run a single instance of the branching process model
# This function is an extension of the 'outbreak_model' in the 'ringbp' library
# We include counting of contacts and depletion of susceptibles
outbreak_model_ext<-function (num.initial.cases = NULL, prop.ascertain = NULL, cap_max_days = NULL, 
                              cap_cases = NULL, cap_max_hhcon=NULL, cap_max_nhhcon=NULL, 
                              r0isolated = NULL, r0community = NULL, 
                              disp.iso = NULL, disp.com = NULL, k = NULL, delay_shape = NULL, 
                              delay_scale = NULL, prop.asym = NULL, quarantine = NULL,
                              HH_trace=NULL,
                              contact_data=NULL,t=NULL,p.urban=NULL,freq =NULL,
                              pop=NULL) 
{
  # Bootstrap contacts for t=14 days
  # CD <- count.contacts.full(contact_data,t=14)
  # time.cd <- system.time(CD <- count.contacts3(contact_data,t=t,p.urban = p.urban,freq=freq,pop=pop))
  CD <- contact_data
  study_participants <- CD[[1]]
  if ("caseid"%in%names(CD[[1]])) {
    setnames(CD[[1]],old = "caseid",new = "csid")
  }
  # Set up params
  incfn <- dist_setup(dist_shape = 2.322737, dist_scale = 6.492272)
  delayfn <- dist_setup(delay_shape, delay_scale)
  total.cases <- num.initial.cases
  total.hh.con <- 0
  total.nhh.con <- 0
  latest.onset <- 0
  extinct <- FALSE
  case_data <- outbreak_setup_ext(num.initial.cases = num.initial.cases,case_ids = sample(CD[[1]]$csid,num.initial.cases,replace = F), 
                                  incfn = incfn, prop.asym = prop.asym, delayfn = delayfn, 
                                  k = k)
  # effective_r0_vect <- 0
  effective_r0_vect <- c()
  dep.iso_vect <- c()
  dep.con_vect <- c()
  dep.dup_vect <- c()
  dep.HH_vect <- c()
  dep.nHH_vect <- c()
  cases_in_gen_vect <- c()
  # gen <- 1
  while (!extinct & latest.onset < cap_max_days & total.cases < cap_cases & 
         total.hh.con < cap_max_hhcon &  total.nhh.con < cap_max_nhhcon) {
    out <- outbreak_step_ext(case_data = case_data, disp.iso = disp.iso, 
                             disp.com = disp.com, r0isolated = r0isolated, r0community = r0community, 
                             incfn = incfn, delayfn = delayfn, prop.ascertain = prop.ascertain, 
                             k = k, quarantine = quarantine, prop.asym = prop.asym,
                             HH_trace=HH_trace,
                             contact_data=CD)
    
    
    ##Can remove all the cases that were exposed after cap_max_days, as these can only infect cases we're not interested in
    case_data <- out$cases[exposure<cap_max_days]
    # case_data <- out$cases
    
    #Add to our vectors
    effective_r0_vect <- c(effective_r0_vect, out$effective_r0)
    dep.iso_vect <- c(dep.iso_vect, out$dep.iso)
    dep.con_vect <- c(dep.con_vect, out$dep.con)
    dep.dup_vect <- c(dep.dup_vect, out$dep.dup)
    dep.HH_vect <- c(dep.HH_vect, out$dep.HH)
    dep.nHH_vect <- c(dep.nHH_vect, out$dep.nHH)
    cases_in_gen_vect <- c(cases_in_gen_vect, nrow(out$cases[is.na(new_cases)]))
    
    # total.cases <- nrow(case_data)
    total.cases <- nrow(case_data[onset<cap_max_days]) ## Want total cases up until cap time!
    # total.cases <- nrow(case_data)
    # latest.onset <- max(case_data$onset)
    # latest.onset <- min(case_data[is.na(new_cases),onset]) #earliest onset of latest generation, to make sure all the cases are included??
    latest.onset <- min(case_data[is.na(new_cases),exposure]) #earliest onset of latest generation, to make sure all the cases are included??
    extinct <- all(out$cases$isolated)
    # latest.onset
    # total.hh.con <- sum(case_data$hh.con,na.rm =TRUE)*HH_trace
    # total.nhh.con <- sum(case_data$nhh.con,na.rm =TRUE)*prop.ascertain
    # case_data[order(exposure)]
    # gen <- gen+1
  }
  # gen
  extinct
  cases_in_gen_vect
  cumsum(cases_in_gen_vect)
  length(effective_r0_vect)
  
  if (nrow(case_data)==num.initial.cases) {
    case_data[,`:=`(hh_cases=0,type=NA)]
  }
  
  # Assign HH and non-HH contacts
  setnames(study_participants,old = "csid",new = "caseid")
  case_data <- case_data[study_participants,nomatch=0L,on="caseid"]
  ##If a contact was a HH contact, then tracing/quarantine doesn't need to be reapeated for hh contacts
  case_data[type=="HH",`:=`(hh.con=0)]
  
  #Summing by week
  weekly_cases <- case_data[, `:=`(week, floor(onset/7))][, 
                                                          .(weekly_cases = .N,
                                                            hh.con=sum(hh.con,na.rm = TRUE),
                                                            nhh.con=sum(nhh.con,na.rm = TRUE),
                                                            hh_cases=sum(type=='HH',na.rm = TRUE)), by = week]
  max_week <- floor(cap_max_days/7)
  missing_weeks <- (0:max_week)[!(0:max_week %in% weekly_cases$week)]
  if (length(missing_weeks > 0)) {
    weekly_cases <- data.table::rbindlist(list(weekly_cases, 
                                               data.table(week = missing_weeks, weekly_cases = 0)),fill=TRUE)
  }
  weekly_cases <- weekly_cases[order(week)][, `:=`(cumulative, 
                                                   cumsum(weekly_cases))]
  weekly_cases <- weekly_cases[week <= max_week]
  data.table(gen=1:length(effective_r0_vect),
             cases_in_gen=cases_in_gen_vect,
             effective_r0=effective_r0_vect,
             dep.iso=dep.iso_vect,
             dep.con=dep.con_vect,
             dep.dup=dep.dup_vect,
             dep.HH=dep.HH_vect,
             dep.nHH=dep.nHH_vect)
  
  # weekly_cases$mean_effective_r0 <- effective_r0_vect/gen
  # weekly_cases <- weekly_cases[, `:=`(mean_effective_r0 = effective_r0_vect/gen,
  #                                     effective_r0_per_gen=list(effective_r0_vect),
  #                                     cases_per_gen = list(cases_in_gen_vect))]
  return(list(weekly_cases=weekly_cases,
              gen_vect=data.table(gen=1:length(effective_r0_vect),
                                  cases_in_gen=cases_in_gen_vect,
                                  effective_r0=effective_r0_vect,
                                  dep.iso=dep.iso_vect,
                                  dep.con=dep.con_vect,
                                  dep.dup=dep.dup_vect,
                                  dep.HH=dep.HH_vect,
                                  dep.nHH=dep.nHH_vect)))
}


# --- Running multiple simulations in parallel
scenario_sim_ext_parallel<-function (n.sim = NULL, prop.ascertain = NULL, cap_max_days = NULL, 
                                     cap_cases = NULL, cap_max_hhcon=NULL, cap_max_nhhcon=NULL, 
                                     r0isolated = NULL, r0community = NULL, 
                                     disp.iso = NULL, disp.com = NULL, k = NULL, delay_shape = NULL, 
                                     delay_scale = NULL, num.initial.cases = NULL, prop.asym = NULL, 
                                     quarantine = NULL,HH_trace = NULL,
                                     contact_data=NULL,t=NULL,p.urban=NULL,freq = NULL,
                                     pop=NULL) 
{
  
  res <- try(parallel::mclapply(1:n.sim,
                            function(x) {
                              outbreak_model_ext(num.initial.cases = num.initial.cases, 
                                                 prop.ascertain = prop.ascertain, cap_max_days = cap_max_days, 
                                                 cap_cases = cap_cases, cap_max_hhcon=cap_max_hhcon, cap_max_nhhcon=cap_max_nhhcon, 
                                                 r0isolated = r0isolated, r0community = r0community, 
                                                 disp.iso = disp.iso, disp.com = disp.com, delay_shape = delay_shape, 
                                                 delay_scale = delay_scale, k = k, prop.asym = prop.asym, 
                                                 quarantine = quarantine,HH_trace = HH_trace,
                                                 contact_data=contact_data,p.urban = p.urban,
                                                 t=t,freq = freq,
                                                 pop=pop)
                            }
                            ,mc.cores = parallel::detectCores())
             )
  
  
  res.week <- data.table(bind_rows(lapply(res,function(x) x$weekly_cases)))
  res.week
  res.week$sim <- rep(1:n.sim,each=max(res.week$week)+1)
  
  res.gen <- data.table(bind_rows(lapply(res,function(x) x$gen_vect)))
  
  res.gen$sim <- 1
  sim <- 1
  for (i in 2:nrow(res.gen)) {
    if (res.gen$gen[i]<=res.gen$gen[i-1]) sim=sim+1
    res.gen$sim[i] <- sim
  }
  
  return(list(res.week=res.week,res.gen=res.gen))
}
