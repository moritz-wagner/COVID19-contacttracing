# Outbreak functions
# as informed by contact data
# The popolation being modelled is closed and has household structure
require(ringbp)


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


# --- Create the next generation of cases
# This function is an extension of the 'outbreak_step' in the 'ringbp' library
# the tries to account for clustering of contacts and hence depletion of local 
# pool of susceptibles
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
  two_week_contacts_subset <- contact_data[[2]]
  # study_participants_demo <- contact_data[[3]]
  
  
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
    out <- list(case_data, effective_r0,0,0,0,0)
    names(out) <- c("cases", "effective_r0","dep.con","dep.dup","dep.HH","dep.nHH")
    return(out)
  }
  
  # new_case_data$exposure
  # unlist(purrr::map2(new_case_data$new_cases,
  #                    new_case_data$onset,
  #                    function(x, y) {
  #                      inf_fn(rep(y, x), k)
  #                    })) %>% sort
  # unlist(purrr::pmap(list(new_case_data$new_cases,
  #                         new_case_data$onset,
  #                         new_case_data$exposure),
  #                    function(x, y, z) {
  #                      inf_fn(rep(y-z, x), k)+z
  #                    })) %>% sort
  # 
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
  # If no secondary case exposures occured pre-isolation
  if (nrow(prob_samples)==0) {
    case_data$isolated <- TRUE
    effective_r0 <- 0
    # cases_in_gen <- 0
    out <- list(case_data, effective_r0,0,0,0,0)
    names(out) <- c("cases", "effective_r0","dep.con","dep.dup","dep.HH","dep.nHH")
    return(out)
  }
  # Assign case id according to local contact network
  #-----------------------------------------------------------------------------
  # Function for limiting secondary cases according to local contact network
  case.id.assignement<-function(infector.i,inf){
    # --- Contacts
    # if case is a study participant
    
    if(any(study_participants$csid==infector.i)){
      # contacts
      hh.con <- two_week_contacts_subset[csid==infector.i & share_hh=="HH",contact]
      nhh.con <- two_week_contacts_subset[csid==infector.i & share_hh!="HH",contact]
      con=two_week_contacts_subset$contact[which(two_week_contacts_subset$csid==infector.i)]
      # Contact type
      # con.type=two_week_contacts_subset$share_hh[which(two_week_contacts_subset$csid==infector.i)]
      # Weight of contact
      weight=two_week_contacts_subset$weight[which(two_week_contacts_subset$csid==infector.i)]
      
    }else{
      # --- contacts
      # study participant linked to contact
      part.link=floor(as.numeric(infector.i))
      # - HH contacts
      hh.con=two_week_contacts_subset[csid==part.link & share_hh=="HH",contact]
      hh.con=c(part.link,hh.con[hh.con!=infector.i])
      weight.hh.con = sample(two_week_contacts_subset$weight[two_week_contacts_subset$share_hh=='HH'],size = length(hh.con))
      # - Non-HH contacts (sample from pool of entire population)
      # Age of infector
      age.infector=two_week_contacts_subset$age_class_cont[which(two_week_contacts_subset$contact==infector.i)]
      
      # # Get a random participant id of same age and use their non-HH contacts and weights
      # system.time({
      # x=two_week_contacts_subset[csid==two_week_contacts_subset[
      #   age_class_part==age.infector,sample(csid,1)] & share_hh=="non-HH",.(contact,weight)]
      # 
      # nhh.con=x$contact
      # weight.nhh.con=x$weight
      # })
      # system.time({
      # Study participants of similar age
      age.similar=which(study_participants$age==age.infector & study_participants$csid!=part.link)
      # get a random participant id of same age
      nHH.con.infector=sample(study_participants$csid[age.similar],size = 1)
      # Count their non-HH contacts and weights
      x=which(two_week_contacts_subset$csid==nHH.con.infector & two_week_contacts_subset$share_hh=='non-HH')
      nhh.con=two_week_contacts_subset$contact[x]
      weight.nhh.con=two_week_contacts_subset$weight[x]
      # })
      # All contacts
      con=c(hh.con,nhh.con)
      # con.type=c(rep("HH",length(hh.con)),rep("non-HH",length(hh.con)))
      weight=c(weight.hh.con,weight.nhh.con)
    }
    # Keep track of size of the contact network
    hh.con.count <- length(hh.con)
    nhh.con.count <- length(nhh.con)
    
    if (inf) {
      
      # --- Required number of secondary cases
      sec.cases=sum(prob_samples$infector==infector.i)
      
      # --- Weighted sampling of which of these contacts gets the infection
      # if there are susceptible contacts
      if(length(con)>0){
        # if pool of susceptible is greater than required secondary contacts
        if(sec.cases<length(con)){
          sec.case.id<-sample(x = con ,size = sec.cases,replace = F,
                              prob = weight)
          #sec.case.type<-con.type[which(con %in% sec.case.id)]
        }else{
          x=which(prob_samples$infector==infector.i)
          sec.case.id<-rep(NaN,sec.cases);sec.case.id[1:length(con)]<-con
          #sec.case.type<-rep(NaN,sec.cases);sec.case.type[1:length(con)]<-con.type
        }
      }else{
        sec.case.id<-rep(NaN,sec.cases); #sec.case.type<-rep(NaN,sec.cases)
      }
      
      # --- Susceptibility of contacts, remove those that are already infected
      # hh.dep=0
      # nhh.dep=0
      # if(any(sec.case.id %in% case_data$caseid)){
      #   # hh.dep=sum(hh.con[hh.con%in%sec.case.id] %in% case_data$caseid)
      #   # nhh.dep=sum(nhh.con[nhh.con%in%sec.case.id] %in% case_data$caseid)
      #   # sec.case.id=sec.case.id[!sec.case.id%in%case_data$caseid]
      #   sec.case.id[sec.case.id%in%case_data$caseid] <- NaN
      #   
      #   # weight=weight[-which(con %in% case_data$caseid)]
      #   # # con.type=con.type[-which(con %in% case_data$caseid)]
      #   # con=con[-which(con %in% case_data$caseid)]
      # }
      # Keep track of the size of the susceptible network
      #res=cbind(caseid=sec.case.id,type=sec.case.type)
      return(list(infector=infector.i,hh.con=hh.con.count,nhh.con=nhh.con.count,sec.case.id=sec.case.id) )
    }else{
      return(list(hh.con=hh.con.count,nhh.con=nhh.con.count)) 
    }
  }
  
  ## Get contact counts of non-infectors
  if (nrow(case_data[isolated==FALSE & !caseid%in%unique(prob_samples$infector)])>0) {
    case.assign.no.inf <- simplify2array(parallel::mclapply(case_data[isolated==FALSE & !caseid%in%unique(prob_samples$infector),caseid],
                                                            case.id.assignement,inf=FALSE))
    #Add contact counts 
    if (length(case.assign.no.inf)>0) {
      case_data[isolated==FALSE & !caseid%in%unique(prob_samples$infector),`:=`(hh.con,unlist(case.assign.no.inf[1,]))]
      case_data[isolated==FALSE & !caseid%in%unique(prob_samples$infector),`:=`(nhh.con,unlist(case.assign.no.inf[2,]))]
    }
  }
  
  ## Get ids and contact-contacts of infectors
  case.assign.inf <- parallel::mclapply(unique(prob_samples$infector),case.id.assignement,inf=TRUE)
  #Add contact counts to infectors
  case_data[isolated==FALSE & caseid%in%unique(prob_samples$infector),`:=`(hh.con,sapply(case.assign.inf,function(x) x$hh.con))]
  case_data[isolated==FALSE & caseid%in%unique(prob_samples$infector),`:=`(nhh.con,sapply(case.assign.inf,function(x) x$nhh.con))]
  
  # case_data[isolated==FALSE & caseid%in%unique(prob_samples$infector),`:=`(hh.dep,sapply(case.assign.inf,function(x) x$hh.dep))]
  # case_data[isolated==FALSE & caseid%in%unique(prob_samples$infector),`:=`(nhh.dep,sapply(case.assign.inf,function(x) x$nhh.dep))]
  
  # 
  # Sec case id
  # prob_samples$caseid <- NA
  #Assign case ids
  if (nrow(prob_samples)==sum(sapply(case.assign.inf,function(x) length(x$sec.case.id)))) {
    prob_samples$caseid <- unlist(sapply(case.assign.inf,function(x) {x$sec.case.id}))
  } else {
    prob_samples$caseid <- unlist(sapply(case.assign.inf,function(x) {
      ids <- c(x$sec.case.id,
               rep(NA,nrow(prob_samples[infector==x$infector])-length(x$sec.case.id)))
      if (length(ids)==1) {
        ids
      } else {
        sample(ids)
      }
    }
    ))
  }
  
  # system.time(unlist(sapply(unique(prob_samples$infector),case.id.assignement)))
  # system.time(unlist(unique(prob_samples$infector) %>% purrr::map(~case.id.assignement(.x))))
  
  # prob_samples$caseid<-unlist(parallel::mclapply(unique(prob_samples$infector),case.id.assignement))
  
  # profvis({
  #   for (id in unique(prob_samples$infector)) {
  #     x <- case.id.assignement(id)
  #   }
  # })
  #  Sec case contact type
  prob_samples$case.link=unlist(sapply((prob_samples$caseid),function(x)floor(as.numeric(x))))
  prob_samples$infector.link=unlist(sapply((prob_samples$infector),function(x)floor(as.numeric(x))))
  
  prob_samples[,`:=`(type,
                     ifelse(infector==infector.link,
                            as.character(two_week_contacts_subset$share_hh[which(two_week_contacts_subset$contact==caseid)]),
                            ifelse(case.link==infector.link,"HH","non-HH"))),
               by=.(caseid)]
  
  
  # If after considering local network of susceptibles there are no secondary cases
  if (all(is.na(prob_samples$caseid))) {
    case_data$isolated <- TRUE
    effective_r0 <- 0
    # cases_in_gen <- 0
    out <- list(case_data, effective_r0,0,0,0,0)
    names(out) <- c("cases", "effective_r0","dep.con","dep.dup","dep.HH","dep.nHH")
    # out <- list(case_data, effective_r0)
    # names(out) <- c("cases", "effective_r0")
    return(out)
  }
  
  ##Different depletion scenarios
  #Contacts
  dep.con = sum(is.na(as.numeric(prob_samples$caseid)))/nrow(prob_samples)
  #Dupplicate assignment
  dep.dup = sum(duplicated(prob_samples$caseid[!is.na(as.numeric(prob_samples$caseid))]))/nrow(prob_samples)
  # HH already infected
  dep.HH = nrow(prob_samples[type=="HH" & caseid%in%case_data$caseid])/nrow(prob_samples)
  # nHH already infected
  dep.nHH = nrow(prob_samples[type=="non-HH" & caseid%in%case_data$caseid])/nrow(prob_samples)
  
  ##Now get rid of them
  #Contacts
  prob_samples=prob_samples[!is.na(as.numeric(prob_samples$caseid)),]
  # If caseIDs were repeated, i.e. allocated to more than 1 source, remove the 
  # duplicate sec case.
  prob_samples<-prob_samples[!duplicated(prob_samples$caseid),]
  
  #Remove assigned ids that are already infected (HH and non-HH)
  prob_samples <- prob_samples[!caseid%in%case_data$caseid]
  
  
  #-----------------------------------------------------------------------------
  
  # --- If HH_tracing, change missed status for HH contacts to FALSE. Assuming all HH contacts
  # can be traced
  if (HH_trace) prob_samples$missed[prob_samples$type=='HH']<-FALSE
  # If the infector was asymptomatic, then all the secondary cases were missed (i.e. not traced)
  prob_samples$missed[vect_isTRUE(prob_samples$infector_asym)] <- TRUE
  # Assign an isolation time for each new case
  prob_samples[, `:=`(isolated_time, ifelse(vect_isTRUE(asym), 
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
  out <- list(case_data, effective_r0,dep.con,dep.dup,dep.HH,dep.nHH)
  names(out) <- c("cases", "effective_r0","dep.con","dep.dup","dep.HH","dep.nHH")
  return(out)
  
}

# --- Run a single instance of the branching process model
# This function is more of less the same function in the 'ringbp' library
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
  CD <- count.contacts(contact_data,t=t,p.urban = p.urban,freq=freq,pop=pop) 
  
  # Set up params
  incfn <- dist_setup(dist_shape = 2.322737, dist_scale = 6.492272)
  delayfn <- dist_setup(delay_shape, delay_scale)
  total.cases <- num.initial.cases
  total.hh.con <- 0
  total.nhh.con <- 0
  latest.onset <- 0
  extinct <- FALSE
  case_data <- outbreak_setup_ext(num.initial.cases = num.initial.cases,case_ids = CD[[1]]$csid, 
                                  incfn = incfn, prop.asym = prop.asym, delayfn = delayfn, 
                                  k = k)
  # effective_r0_vect <- 0
  effective_r0_vect <- c()
  dep.con_vect <- c()
  dep.dup_vect <- c()
  dep.HH_vect <- c()
  dep.nHH_vect <- c()
  cases_in_gen_vect <- c()
  # gen <- 1
  cases_in_gen_vect <- c()
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
    total.hh.con <- sum(case_data$hh.con,na.rm =TRUE)*HH_trace
    total.nhh.con <- sum(case_data$nhh.con,na.rm =TRUE)*prop.ascertain
    # case_data[order(exposure)]
    # gen <- gen+1
  }
  # gen
  extinct
  cases_in_gen_vect
  cumsum(cases_in_gen_vect)
  length(effective_r0_vect)
  
  if (nrow(case_data)==num.initial.cases) {
    case_data[,`:=`(hh.con=0,nhh.con=0,hh_cases=0,type=NA)]
  }
  
  # Can have all the HH contacts where infector was infected through the HH equal to zero, as they will have already been traced
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
  
  # weekly_cases$mean_effective_r0 <- effective_r0_vect/gen
  # weekly_cases <- weekly_cases[, `:=`(mean_effective_r0 = effective_r0_vect/gen,
  #                                     effective_r0_per_gen=list(effective_r0_vect),
  #                                     cases_per_gen = list(cases_in_gen_vect))]
  return(list(weekly_cases=weekly_cases,
              gen_vect=data.table(gen=1:length(effective_r0_vect),
                                  cases_in_gen=cases_in_gen_vect,
                                  effective_r0=effective_r0_vect,
                                  dep.con=dep.con_vect,
                                  dep.dup=dep.dup_vect,
                                  dep.HH=dep.HH_vect,
                                  dep.nHH=dep.nHH_vect)))
}


##Change up LSHTM's outbreak model to match ours with the same method for capping days
outbreak_model <- function (num.initial.cases = NULL, prop.ascertain = NULL, cap_max_days = NULL, 
                            cap_cases = NULL, r0isolated = NULL, r0community = NULL, 
                            disp.iso = NULL, disp.com = NULL, k = NULL, delay_shape = NULL, 
                            delay_scale = NULL, prop.asym = NULL, quarantine = NULL) 
{
  
  incfn <- dist_setup(dist_shape = 2.322737, dist_scale = 6.492272)
  delayfn <- dist_setup(delay_shape, delay_scale)
  total.cases <- num.initial.cases
  latest.onset <- 0
  extinct <- FALSE
  case_data <- outbreak_setup(num.initial.cases = num.initial.cases, 
                              incfn = incfn, prop.asym = prop.asym, delayfn = delayfn, 
                              k = k)
  effective_r0_vect <- c()
  cases_in_gen_vect <- c()
  while (latest.onset < cap_max_days & total.cases < cap_cases & 
         !extinct) {
    out <- outbreak_step(case_data = case_data, disp.iso = disp.iso, 
                         disp.com = disp.com, r0isolated = r0isolated, r0community = r0community, 
                         incfn = incfn, delayfn = delayfn, prop.ascertain = prop.ascertain, 
                         k = k, quarantine = quarantine, prop.asym = prop.asym)
    case_data <- out[[1]]
    effective_r0_vect <- c(effective_r0_vect, out[[2]])
    cases_in_gen_vect <- c(cases_in_gen_vect, out[[3]])
    total.cases <- nrow(case_data)
    latest.onset <- max(case_data$onset)
    # latest.onset <- min(case_data[is.na(new_cases),onset]) #earliest onset of latest generation, to make sure all the cases are included
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
  weekly_cases <- weekly_cases[, `:=`(effective_r0 = mean(effective_r0_vect, 
                                                          na.rm = TRUE), cases_per_gen = list(cases_in_gen_vect))]
  return(weekly_cases)
}





# --- Running multiple simulations
# This function is the same as the LSHTM one, just changed the name of the model
# function 
scenario_sim_ext<-function (n.sim = NULL, prop.ascertain = NULL, cap_max_days = NULL, 
                            cap_cases = NULL, r0isolated = NULL, r0community = NULL, 
                            disp.iso = NULL, disp.com = NULL, k = NULL, delay_shape = NULL, 
                            delay_scale = NULL, num.initial.cases = NULL, prop.asym = NULL, 
                            quarantine = NULL,HH_trace=NULL,contact_data=NULL) 
{
  res <- purrr::map(.x = 1:n.sim, ~outbreak_model_ext(num.initial.cases = num.initial.cases, 
                                                      prop.ascertain = prop.ascertain, cap_max_days = cap_max_days, 
                                                      cap_cases = cap_cases, cap_max_hhcon=cap_max_hhcon, cap_max_nhhcon=cap_max_nhhcon,
                                                      r0isolated = r0isolated, r0community = r0community, 
                                                      disp.iso = disp.iso, disp.com = disp.com, delay_shape = delay_shape, 
                                                      delay_scale = delay_scale, k = k, prop.asym = prop.asym, 
                                                      quarantine = quarantine,HH_trace = HH_trace,
                                                      contact_data=contact_data))
  
  res.week <- data.table(bind_rows(lapply(res,function(x) x$weekly_cases)))
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

# Same function but using foreach packacge to increase speed
scenario_sim_ext_parallel<-function (n.sim = NULL, prop.ascertain = NULL, cap_max_days = NULL, 
                                     cap_cases = NULL, cap_max_hhcon=NULL, cap_max_nhhcon=NULL, 
                                     r0isolated = NULL, r0community = NULL, 
                                     disp.iso = NULL, disp.com = NULL, k = NULL, delay_shape = NULL, 
                                     delay_scale = NULL, num.initial.cases = NULL, prop.asym = NULL, 
                                     quarantine = NULL,HH_trace = NULL,
                                     contact_data=NULL,t=NULL,p.urban=NULL,freq = NULL,
                                     pop=NULL) 
{
  # doParallel::registerDoParallel(cores=4)
  # require(foreach)
  # 
  # res <- foreach(i=1:n.sim,.combine = 'bind_rows') %dopar% {
  #   outbreak_model_ext(num.initial.cases = num.initial.cases, 
  #                      prop.ascertain = prop.ascertain, cap_max_days = cap_max_days, 
  #                      cap_cases = cap_cases, r0isolated = r0isolated, r0community = r0community, 
  #                      disp.iso = disp.iso, disp.com = disp.com, delay_shape = delay_shape, 
  #                      delay_scale = delay_scale, k = k, prop.asym = prop.asym, 
  #                      quarantine = quarantine,HH_trace = HH_trace,
  #                      contact_data=contact_data)
  # }
  
  # CD <- count.contacts.full(contact_data,t=14)
  # CD <- contact_data
  
  res <- parallel::mclapply(1:n.sim,
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
  
  
  res.week <- data.table(bind_rows(lapply(res,function(x) x$weekly_cases)))
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

scenario_sim <- function (n.sim = NULL, prop.ascertain = NULL, cap_max_days = NULL, 
                          cap_cases = NULL, r0isolated = NULL, r0community = NULL, 
                          disp.iso = NULL, disp.com = NULL, k = NULL, delay_shape = NULL, 
                          delay_scale = NULL, num.initial.cases = NULL, prop.asym = NULL, 
                          quarantine = NULL) 
{
  
  res <- parallel::mclapply(1:n.sim,
                            function(x) {
                              outbreak_model(num.initial.cases = num.initial.cases, 
                                             prop.ascertain = prop.ascertain, cap_max_days = cap_max_days, 
                                             cap_cases = cap_cases, r0isolated = r0isolated, r0community = r0community, 
                                             disp.iso = disp.iso, disp.com = disp.com, delay_shape = delay_shape, 
                                             delay_scale = delay_scale, k = k, prop.asym = prop.asym, 
                                             quarantine = quarantine)
                            }
                            ,mc.cores = parallel::detectCores())
  
  
  res <- bind_rows(res)
  res$sim <- rep(1:n.sim,each=max(res$week)+1)
  # res[, `:=`(sim, rep(1:n.sim, rep(floor(cap_max_days/7) + 
  #                                    1, n.sim))), ]
  return(res)
}

