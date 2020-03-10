##Contact data functions
require(tidyverse)
require(data.table)

##Get and tidy contact data
contact_data <- as.data.frame(read_csv('data/journal.pone.0104786.s007.CSV'))
age_groups <- c('<1','1-5','6-14','15-19','20-49','>=50')
contact_data %<>% 
  mutate_all(as.factor) %>%
  mutate(contact=1) 
levels(contact_data$loc_stat) <- c("rural","urban")
levels(contact_data$season) <- c("dry","wet")
##Fix the one missing share_HH to nHH (based on similar part and contact age, majority were non-HH)
levels(contact_data$share_hh) <- c("Missing","non-HH","HH")
contact_data$share_hh[contact_data$share_hh=="Missing"] <- "non-HH"

levels(contact_data$age_class_cont) <- c(5,1,6,2,3,4)
contact_data$age_class_cont <- factor(contact_data$age_class_cont,levels = sort(levels(contact_data$age_class_cont)))
levels(contact_data$age_class_cont) <- age_groups
levels(contact_data$age_class_part) <- age_groups
contact_data$csid <- as.integer(as.character(contact_data$csid))

contact_data <- data.table(contact_data)
##Do some fixing of missing/inconsistent data
#Fix missing howoften that are in HH to daily
contact_data[share_hh=="HH" & how_often=='Missing',how_often := 'Daily']
##Fix how_often missing and evermet==no to never
contact_data[ever_met=="No" & how_often=='Missing',how_often := 'Never']
#Fix missing how_often and ever_met==Missing to sample based on rest of data 
missing.contacts <- contact_data[how_often=='Missing']

contact_data[how_often=='Missing',"how_often"] <- sapply(1:nrow(missing.contacts),
                                                         function(i) sample(contact_data[how_often!='Missing'][
                                                           missing.contacts[i],how_often,
                                                           on=.(age_class_part,age_class_cont,loc_stat,share_hh)],1))

##Bootstrapping and sampling from original contact data set to form a study population size of n---------
contact_boot <- function(contact_data,n,p.urban) {
  
  contact_data %>%   
    select(csid,age_class_part,day,age_class_cont,loc_stat,share_hh,how_often,ever_met) -> contact_data_simple
  
  contact_data_simple$contact_id <- seq(1,length.out = nrow(contact_data_simple))
  contact_data_simple <- as.data.table(contact_data_simple)
  # contacts <- contact_data_simple[,.(contacts=.N),by=csid]
  # 
  # participants.urban <- unique(contact_data_simple[loc_stat=="urban",csid])
  # participants.rural <- unique(contact_data_simple[loc_stat=="rural",csid])
  # 
  # ## Sample participants according to urban/rural prop
  # participants <- data.table(csid.new=1:n,csid=c(sample(participants.urban,round(n*p.urban),replace = TRUE),
  #                                                sample(participants.urban,n-round(n*p.urban),replace = TRUE)))
  # 
  
  # contact.sample <- merge(contact_data_simple[,.(contact_id=list(c(contact_id))),by=csid],participants,by='csid')
  # merge(contact.sample[,.(contact_id=unlist(contact_id)),by=csid.new],contact_data_simple,by="contact_id")[order(csid.new)]
  
  
  part.sample <- contact_data_simple[,.(csid = sample(csid,replace=TRUE))][,.(weight=.N),keyby=csid]
  contact.sample <- merge(contact_data_simple[,.(contact_id=list(c(contact_id))),by=csid],part.sample,by='csid')
  contact.sample <- contact.sample[,.(contact_id=as.double(sample(unlist(contact_id),size = weight,replace = TRUE))),by=csid]
  
  contact_data_boot <- contact_data_simple[,-c("csid")][contact.sample, on = 'contact_id']
  contact_data_boot$contact_id <- seq(1,length.out = nrow(contact_data_boot))
  
  # ##Do some fixing of missing/inconsistent data
  # #Fix missing howoften that are in HH to daily
  # contact_data_boot[share_hh=="HH" & how_often=='Missing',how_often := 'Daily']
  # ##Fix how_often missing and evermet==no to never
  # contact_data_boot[ever_met=="No" & how_often=='Missing',how_often := 'Never']
  # #Fix missing how_often and ever_met==Missing to sample based on rest of data 
  # missing.contacts <- contact_data_boot[how_often=='Missing']
  # 
  # contact_data_boot[how_often=='Missing',"how_often"] <- sapply(1:nrow(missing.contacts),
  #                                                               function(i) sample(contact_data_boot[how_often!='Missing'][
  #                                                                 missing.contacts[i],how_often,
  #                                                                 on=.(age_class_part,age_class_cont,loc_stat,share_hh)],1))
  # 
  
  # missing.contacts[.(how_often=sample(contact_data_boot[,how_often]))]
  # merge(missing.contacts,contact_data_boot[how_often!='Missing'],
  #       by=c("age_class_part","age_class_cont","loc_stat","share_hh"))
  # 
  return(contact_data_boot)
}

##Bootstrapping and sampling from original contact data set to form a study population size of n---------
contact_boot <- function(contact_data,n,p.urban) {
  
  contact_data %>%   
    select(csid,age_class_part,day,age_class_cont,loc_stat,share_hh,how_often,ever_met) -> contact_data_simple
  
  contact_data_simple$contact_id <- seq(1,length.out = nrow(contact_data_simple))
  contact_data_simple <- as.data.table(contact_data_simple)
  contacts <- contact_data_simple[,.(contacts=.N),by=csid]
  
  participants.urban <- unique(contact_data_simple[loc_stat=="urban",csid])
  participants.rural <- unique(contact_data_simple[loc_stat=="rural",csid])
  
  ## Sample participants according to urban/rural prop
  participants <- data.table(csid.new=1:n,csid=c(sample(participants.urban,round(n*p.urban),replace = TRUE),
                    sample(participants.urban,n-round(n*p.urban),replace = TRUE)))
  
  
  contact.sample <- merge(contact_data_simple[,.(contact_id=list(c(contact_id))),by=csid],participants,by='csid')
  contact_data_boot <- merge(contact.sample[,.(contact_id=unlist(contact_id)),by=csid.new],contact_data_simple,by="contact_id")[order(csid.new)]
  contact_data_boot %<>% select(-csid) %>% mutate(csid=csid.new) %>% select(-csid.new) %>% data.table()
  
  # 
  # part.sample <- contact_data_simple[,.(csid = sample(csid,replace=TRUE))][,.(weight=.N),keyby=csid]
  # contact.sample <- merge(contact_data_simple[,.(contact_id=list(c(contact_id))),by=csid],part.sample,by='csid')
  # contact.sample <- contact.sample[,.(contact_id=as.double(sample(unlist(contact_id),size = weight,replace = TRUE))),by=csid]
  # 
  # contact_data_boot <- contact_data_simple[,-c("csid")][contact.sample, on = 'contact_id']
  contact_data_boot$contact_id <- seq(1,length.out = nrow(contact_data_boot))
  
  # ##Do some fixing of missing/inconsistent data
  # #Fix missing howoften that are in HH to daily
  # contact_data_boot[share_hh=="HH" & how_often=='Missing',how_often := 'Daily']
  # ##Fix how_often missing and evermet==no to never
  # contact_data_boot[ever_met=="No" & how_often=='Missing',how_often := 'Never']
  # #Fix missing how_often and ever_met==Missing to sample based on rest of data 
  # missing.contacts <- contact_data_boot[how_often=='Missing']
  # 
  # contact_data_boot[how_often=='Missing',"how_often"] <- sapply(1:nrow(missing.contacts),
  #                                                               function(i) sample(contact_data_boot[how_often!='Missing'][
  #                                                                 missing.contacts[i],how_often,
  #                                                                 on=.(age_class_part,age_class_cont,loc_stat,share_hh)],1))
  # 
  
  # missing.contacts[.(how_often=sample(contact_data_boot[,how_often]))]
  # merge(missing.contacts,contact_data_boot[how_often!='Missing'],
  #       by=c("age_class_part","age_class_cont","loc_stat","share_hh"))
  # 
  return(contact_data_boot)
}


##Given a single contact with probability p of being repeated each dat (based on frequency)
# how many new such contacts will be made in t days (with some randomness)
new.contacts <- function(t,p) {
  n.c <- 1
  P <- 1-p
  for (i in 2:t) {
    if (runif(1) <= P) {
      P <- P*(1-p)
      n.c <- n.c+1
    }
  }
  return(n.c)
}


##Count number of contacts for a given period weighting by "how_often" and "ever_met" measure
##output: list of contacts and a weight for becoming infected across time t
count.contacts <- function(id,CD,t=14) {
  
  ## Add weights to count how many contacts are repeated during t
  weights <- c(1,0,1.5/7,1.5/30,.5/30,0)
  CD$weights <- weights[as.numeric(CD$how_often)]
  
  id.contacts <- CD[csid==id]
  new.id.contacts <- id.contacts
  
  if (nrow(id.contacts)!=0) {
    
    ## Count how many new contacts of each type were made during t days
    id.contacts <- id.contacts[,':='(new.count = new.contacts(t,weights)),by=contact_id]
    
    ##Extend by new contacts and give them new ids
    new.id.contacts <- id.contacts[,.(new.contact_id=rep(contact_id,new.count)),by=c(colnames(id.contacts)[1:9])][,-"new.contact_id"]
    new.id.contacts$contact_id <- 1:nrow(new.id.contacts)  
    
  }
  
  #Add final weights for how often each contact was made during t
  t.weights <- c(1,1/t,1.5/7,1.5/30,.5/30,1/t)*t
  new.id.contacts$weight <- t.weights[as.numeric(new.id.contacts$how_often)]
  
  return(new.id.contacts)
}

##Given a contact data set, bootstrap and output all the contacts across time t
count.contacts.full <- function(CD,t=14) {
  
  CD_boot <- contact_boot(CD,n = 10000,p.urban = .3)
  
  # participants <- unique(CD_boot$csid)
  
  participants <- CD_boot[,.(age=unique(age_class_part)),by=csid]
  
  # system.time(t_contacts <- bind_rows(participants %>% purrr::map(~count.contacts(.x,CD_boot,t))))
  # 
  # system.time(t_contacts <- bind_rows(lapply(participants,function(x) count.contacts(x,CD_boot,t))))
  
  t_contacts <- bind_rows(parallel::mclapply(participants$csid,function(x) count.contacts(x,CD_boot,t),mc.cores = 4))[,c(9,1:8,10)]
  
  t_contacts$contact=paste(t_contacts$csid,t_contacts$contact_id,sep='.')
  
  
  # participants_demo <- dcast(t_contacts[,.(count=.N),by=.(csid,age_class_part,share_hh)],
  #                            csid+age_class_part~share_hh)
  # colnames(participants_demo) <- c("id","age","no.nHH.con","no.con")
  # participants_demo[is.na(participants_demo[,no.nHH.con]),`:=`(no.nHH.con=0)]
  # participants_demo[is.na(participants_demo[,no.con]),`:=`(no.con=0)]
  # participants_demo[,`:=`(no.con=no.con+no.nHH.con)]
  
  
  # participants_demo<-data.frame(id=participants,
  #                               age=t_contacts$age_class_part[match(participants,t_contacts$csid)],
  #                               no.con= sapply(participants, function(p)sum(t_contacts$csid==p)),
  #                               no.nHH.con=sapply(participants, function(p)sum(t_contacts$csid==p & t_contacts$share_hh=="non-HH")))
  
  # return(list(participants,t_contacts,participants_demo))
  return(list(participants,t_contacts))
  
}
