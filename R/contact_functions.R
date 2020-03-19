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


##Bootstrapping and sampling from original contact data set by proportion urban vs rural by frequencies---------
contact_sampling <- function(contact_data,p.urban,n=1000,freq="Daily") {
  
  contact_data %>% subset(share_hh=="HH" | how_often%in%freq) %>% 
    select(csid,age_class_part,day,age_class_cont,loc_stat,share_hh,how_often,ever_met) -> contact_data_simple
  
  contact_data_simple$contact_id <- seq(1,length.out = nrow(contact_data_simple))
  contact_data_simple <- as.data.table(contact_data_simple)
  contacts <- contact_data_simple[,.(contacts=.N),by=csid]
  
  participants.urban <- unique(contact_data_simple[loc_stat=="urban",csid])
  participants.rural <- unique(contact_data_simple[loc_stat=="rural",csid])
  
  ## Sample participants according to urban/rural prop
  participants <- data.table(csid.new=1:n,csid=c(sample(participants.urban,round(n*p.urban),replace = TRUE),
                                                 sample(participants.rural,n-round(n*p.urban),replace = TRUE)))
  
  
  contact.sample <- merge(contact_data_simple[,.(contact_id=list(c(contact_id))),by=csid],participants,by='csid')
  contact_data_boot <- merge(contact.sample[,.(contact_id=unlist(contact_id)),by=csid.new],contact_data_simple,by="contact_id")[order(csid.new)]
  contact_data_boot %<>% select(-csid) %>% mutate(csid=csid.new) %>% select(-csid.new) %>% data.table()
  
  contact_data_boot$contact_id <- seq(1,length.out = nrow(contact_data_boot))
  
  return(contact_data_boot)
}

##Bootstrapping and sampling from original contact data set to form a study population size of n---------
contact_boot <- function(contact_data,n,p.urban,freq) {
  
  contact_data %>% subset(share_hh=="HH" | how_often%in%freq) %>% 
    select(csid,age_class_part,day,age_class_cont,loc_stat,share_hh,how_often,ever_met) -> contact_data_simple
  
  contact_data_simple$contact_id <- seq(1,length.out = nrow(contact_data_simple))
  contact_data_simple <- as.data.table(contact_data_simple)
  contacts <- contact_data_simple[,.(contacts=.N),by=csid]
  
  participants.urban <- unique(contact_data_simple[loc_stat=="urban",csid])
  participants.rural <- unique(contact_data_simple[loc_stat=="rural",csid])
  
  ## Sample participants according to urban/rural prop
  participants <- data.table(csid.new=1:n,csid=c(sample(participants.urban,round(n*p.urban),replace = TRUE),
                                                 sample(participants.rural,n-round(n*p.urban),replace = TRUE)))
  
  
  contact.sample <- merge(contact_data_simple[,.(contact_id=list(c(contact_id))),by=csid],participants,by='csid')
  contact_data_boot <- merge(contact.sample[,.(contact_id=unlist(contact_id)),by=csid.new],contact_data_simple,by="contact_id")[order(csid.new)]
  contact_data_boot %<>% select(-csid) %>% mutate(csid=csid.new) %>% select(-csid.new) %>% data.table()
  
  
  contact_data_boot$contact_id <- seq(1,length.out = nrow(contact_data_boot))
  
  return(contact_data_boot)
}

##Given a single contact with probability p of being repeated each dat (based on frequency)
# how many new such contacts will be made in t days (with some randomness)
new.contacts <- function(t,p) {
  n.c <- t
  if (p>0 & p<1) {
    n.c <- 1
    P <- 1-p
    for (i in 2:t) {
      if (runif(1) <= P) {
        P <- P*(1-p)
        n.c <- n.c+1
      }
    }
  }
  return(n.c)
}

##Given a contact data set, sample n individuals from it and output all the contacts across time t
#Then take a sampled subset to roughly match the desired population size
count.contacts <- function(CD,t=14,p.urban = .3,freq="Daily",pop=200000) {
  
  #Bootstrap pop/20 participants to get a large enough group, since on average people have 20 contacts
  
  # p.rep <- c(Daily=1,Often=1.5/7,Regularly=1.5/30,Rarely=.5/30,Never=0)
  # mean.new <- 1+(1-p.rep)*(t-1)
  # mean_ct <- cbind(CD[,.N,by=.(how_often,csid)][,.(ct=mean(N)),by=how_often],mean.new)[,(ct=sum(ct+mean.new))]
  # 
  ##Get sample from contct data with urban prop. n = desired population size divded by the mean number of contacts per participant
  CD_boot <- contact_sampling(CD,p.urban,freq = freq,n=round(pop/mean(CD[how_often%in%freq,.N,by=csid]$N)))
  
  #Split into daily and non-daily to save us some work
  CD_boot_daily <- CD_boot[how_often=="Daily",]
  CD_boot_non_daily <- CD_boot[how_often!="Daily",]
  
  ## Add weights to count how many contacts are repeated during t (for non-daily)
  weights <- c(1,0,1.5/7,1.5/30,.5/30,0)
  CD_boot_non_daily$weights <- weights[as.numeric(CD_boot_non_daily$how_often)]
  
  ## Count how many new contacts of each type were made during t days (for non-daily)
  CD_boot_non_daily[,':='(new.count = new.contacts(t,weights)),by=contact_id]
  
  ##Extend by new contacts and give them new ids (for non-daily)
  CD_boot_non_daily <- CD_boot_non_daily[how_often!="Daily",][,.(contact_id=rep(contact_id,new.count)),
                                                              by=c(colnames(CD_boot_non_daily)[1:9])][,-"contact_id"]
  
  ##Combine daily and new non-daily
  CD_new <- bind_rows(CD_boot_daily,CD_boot_non_daily)[order(csid),]
  ##Sample to roughly get desired population size
  parts <- sample(unique(CD_new$csid),round(pop/mean(CD_new[how_often%in%freq,.N,by=csid]$N)))
  parts <- data.table(csid=parts,age=CD_new[csid%in%parts,.N,by=.(csid,age_class_part)][,age_class_part])
  CD_new <- CD_new[csid%in%parts$csid]
  
  ##Add contacts ids
  CD_new$contact <- 1:nrow(CD_new)  
  
  #Add final weights for how often each contact was made during t
  t.weights <- c(1,1/t,1.5/7,1.5/30,.5/30,1/t)*t
  CD_new$weight <- t.weights[as.numeric(CD_new$how_often)]
  
  CD_new$contact=paste(CD_new$csid,CD_new$contact,sep='.')
  
  return(list(parts,CD_new))
}

##Given a contact data set, sample n individuals from it and output all the contacts across time t
#Then take a sampled subset to roughly match the desired population size
count.contacts2 <- function(CD,t=14,p.urban = .3,freq="Daily",pop=200000) {
  
  #Bootstrap pop/20 participants to get a large enough group, since on average people have 20 contacts
  
  # p.rep <- c(Daily=1,Often=1.5/7,Regularly=1.5/30,Rarely=.5/30,Never=0)
  # mean.new <- 1+(1-p.rep)*(t-1)
  # mean_ct <- cbind(CD[,.N,by=.(how_often,csid)][,.(ct=mean(N)),by=how_often],mean.new)[,(ct=sum(ct+mean.new))]
  # 
  ##Get sample from contct data with urban prop. n = desired population size divded by the mean number of contacts per participant
  CD_boot <- contact_sampling(CD,p.urban,freq = freq,n=round(pop/mean(CD[how_often%in%freq,.N,by=csid]$N)))
  
  #Split into daily and non-daily to save us some work
  CD_boot_daily <- CD_boot[how_often=="Daily",]
  CD_boot_non_daily <- CD_boot[how_often!="Daily",]
  
  ## Add weights to count how many contacts are repeated during t (for non-daily)
  weights <- c(1,0,1.5/7,1.5/30,.5/30,0)
  CD_boot_non_daily$weights <- weights[as.numeric(CD_boot_non_daily$how_often)]
  
  ## Count how many new contacts of each type were made during t days (for non-daily)
  CD_boot_non_daily[,':='(new.count = new.contacts(t,weights)),by=contact_id]
  
  ##Extend by new contacts and give them new ids (for non-daily)
  CD_boot_non_daily <- CD_boot_non_daily[how_often!="Daily",][,.(contact_id=rep(contact_id,new.count)),
                                                              by=c(colnames(CD_boot_non_daily)[1:9])][,-"contact_id"]
  
  ##Combine daily and new non-daily
  CD_new <- bind_rows(CD_boot_daily,CD_boot_non_daily)[order(csid),]
  ##Sample to roughly get desired population size
  parts <- sample(unique(CD_new$csid),round(pop/mean(CD_new[how_often%in%freq,.N,by=csid]$N)))
  parts <- data.table(csid=parts,age=CD_new[csid%in%parts,.N,by=.(csid,age_class_part)][,age_class_part])
  CD_new <- CD_new[csid%in%parts$csid]
  
  ##Add contacts ids
  CD_new$contact <- 1:nrow(CD_new)  
  
  #Add final weights for how often each contact was made during t
  t.weights <- c(1,1/t,1.5/7,1.5/30,.5/30,1/t)*t
  CD_new$weight <- t.weights[as.numeric(CD_new$how_often)]
  
  CD_new$contact=paste(CD_new$csid,CD_new$contact,sep='.')
  
  
  ### Now create full set of contacts for participants + contacts to form a close network
  parts
  parts.age <- parts[,.(csid=list(csid)),by=age]
  parts_full <- bind_rows(parts[,`:=`(csid=as.character(csid))],CD_new[,.(csid=contact,age=age_class_cont)])
  
  # x <- CD_new[,.(hh.contact_id=list(contact[share_hh=='HH']),
  #                contact_id=list(contact[share_hh!='HH']),
  #                hh.weight=list(weight[share_hh=='HH']),
  #                nhh.weight=list(weight[share_hh!='HH'])),by=csid]
  
  # x <- CD_new[,.(hh.contact_id=list(contact[share_hh=='HH']),
  #                nhh.contact_id=list(contact[share_hh!='HH']),
  #                contact_id=list(contact),
  #                weight=list(weight)),by=csid]
  
  x <- CD_new[,.(hh.contact_id=list(contact[share_hh=='HH']),
                 hh.weight=list(weight[share_hh=='HH']),
                 nhh.contact_id=list(contact[share_hh!='HH']),
                 nhh.weight=list(weight[share_hh!='HH'])),
              by=csid][,`:=`(csid=csid)]
  
  ## non-HH contacts
  y.nhh <- CD_new[share_hh!="HH",.(csid=contact,age=age_class_cont)]
  # y.nhh[,parts.age[age==age,sample(unlist(csid),1)],by=csid]
  y <- merge(parts.age,y.nhh[,.(.N,csid.orig=list(csid)),by=age])
  y.nhh <- data.table(bind_cols(csid=y.nhh[order(age)]$csid,
                                y[,.(csid.match=as.character(sample(unlist(csid),N,replace = TRUE))),
                                  by=age]))
  y.nhh <- merge(y.nhh,x[,`:=`(csid.match=as.character(csid))][,-"csid"])[order(csid),-"csid.match"]
  
  ## HH contacts
  y.hh <- CD_new[share_hh=="HH",.(csid=contact,age=age_class_cont)]
  y <- merge(parts.age,y.hh[,.(.N,csid.orig=list(csid)),by=age])
  y.hh <- data.table(bind_cols(csid=y.hh[order(age)]$csid,
                                y[,.(csid.match=as.character(sample(unlist(csid),N,replace = TRUE))),
                                  by=age]))
  #Get non-HH contacts by taking all the contacts of a random individual of same age
  y.hh <- merge(y.hh,x[,`:=`(csid.match=as.character(csid))][,-c("csid","hh.contact_id","hh.weight")])[order(csid),-"csid.match"]
  
  #Get their HH contacts by using the same ones as their original HH contactee (and swapping the csid with the contactees csid)
  y.hh <- merge(y.hh[,`:=`(csid.match=floor(as.numeric(csid)))],
                x[,`:=`(csid.match=csid)][,c("csid.match","hh.contact_id","hh.weight")])[,-"csid.match"]
  y.hh2 <- y.hh[,.(hh.contact_id=unlist(hh.contact_id)),by=csid][,`:=`(hh.contact_id=ifelse(csid==hh.contact_id,
                                                                                  floor(as.numeric(hh.contact_id)),
                                                                                  hh.contact_id))]
  y.hh$hh.contact_id <- y.hh2[,.(hh.contact_id=list(hh.contact_id)),by=csid][,-"csid"]
  
  ## Put them together
  y <- bind_rows(y.hh,y.nhh)[order(csid)]
  
  ##now put the full contact dataset set
  
  parts_full <- bind_rows(parts[,`:=`(csid=as.character(csid))],CD_new[,.(csid=contact,age=age_class_cont)])
  CD_new_full <- bind_rows(x[,`:=`(csid=as.numeric(csid))][,-"csid.match"],
                           y[,`:=`(csid=as.numeric(csid))][,-"age"])[
                             order(csid)]
  
  return(list(parts_full,CD_new_full))
}

# cd <- count.contacts2(contact_data,t=14,p.urban = .3,freq=c("Daily","Often","Regularly","Rarely","Never"),pop=200000)


##Given a contact data set, sample n individuals from it and output all the contacts across time t
#Then take a sampled subset to roughly match the desired population size
count.contacts3 <- function(CD,t=14,p.urban = .3,freq="Daily",pop=200000) {

  ##Get sample from contct data with urban prop. n = desired population size divded by the mean number of contacts per participant
  # CD_boot <- contact_sampling(CD,p.urban,freq = freq,n=round(pop/mean(CD[how_often%in%freq,.N,by=csid]$N)))
  CD_boot <- contact_sampling(CD,p.urban,freq = freq,n=round(pop/mean(CD[share_hh=='HH',.N,by=csid]$N)))
  
  #Split into daily and non-daily to save us some work
  CD_boot_daily <- CD_boot[how_often=="Daily",]
  CD_boot_non_daily <- CD_boot[how_often!="Daily",]
  
  ## Add weights to count how many contacts are repeated during t (for non-daily)
  weights <- c(1,0,1.5/7,1.5/30,.5/30,0)
  CD_boot_non_daily$weights <- weights[as.numeric(CD_boot_non_daily$how_often)]
  
  ## Count how many new contacts of each type were made during t days (for non-daily)
  CD_boot_non_daily[,':='(new.count = new.contacts(t,weights)),by=contact_id]
  
  ##Extend by new contacts and give them new ids (for non-daily)
  CD_boot_non_daily <- CD_boot_non_daily[how_often!="Daily",][,.(contact_id=rep(contact_id,new.count)),
                                                              by=c(colnames(CD_boot_non_daily)[1:9])][,-"contact_id"]
  
  ##Combine daily and new non-daily
  CD_new <- bind_rows(CD_boot_daily,CD_boot_non_daily)[order(csid),]
  
  ##Sample to roughly get desired population size
  # parts <- sample(unique(CD_new$csid),round(pop/mean(CD_new[how_often%in%freq,.N,by=csid]$N)))
  parts <- unique(CD_new$csid)
  parts <- data.table(csid=as.character(parts),
                      age=CD_new[,.N,by=.(csid,age_class_part)][,age_class_part],
                      CD_new[,.(hh.con=sum(share_hh=="HH"),nhh.con=sum(share_hh!="HH")),by=csid][,-"csid"])
  parts.age <- parts[,.(csid=list(csid)),by=age]
  # CD_new <- CD_new[csid%in%parts$csid]
  CD_new$contact_id <- as.character(CD_new$contact_id)
  CD_new$contact_id <- NA
  CD_new$csid <- as.character(CD_new$csid)
  
  ##Add contact_ids to HH members: These can be new individuals
  CD_new[share_hh=="HH"]$contact_id <- 1:nrow(CD_new[share_hh=="HH"])
  CD_new[share_hh=="HH",`:=`(contact_id=as.character(paste(csid,contact_id,sep = '.')))]
  
  ##Add contact ids to non-HH members so they can be matched by age later
  CD_new[share_hh!="HH"]$contact_id <- 1:nrow(CD_new[share_hh!="HH"])
  CD_new[share_hh!="HH",`:=`(contact_id=as.character(paste("nhh",csid,contact_id,sep = '.')))]
  
  #Add final weights for how often each contact was made during t
  t.weights <- c(1,1/t,1.5/7,1.5/30,.5/30,1/t)*t
  CD_new$weight <- t.weights[as.numeric(CD_new$how_often)]
  
  CD_new <- CD_new[,.(csid,age_class_part,share_hh,age_class_cont,contact_id,weight)]
  
  
  x <- CD_new[,.(hh.contact_id=list(contact_id[share_hh=='HH']),
                 hh.weight=list(weight[share_hh=='HH']),
                 nhh.contact_id=list(contact_id[share_hh!='HH']),
                 nhh.weight=list(weight[share_hh!='HH'])),
              by=csid][,`:=`(csid=csid)]
  
  # x <- CD_new[,.(contact_id=list(contact_id),
  #                weight=list(weight[share_hh=='HH'])),
  #             by=csid][,`:=`(csid=csid)]
  # 
  ## HH contacts
  y.hh <- CD_new[share_hh=="HH",.(csid=contact_id,age=age_class_cont)]
  y <- merge(parts.age,y.hh[,.(.N,csid.orig=list(csid)),by=age])
  y.hh <- data.table(bind_cols(csid=y.hh[order(age)]$csid,
                               y[,.(csid.match=as.character(sample(unlist(csid),N,replace = TRUE))),
                                 by=age]))
  #Get non-HH contacts by taking all the contacts of a random individual of same age
  y.hh <- merge(y.hh,x[,`:=`(csid.match=as.character(csid))][,-c("csid","hh.contact_id","hh.weight")])[order(csid),-"csid.match"]
  
  #Get their HH contacts by using the same ones as their original HH contactee (and swapping the csid with the contactees csid)
  y.hh <- merge(y.hh[,`:=`(csid.match=as.character(floor(as.numeric(csid))))],
                x[,`:=`(csid.match=csid)][,c("csid.match","hh.contact_id","hh.weight")])[,-"csid.match"]
  y.hh2 <- y.hh[,.(hh.contact_id=unlist(hh.contact_id)),by=csid][,`:=`(hh.contact_id=ifelse(csid==hh.contact_id,
                                                                                            floor(as.numeric(hh.contact_id)),
                                                                                            hh.contact_id))]
  y.hh$hh.contact_id <- y.hh2[,.(hh.contact_id=list(hh.contact_id)),by=csid][,-"csid"]
  
  
  # parts_full <- bind_rows(parts2,parts3)
  CD_new_full <- bind_rows(x[,`:=`(csid=as.character(csid))][,-"csid.match"],
                           y.hh[,`:=`(csid=as.character(csid))][,-"age"])[
                             order(csid)]
  
  CD_new_full <- CD_new_full[,.(contact_id=list(c(unlist(hh.contact_id),unlist(nhh.contact_id))),
                 weight=list(c(unlist(hh.weight),unlist(nhh.weight)))),by=csid]
  
  parts2 <- bind_rows(parts,
                      y.hh[,.(hh.con=length(unlist(hh.contact_id)),
                              nhh.con=length(unlist(nhh.contact_id))),by=.(csid,age)])
  parts3 <- CD_new[share_hh!="HH",.(csid=contact_id,age=age_class_cont)]
  
  return(list(parts2,parts3,CD_new_full))
}

# cd <- count.contacts3(contact_data,t=7,p.urban = .3,freq=c("Daily","Often","Regularly","Rarely","Never"),pop=200000)








