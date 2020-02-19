##Main functions

##Bootstrapping contact data---------
contact_boot <- function(contact_data) {
  
  contact_data %>%   
    select(csid,age_class_part,day,age_class_cont,loc_stat,share_hh) -> contact_data_simple
  
  # contact_data_simple %>% head
  contact_data_simple$contact_id <- seq(1,length.out = nrow(contact_data_simple))
  contact_data_simple <- as.data.table(contact_data_simple)
  
  part.sample <- contact_data_simple[,.(csid = sample(csid,replace=TRUE))][,.(weight=.N),keyby=csid]
  contact.sample <- merge(contact_data_simple[,.(contact_id=list(c(contact_id))),by=csid],part.sample,by='csid')
  contact.sample <- contact.sample[,.(contact_id=as.double(sample(unlist(contact_id),size = weight,replace = TRUE))),by=csid]
  
  contact_data_boot <- contact_data_simple[,-c("csid")][contact.sample, on = 'contact_id']
  contact_data_boot$contact_id <- seq(1,length.out = nrow(contact_data_boot))
  
  return(contact_data_boot)
}


#New infection step from single infected id-----------
new.inf <- function(id,contact_data,r0=2,HH=FALSE,n.inf.HH=1,HH.herd.imm=FALSE) {
  
  id_new_inf_nHH <- NULL
  id_new_inf_HH <- NULL
  inds.contacted <- NULL
  
  if (!all(id%in%contact_data$csid)) stop('id is not in contact data')
  
  if (!HH) {
    n.inf.HH <- 1
    HH.herd.imm <- FALSE
  }
  if (HH & n.inf.HH<1) n.inf.HH <- 1
  
  ##Get all the contacts of id
  inds.contacted <- contact_data[csid==id,]
  
  #Total number of contacts
  c <- nrow(inds.contacted)
  #Number of new contacts (if HH member) & Choose which are already infected
  c.new <- c
  contact.id.inf <- NULL
  if (HH) {
    c.new <- nrow(inds.contacted[share_hh=="non-HH"])
    #Sample the HH memebers already infected
    if (HH.herd.imm) {
      contact.id.inf <- sample(inds.contacted[share_hh=="HH",contact_id],
                               size = min(n.inf.HH,nrow(inds.contacted[share_hh=="HH"])))   #Can skip things when all HH members are already infected...
    }
  }  
  
  c_inf <- data.table()
  
  if (r0>0) {
    
    #Which of these contacts become infected
    c_inf <- inds.contacted[sample(1:c,size = min(r0,c))] ###Is this an issue mathematically (skews R0)????
    #Remove already infected (for herd immunity??)
    if (HH.herd.imm) {
      c_inf <- c_inf[!(contact_id%in%contact.id.inf)]   ## may be a problem if same id sampled twice... fix!
    }
    
    #Get age groups of non-HH infections
    c_inf_ages_nHH <- c_inf[share_hh=='non-HH',.(age_class_cont)]
    c_inf_ages_nHH <- c_inf_ages_nHH[,.N,keyby=age_class_cont]
    
    #Get age groups of HH infections
    c_inf_ages_HH <- c_inf[share_hh=='HH',.(age_class_cont)]
    c_inf_ages_HH <- c_inf_ages_HH[,.N,keyby=age_class_cont]
    
    ##csids for newly infecteds
    # new non-HH contacts of infected individual (sampling from contacts in the same age group as infected individual)
    if (nrow(c_inf_ages_nHH)>0) {
      id_new_inf_nHH <- c_inf_ages_nHH[,contact_boot(contact_data)[age_class_part%in%age_class_cont],
                                       by=.(age_class_cont,N)
                                       ][,.(csid=list(unique(csid))),by=.(age_class_cont,N)
                                         ][,.(csid=as.double(sample(unlist(csid),size = N))),by=age_class_cont][,sort(csid)]
      
    }
    
    # new HH contacts of infected individual (sampling from contacts in the same age group as infected individual)
    if (dim(c_inf_ages_HH)[1]>0) {
      id_new_inf_HH <- c_inf_ages_HH[,contact_boot(contact_data)[age_class_part%in%age_class_cont],
                                     by=.(age_class_cont,N)
                                     ][,.(csid=list(unique(csid))),by=.(age_class_cont,N)
                                       ][,.(csid=as.double(sample(unlist(csid),size = N))),by=age_class_cont][,sort(csid)]
    }
  }
  
  return(list(total.contacts=c,
              total.new.contacts= c.new,
              inf.contacts=length(c(id_new_inf_nHH,id_new_inf_HH)),
              id.inf.HH=id_new_inf_HH,
              id.inf.nHH=id_new_inf_nHH,
              n.inf.HH=if (length(id_new_inf_HH)>0) {
                rep(n.inf.HH,times=length(id_new_inf_HH))+length(id_new_inf_HH)
                } else NULL
  )
  )
  
}


##Simulate a single outbreak simulation-----------
bp.inf <- function(init_inf=1,n=10,contact_data,HH.herd.imm=TRUE) {
  
  ##To track
  total.contacts <- rep(NA,n)
  total.new.contacts <- rep(NA,n)
  inf.contacts <- rep(NA,n)
  max.HH.inf <- rep(NA,n)
  mean.HH.inf <- rep(NA,n)
  r0s <- rep(NA,n)
  r0_eff <- rep(NA,n)
  n.inf.HH <- 0
  
  ##Set up initial random infected seed ids
  ids.nHH <- sample(unique(contact_data$csid),init_inf)   
  ids.HH <- NULL
  
  ##Get some r0 values (so you can skip loop if all==0)
  r0.nHH <- rnbinom(length(ids.nHH),size = .16, mu = 2.5)
  r0.HH <- NULL
  
  if (!all(c(r0.nHH,r0.HH)==0)) {
    
    for (i in 1:n) {
      
      # X.HH <- NULL
      # X.nHH <- NULL
      
      if (i>1) {
        ## Get some r0 values
        r0.nHH <- rnbinom(length(ids.nHH),size = .16, mu = 2.5)
        r0.HH <- rnbinom(length(ids.HH),size = .16, mu = 2.5)
      }
      r0s[i] <- mean(c(r0.nHH,r0.HH))
      
      ##Get contact data sample, make sure all ids are present
      CD <- NULL
      while (!all(c(ids.nHH,ids.HH)%in%CD$csid)) {
        CD <- contact_boot(contact_data)
      }
      
      ##Count the contacts
      total.contacts[i] <- nrow(CD[csid%in%c(ids.nHH,ids.HH)])
      total.new.contacts[i] <- nrow(CD[csid%in%c(ids.nHH)]) +
        nrow(CD[csid%in%c(ids.HH)][share_hh=="non-HH"])
      
      if (all(c(r0.nHH,r0.HH)==0)) break
      
      
      ##Only need to simulate ids with R0>0
      ids.nHH <- ids.nHH[r0.nHH>0]
      ids.HH <- ids.HH[r0.HH>0]
      r0.nHH <- r0.nHH[r0.nHH>0]
      n.inf.HH <- n.inf.HH[r0.HH>0]
      r0.HH <- r0.HH[r0.HH>0]
      
      ##Infection step for HH contacts
      if (length(ids.HH) > 0) {
        X.HH <- foreach(i=1:length(ids.HH)) %dopar% {
          new.inf(id = ids.HH[i],contact_data = CD,r0 = r0.HH[i],HH = TRUE,n.inf.HH = n.inf.HH[i],HH.herd.imm=HH.herd.imm)
        }
      }
      
      ##Infection step for non-HH contacts
      if (length(ids.nHH) > 0) {
        X.nHH <- foreach(i=1:length(ids.nHH)) %dopar% {
          new.inf(id = ids.nHH[i],contact_data = CD,r0 = r0.nHH[i],HH = FALSE)
        }
      }
      
      ##New infections
      inf.contacts[i] <- ifelse(is.null(X.nHH),0,sum(sapply(X.nHH,function(x) sum(x$inf.contacts)))) + 
        ifelse(is.null(X.HH),0,sum(sapply(X.HH,function(x) sum(x$inf.contacts))))
      
      ##Effective R
      ##Check this, not quite right, maybe put it inside the new infection step!!!!
      r0_eff[i] <- mean(c(ifelse(is.null(X.nHH),NA,mean(sapply(X.nHH,function(x) c(x$inf.contacts)))),
                          ifelse(is.null(X.HH),NA,mean(sapply(X.HH,function(x) c(x$inf.contacts))))),na.rm = TRUE)
      
      ##Get the largest infected HH
      max.HH.inf[i] <- max(c(unlist(sapply(X.HH,function(x) c(x$n.inf.HH))),
                             unlist(sapply(X.nHH,function(x) c(x$n.inf.HH)))))
      mean.HH.inf[i] <- mean(c(unlist(sapply(X.HH,function(x) c(x$n.inf.HH))),
                               unlist(sapply(X.nHH,function(x) c(x$n.inf.HH)))))
      
      #Newly infected HH ids
      ids.HH <- c(unlist(sapply(X.HH,function(x) c(x$id.inf.HH))),
                  unlist(sapply(X.nHH,function(x) c(x$id.inf.HH))))
      
      #Newly infected non-HH ids
      ids.nHH <- c(unlist(sapply(X.HH,function(x) c(x$id.inf.nHH))),
                   unlist(sapply(X.nHH,function(x) c(x$id.inf.nHH))))
      
      ##For each HH id, number of people infected in that HH
      n.inf.HH <- c(unlist(sapply(X.HH,function(x) c(x$n.inf.HH))),
                    unlist(sapply(X.nHH,function(x) c(x$n.inf.HH))))
      
      if (length(c(ids.nHH,ids.HH))<1) break
      
    }
  }
  return(data.frame(gen=1:n,
                    contacts=total.contacts,
                    new.contacts=total.new.contacts,
                    infections=inf.contacts,
                    r0_mean=r0s,
                    r_eff=r0_eff,
                    max.HH.inf,
                    mean.HH.inf)) ##Track "avoided" HH infections as well?
  
}