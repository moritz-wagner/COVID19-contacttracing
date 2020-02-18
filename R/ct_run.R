library(readr)
library(tidyverse)
library(ggplot2)
library(data.table)

contact_data <- as.data.frame(read_csv('data/journal.pone.0104786.s007.CSV'))
age_groups <- c('<1','1-5','6-14','15-19','20-49','>=50')
contact_data %<>% 
  mutate_all(as.factor) %>%
  mutate(contact=1) 
levels(contact_data$loc_stat) <- c("rural","urban")
levels(contact_data$season) <- c("dry","wet")
levels(contact_data$share_hh) <- c("Missing","non-HH","HH")

levels(contact_data$age_class_cont) <- c(5,1,6,2,3,4)
contact_data$age_class_cont <- factor(contact_data$age_class_cont,levels = sort(levels(contact_data$age_class_cont)))
levels(contact_data$age_class_cont) <- age_groups
levels(contact_data$age_class_part) <- age_groups
contact_data$csid <- as.integer(as.character(contact_data$csid))

##Compare bootstrap to original data
#Hist of contacts
ggplot()+
  geom_histogram(data=contact_data %>% 
                   group_by(csid,day) %>% 
                   summarise(contacts=sum(contact)),
                 aes(x=contacts),alpha=.3,
                 binwidth = 1)+
  geom_vline(data=contact_data %>% 
               group_by(csid,day) %>% 
               summarise(contacts=sum(contact)),
             aes(xintercept=mean(contacts)),color="red")+
  geom_histogram(data=contact_boot(contact_data) %>% 
                   mutate(contact=1) %>% 
                   group_by(csid,day) %>% 
                   summarise(contacts=sum(contact)),
                 aes(x=contacts),fill="green",alpha=.3,
                 binwidth = 1)+
  geom_vline(data=contact_boot(contact_data) %>% 
               mutate(contact=1) %>% 
               group_by(csid,day) %>% 
               summarise(contacts=sum(contact)),
             aes(xintercept=mean(contacts)),color="blue")


##New Infection check
id <- sample(contact_data$csid,1)
cd <- contact_boot(contact_data)
new.inf(sample(contact_data$csid,1),cd,r0 = 30,HH = TRUE,n.inf.HH=3)

library(doParallel)
pkgs <- c('doParallel', 'foreach')
lapply(pkgs, require, character.only = T)
registerDoParallel(cores = 4)
getDoParWorkers()

##Single sim check
bp.inf(init_inf=1,n=10,contact_data)
bp.inf(init_inf = 5,n=5,contact_data,rep_HH=FALSE)

##Multiple sims
n_gen <- 5
n_sims <- 10
# init_inf_v <- c(1,3,5,10)
init_inf_v <- 5
Y <- NULL
for (init_inf in init_inf_v) {
  tic <- Sys.time()
  X <- foreach(i=1:n_sims,.combine = 'rbind') %dopar% {
    x <- bp.inf(init_inf = init_inf,n=n_gen ,contact_data,rep_HH=FALSE)
    x$sim <- i
    x
  }
  toc <- Sys.time()
  toc-tic
  X$init_inf <- init_inf
  Y <- rbind(Y,X)
}

Y %>% 
  group_by(gen,init_inf) %>% 
  summarise(faded=sum(is.na(infections))/max(sim)) %>% 
  ggplot(aes(x=gen,y=faded,fill=factor(init_inf)))+
  geom_bar(stat="identity",position = "dodge")+
  geom_text(aes(label=round(faded,2),vjust=-.2))+
  ylim(0,1)

Y %>% subset(!is.na(infections)) %>% 
  # subset(gen<4) %>% head
  mutate(prop.new.contacts=new.contacts/contacts) %>% 
  gather("type","value",c(contacts,prop.new.contacts,infections)) %>% 
  ggplot(aes(x=gen,y=value,color=factor(init_inf),group=init_inf))+
  geom_line(aes(group=sim),alpha=.2)+
  stat_summary(geom = "line",fun.y = "mean",color="black")+
  facet_grid(type~init_inf,scales = "free_y")

Y %>% 
  mutate(prop_inf=infections/contacts) %>% 
  gather("type","value",c(contacts,infections,prop_inf,r0_mean,r_eff)) %>% 
  ggplot(aes(x=gen,y=value))+
  geom_line(aes(group=sim),alpha=.2)+
  stat_summary(geom = "line",fun.y = "mean",color="red")+
  facet_grid(type~init_inf,scales = "free_y")

Y %>% 
  mutate(prop_inf=infections/contacts) %>% 
  subset(gen==max(gen)) %>% 
  gather("type","value",c(contacts,infections,prop_inf,r0_mean,r_eff)) %>% 
  group_by(type) %>% 
  summarise(mean(value,na.rm = TRUE))

Y %>% 
  group_by(init_inf) %>% 
  summarise(sum(infections,na.rm = TRUE))

Y %>% 
  subset(gen==max(gen)) %>% 
  # subset(infections!=0) %>%
  mutate(prop.new.contacts=new.contacts/contacts) %>% 
  gather("type","value",c(contacts,new.contacts,infections,prop.new.contacts)) %>% 
  ggplot(aes(x=factor(init_inf),y=value,color=factor(init_inf)))+
  geom_boxplot()+
  # stat_summary(aes(x=.1,y=value,xintercept=stat(y)),
  #              geom = "vline",fun.y = mean,color="red")+
  # stat_summary(aes(x=.1,y=value,xintercept=stat(y)),
  #              geom = "vline",fun.y = median,color="green")+
  facet_wrap(type~.,scales = "free",nrow = 3)+
  coord_flip()
# coord_cartesian(ylim=c(0,n_sims))

Y %>% 
  subset(gen==max(gen)) %>% 
  subset(infections!=0) %>%
  # gather("type","value",c(contacts,infections)) %>% 
  ggplot(aes(x=infections,fill=factor(init_inf),group=factor(init_inf)))+
  geom_histogram(position = "dodge")+facet_grid(init_inf~.)

Y %>% group_by(gen,init_inf) %>% 
  summarise(value=mean(infections,na.rm = TRUE)) %>% 
  ggplot(aes(x=gen,y=value,color=factor(init_inf)))+geom_line()

Y %>% 
  subset(!is.na(infections)) %>% 
  gather("R","value",c('r_eff','r0_mean')) %>% 
  group_by(gen,init_inf,R) %>%
  summarise(value=mean(value,na.rm = TRUE)) %>%
  ggplot(aes(x=gen,y=value,color=R))+
  geom_line()+facet_wrap(init_inf~.)