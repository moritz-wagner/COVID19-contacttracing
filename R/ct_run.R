library(readr)
library(tidyverse)
library(ggplot2)
library(data.table)

source('R/ct_functions.R')

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
new.inf(sample(contact_data$csid,1),cd,r0 = 10,HH = TRUE,n.inf.HH=3,HH.herd.imm = FALSE)

library(doParallel)
pkgs <- c('doParallel', 'foreach')
lapply(pkgs, require, character.only = T)
registerDoParallel(cores = 4)
getDoParWorkers()

##Single sim check
bp.inf(init_inf=1,n=10,contact_data)
bp.inf(init_inf = 1,n=5,contact_data = contact_data,HH.herd.imm = FALSE)

##Multiple sims
n_gen <- 5
n_sims <- 100
# init_inf_v <- c(1,3,5,10)
init_inf_v <- 5
Y <- NULL
for (init_inf in init_inf_v) {
  tic <- Sys.time()
  X <- foreach(i=1:n_sims,.combine = 'rbind') %dopar% {
    x <- bp.inf(init_inf = init_inf,n=n_gen ,contact_data,HH.herd.imm=FALSE)
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
  gather("type","value",c(contacts,infections,prop_inf,r0_mean,r_eff,max.HH.inf,mean.HH.inf)) %>% 
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


Y %>% 
  select(max.HH.inf,mean.HH.inf) %>% 
  ggplot(aes(alpha=.2))+
  geom_histogram(aes(max.HH.inf),fill="red")+
  geom_histogram(aes(mean.HH.inf),fill="green")+
  theme(legend.position = "none")+ggtitle("Max and mean number of infecteds in a HH\nafter 5 generations")



###counting contacts
contact_data$ever_met %>% summary() %>% barplot()
contact_data$how_often %>% summary() %>% barplot()
contact_data$how_often <- factor(contact_data$how_often,levels = unique(contact_data$how_often)[c(3,2,4,5,1)])
as.data.table(contact_data)[,.N,by=c('ever_met','how_often')] %>% 
  ggplot(aes(how_often,N,fill=ever_met))+geom_bar(stat="identity")
as.data.table(contact_data)[,.N,by=c('ever_met','how_often')][,prop:=N/sum(N),by=how_often] %>% 
  ggplot(aes(how_often,prop,fill=ever_met))+geom_bar(stat="identity")

C <- NULL
for (i in 1:50) {
  cd <- contact_boot(contact_data)
  id <- sample(cd$csid,1)
  t <- 10
  Ci <- count.contacts(id,cd,t)
  Ci$run <- i
  C <- rbind(C,Ci)
}

C %>% 
  group_by(how_often,ever_met,run) %>% 
  summarise(weight=sum(weight,na.rm = TRUE)) %>% 
  ggplot(aes(x=how_often,y=weight,fill=ever_met,group=factor(run)))+
  geom_bar(stat="identity",position = "dodge")

n <- 500

X <- foreach(t = 1:14,.combine = 'rbind') %dopar% {
  X <- times(n) %dopar% {
    cd <- contact_boot(contact_data)
    id <- sample(cd$csid,1)
    Ci <- count.contacts(id,cd,t)
    nrow(Ci)
  }
  X <- data.frame(contacts=X,days=t)
  X
}

# n <- 100
# Y <- NULL
# for (t in 1:14) {
#   for (i in 1:n) {
#     cd <- contact_boot(contact_data)
#     id <- sample(cd$csid,1)
#     Ci <- count.contacts(id,cd,t) 
#     Y <- rbind(Y,data.frame(contacts=nrow(Ci),days=t))
#   }
# }

X %>% 
  ggplot(aes(x=factor(days),y=contacts))+geom_boxplot()+
  ggtitle("Number of new contacts in 14 days (500 sims)")+
  xlab('days')+ylab("Unique contacts")


