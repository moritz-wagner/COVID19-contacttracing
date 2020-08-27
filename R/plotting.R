require(magrittr)
require(patchwork)
require(tidyverse)
require(data.table)
require(ggsci)

## Setup - which runs to plot
cap_weeks <- 8
sims <- 100
r0.lev <- "medium"
disp <- "medium"
iso_delay <- "short"
pre_symp <- "medium"
rel.inf.asym <- 1
pop <- 230*1000

# Setup - Which interventions to plot
interventions <- c("No Interventions",
                   "Isolation (all symptomatic)",
                   # "Isolation (all symptomatic) + no SSEs" ,
                   # "Shielding elderly (0 non-HH contacts for over 50 year olds)",
                   "School closures (0 non-HH contacts for 6-19 year olds)",
                   "Social Distancing (Max 5 non-HH contacts)",
                   "Tracing HH + 25% non-HH",
                   "Quarantine HH + 25% non-HH")

# setup R0
r0community <- ifelse(r0.lev=="high",3,ifelse(r0.lev=="medium",2.5,2))

# define the summary function
f <- function(x) {
  r <- quantile(x, probs = c(0.05, 0.25, 0.5, 0.75, 0.95))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

# Load pre-run data
load(paste0('R/runs/agesymp_r0',r0.lev,'_disp',disp,'_relasym',rel.inf.asym*100,'_delay_',iso_delay,'_pre_',pre_symp,'_init5_',sims,'sims_',pop/1000,'pop.RData'))

#Create a folder to save figures for the run
dir.create(path = paste0('figures/r0',r0.lev,'_disp',disp,'_relasym',rel.inf.asym*100,'_delay_',iso_delay,'_pre_',pre_symp,'_init5_',sims,'sims_',pop/1000,'pop/'), showWarnings = FALSE)

# Some data wrangling
res %<>% mutate_if(is.numeric, list(~replace_na(., 0))) %>% subset(intervention%in%interventions)
res.gen %<>% subset(intervention%in%interventions)

levels(res$Intervention) <- c("No Interventions",
                              "Isolation (all symptomatic)",
                              "Isolation (all symptomatic)\n+ no SSEs" ,
                              "Shielding elderly\n(0 non-HH contacts for over 50 year olds)",
                              "School closures\n(0 non-HH contacts for 6-19 year olds)",
                              "Physical Distancing\n(Max 5 non-HH contacts)",
                              "Tracing HH + 25% non-HH",
                              "Quarantine HH + 25% non-HH")

levels(res.gen$Intervention) <- c("No Interventions",
                                  "Isolation (all symptomatic)",
                                  "Isolation (all symptomatic)\n+ no SSEs" ,
                                  "Shielding elderly\n(0 non-HH contacts for over 50 year olds)",
                                  "School closures\n(0 non-HH contacts for 6-19 year olds)",
                                  "Physical Distancing\n(Max 5 non-HH contacts)",
                                  "Tracing HH + 25% non-HH",
                                  "Quarantine HH + 25% non-HH")

# Create csv file of outputs---------
res %>% mutate(prop.hh=ifelse(hh_cases>0,hh_cases/weekly_cases,0)) %>% 
  gather("measure","value",-c(week,sim,age,intervention,Intervention,asym)) %>% 
  group_by(week,measure,sim,intervention,Intervention) %>% 
  summarise(value=sum(value)) -> tab
as.data.table(tab)[,as.list(round(f(value))),by=.(week,Intervention,measure)] %>% 
  write.csv(file = paste0('figures/r0',r0.lev,'_disp',disp,'_relasym',rel.inf.asym*100,'_delay_',iso_delay,'_pre_',pre_symp,'_init5_',sims,'sims_',pop/1000,
                          'pop/r0',r0.lev,'_disp',disp,'_relasym',rel.inf.asym*100,'_delay_',iso_delay,'_pre_',pre_symp,'_init5_',sims,'sims_',pop/1000,'pop.csv'))

# Plotting----------
means <- res.gen[,.(mean=mean(effective_r0)),by=.(Intervention)]

### Main plot (equivalent to figure 2 in the MS)-------
# Cumulative cases
res %>% 
  group_by(week,sim,Intervention) %>% 
  summarise_if(is.numeric,sum) %>% 
  subset(week<=cap_weeks) %>% 
  ggplot(aes(x=factor(week),y=cumulative,fill=Intervention))+
  stat_summary(fun.data = f,geom="boxplot",alpha=.5,position = "dodge")+
  theme_minimal()+
  ggtitle("Cumulative Cases")+xlab("week")+ylab("Weekly Cases") -> p.cases.cum

# HH contacts isolated
res %>% 
  group_by(week,sim,Intervention) %>% 
  summarise_if(is.numeric,sum) %>% 
  subset(week<=cap_weeks) %>%
  mutate(hh.con=ifelse(Intervention%in%c('Tracing HH','Tracing HH + 25% non-HH',
                                         'Quarantine HH','Quarantine HH + 25% non-HH',
                                         'Tracing HH + Social Distancing (50% non-HH contacts)',
                                         'Tracing HH + Shielding elderly (25% non-HH contacts for over 50 year olds)',
                                         'Quarantine HH + Social Distancing (50% non-HH contacts)',
                                         'Quarantine HH + Shielding elderly (25% non-HH contacts for over 50 year olds)'),
                       hh.con,0)) %>%  
  ggplot(aes(x=factor(week),y=hh.con,fill=Intervention))+
  stat_summary(fun.data = f,geom="boxplot",alpha=.5,position = "dodge")+
  theme_minimal()+
  ggtitle("HH contacts isolated")+xlab("week")+ylab("Weekly Contacts") -> p.hh.contacts

# non-HH contacts isolated
res %>% 
  group_by(week,sim,Intervention) %>% 
  summarise_if(is.numeric,sum) %>% 
  subset(week<=cap_weeks) %>%
  mutate(nhh.con=ifelse(Intervention%in%c('Tracing HH + 25% non-HH','Quarantine HH + 25% non-HH'),nhh.con*.25,0)) %>%
  ggplot(aes(x=factor(week),y=nhh.con,fill=Intervention))+
  stat_summary(fun.data = f,geom="boxplot",alpha=.5,position = "dodge")+
  theme_minimal()+
  ggtitle("Non-HH contacts isolated")+xlab("week")+ylab("Weekly Contacts")-> p.nhh.contacts

# Effective reproduction number
means <- res.gen[,.(effective_r0=mean(effective_r0)),by=.(Intervention,sim)] 
means %>% 
  ggplot(aes(x=1,y=effective_r0,fill=Intervention))+
  stat_summary(fun.data = f,geom="boxplot",alpha=.5,position="dodge")+
  geom_hline(yintercept = 1,color="black",size=1,alpha=.5,linetype=2)+
  geom_hline(yintercept = r0community,color="black",size=1,alpha=.5,linetype=2)+
  theme_minimal()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ggtitle("Effective Reproduction number")+
  ylab("R effective") -> R_eff

# Proportion of cases that die out
res %>% 
  group_by(week,sim,Intervention) %>% 
  summarise_if(is.numeric,sum) %>%   subset(week==cap_weeks) %>% 
  group_by(Intervention) %>% 
  summarise(extinct=sum(cumulative<control_max_cases & weekly_cases==0)/max(sim)) %>% 
  gather("measure","prop",-1) %>% 
  ggplot(aes(x=Intervention,y=prop,fill=Intervention))+
  geom_bar(stat="identity",position = "dodge",alpha=.5,color="black")+
  geom_text(aes(label=round(prop,2)),vjust=-0.2,size=2)+
  ggtitle(paste0("Outbreaks that go\nextinct in 8\nweeks"))+
  theme_minimal() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none")+
  xlab("")+
  ylab("Proportion")+ylim(0,1) -> p1

# Proportion of cases below 1000
control_max_cases <- 1000
res %>% 
  group_by(week,sim,Intervention) %>% 
  summarise_if(is.numeric,sum) %>% 
  subset(week==cap_weeks) %>% 
  group_by(Intervention) %>% 
  summarise(less=sum(cumulative<control_max_cases)/max(sim)) %>% 
  gather("measure","prop",-1) %>% 
  ggplot(aes(x=Intervention,y=prop,fill=Intervention))+
  geom_bar(stat="identity",position = "dodge",alpha=.5,color="black")+
  geom_text(aes(label=round(prop,2)),vjust=-0.1,size=2)+
  ggtitle(paste0("Outbreaks with\n<",control_max_cases," cases in 8\nweeks"))+
  theme_minimal() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none")+
  xlab("")+
  ylab("Proportion")+ylim(0,1) -> p2

png(paste0('figures/r0',r0.lev,'_disp',disp,'_relasym',rel.inf.asym*100,'_delay_',iso_delay,'_pre_',pre_symp,'_init5_',sims,'sims_',pop/1000,
           'pop/main_r0',r0.lev,'_disp',disp,'_relasym',rel.inf.asym*100,'_delay_',iso_delay,'_pre_',pre_symp,'_init5_',sims,'sims_',pop/1000,'pop.png'),width = 12,height = 7,units = "in",res=300)

((p.cases.cum/(p.hh.contacts/p.nhh.contacts))|(R_eff/(p1|p2))) + plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = c('A'),tag_suffix = ')') & theme(legend.key.size = unit(1.0, 'cm')) &
  scale_fill_nejm()

dev.off()

###Cases by age-------
png(paste0('figures/r0',r0.lev,'_disp',disp,'_relasym',rel.inf.asym*100,'_delay_',iso_delay,'_pre_',pre_symp,'_init5_',sims,'sims_',pop/1000,
           'pop/ageprop_r0',r0.lev,'_disp',disp,'_relasym',rel.inf.asym*100,'_delay_',iso_delay,'_pre_',pre_symp,'_init5_',sims,'sims_',pop/1000,'pop.png'),
    width = 10,height = 7,units = "in",res=300)

as.data.table(res)[!is.na(age) & week==8,] %>% 
  group_by(sim,Intervention,age) %>% 
  summarise_if(is.numeric,sum) %>% 
  group_by(sim,Intervention) %>% 
  mutate(prop=cumulative/sum(cumulative)) %>% 
  ggplot(aes(x=factor(age),y=prop,fill=Intervention))+
  # geom_boxplot(alpha=.5)+
  # geom_boxplot(alpha=.5,outlier.shape = NA)+
  stat_summary(fun.data = f,geom="boxplot",alpha=.5,position = "dodge")+
  theme_minimal()+
  # facet_wrap(.~Intervention,scales="free")
  # coord_cartesian(ylim = c(0,10000))+
  # stat_summary(aes(group=intervention),fun.y = mean,geom = "line")
  ggtitle("Proportion of cumulative cases by age in 8 weeks (all cases)")+xlab("Age")+
  ylab("Proportion")+
  theme(legend.key.size = unit(1.0, 'cm'))+
  scale_fill_nejm()

dev.off()

###Symp Cases by age-------
png(paste0('figures/r0',r0.lev,'_disp',disp,'_relasym',rel.inf.asym*100,'_delay_',iso_delay,'_pre_',pre_symp,'_init5_',sims,'sims_',pop/1000,
           'pop/symp_ageprop_r0',r0.lev,'_disp',disp,'_relasym',rel.inf.asym*100,'_delay_',iso_delay,'_pre_',pre_symp,'_init5_',sims,'sims_',pop/1000,'pop.png'),
    width = 10,height = 7,units = "in",res=300)

as.data.table(res)[!is.na(age) & week==8 & !asym,] %>% 
  group_by(sim,Intervention,age) %>% 
  summarise_if(is.numeric,sum) %>% 
  group_by(sim,Intervention) %>% 
  mutate(prop=cumulative/sum(cumulative)) %>% 
  ggplot(aes(x=factor(age),y=prop,fill=Intervention))+
  # geom_boxplot(alpha=.5)+
  # geom_boxplot(alpha=.5,outlier.shape = NA)+
  stat_summary(fun.data = f,geom="boxplot",alpha=.5,position = "dodge")+
  theme_minimal()+
  # facet_wrap(.~Intervention,scales="free")
  # coord_cartesian(ylim = c(0,10000))+
  # stat_summary(aes(group=intervention),fun.y = mean,geom = "line")
  ggtitle("Proportion of cumulative cases by age in 8 weeks (symptomatic)")+xlab("Age")+
  ylab("Proportion")+
  theme(legend.key.size = unit(1.0, 'cm'))+
  scale_fill_nejm()
dev.off()

### Quarantine rate vs physical distancing (equivalent to for figure 3)------------
load(paste0('R/runs/agesymp_quarantineSD2_r0',r0.lev,'_disp',disp,'_relasym',rel.inf.asym*100,'_delay_',iso_delay,'_pre_',pre_symp,'_init5_',sims,'sims_',pop/1000,'.RData'))

#Data wrangling
res.qu <-  bind_rows(lapply(res.ext,function(x) x$res.week %>% 
                              group_by(week,sim,prop.ascertain,lim.non.HH) %>% 
                              summarise_if(is.numeric,sum) %>% 
                              group_by(sim,prop.ascertain,lim.non.HH) %>%
                              mutate(cumulative=cumsum(weekly_cases))))
res.qu$intervention <- "Quarantine + SD"
res.gen.qu <-  bind_rows(lapply(res.ext,function(x) x$res.gen))
res.gen.qu$intervention <- "Quarantine + SD"

res <- res.qu
res.gen <- res.gen.qu

res %<>% ungroup() %>% mutate(week=week+1) 
res %<>% mutate_if(is.numeric, list(~replace_na(., 0)))
res$lim.non.HH <- factor(res$lim.non.HH,levels=rev(unique(res$lim.non.HH)))
res.gen$lim.non.HH <- factor(res.gen$lim.non.HH,levels=rev(unique(res.gen$lim.non.HH)))

#csv file of outputs
res %>% 
  mutate(nhh.con.iso=nhh.con*prop.ascertain,
         con.iso=hh.con+nhh.con) %>% 
  gather("measure","value",-c(week,sim,prop.ascertain,lim.non.HH,intervention)) -> tab
as.data.table(tab)[,as.list(round(f(value))),by=.(week,prop.ascertain,lim.non.HH,measure,intervention)] %>% 
  write.csv(file = paste0('figures/r0',r0.lev,'_disp',disp,'_relasym',rel.inf.asym*100,'_delay_',iso_delay,'_pre_',pre_symp,'_init5_',sims,'sims_',pop/1000,
                          'pop/fullqu_table_r0',r0.lev,'_disp',disp,'_relasym',rel.inf.asym*100,'_delay_',iso_delay,'_pre_',pre_symp,'_init5_',sims,'sims_',pop/1000,'pop.csv'))
#Main plot (equivalent to figure 3)-------
# Cumulative cases after 8 weeks
res %>% 
  subset(week==max(week)) %>%
  ggplot(aes(x=factor(prop.ascertain),y=cumulative,fill=factor(lim.non.HH),
             group=interaction(lim.non.HH,prop.ascertain,week)))+
  stat_summary(fun.data = f,geom="boxplot",alpha=.5,position="dodge")+
  ylab("")+
  xlab("Proportion of non-HH contacts traced/quarantined")+
  ggtitle("Cumulative number of cases after 8 weeks")+
  labs(fill="Physical Distancing:\nMaximum number of\nnon-HH contacts")+
  theme_minimal() -> p_cases

# HH contacts isolated in week 8
res %>% 
  subset(week==max(week)) %>% 
  ggplot(aes(x=factor(prop.ascertain),y=hh.con ,fill=factor(lim.non.HH),
             group=interaction(lim.non.HH,prop.ascertain,week)))+
  stat_summary(fun.data = f,geom="boxplot",alpha=.5,position="dodge")+
  ylab("")+
  xlab("Proportion of non-HH\ncontacts traced/quarantined")+
  ggtitle("HH contacts\nisolated in week 8")+
  labs(fill="Physical Distancing:\nMaximum number of\nnon-HH contacts")+
  theme_minimal() -> p_hhcon

# non-HH contacts isolated in week 8
res %>% 
  mutate(nhh.con=nhh.con*prop.ascertain,
         con=hh.con+nhh.con) %>% 
  subset(week==max(week)) %>% 
  ggplot(aes(x=factor(prop.ascertain),y=nhh.con ,fill=factor(lim.non.HH),
             group=interaction(lim.non.HH,prop.ascertain,week)))+
  stat_summary(fun.data = f,geom="boxplot",alpha=.5,position="dodge")+
  ylab("")+
  xlab("Proportion of non-HH\ncontacts traced/quarantined")+
  ggtitle("Non-HH contacts\nisolated in week 8")+
  labs(fill="Physical Distancing:\nMaximum number of\nnon-HH contacts")+
  theme_minimal() -> p_nhhcon

# Effective reproduction number
means <- res.gen[,.(effective_r0=mean(effective_r0)),by=.(prop.ascertain,intervention,lim.non.HH,sim)] 
means %>% 
  ggplot(aes(x=factor(prop.ascertain),y=effective_r0,fill=factor(lim.non.HH),
             group=interaction(intervention,lim.non.HH,prop.ascertain)))+
  stat_summary(fun.data = f,geom="boxplot",alpha=.5,position="dodge")+
  geom_hline(yintercept = 1,color="black",size=1,alpha=.5,linetype=2)+
  geom_hline(yintercept = r0community,color="black",size=1,alpha=.5,linetype=2)+
  theme_minimal()+
  ylab("R effective")+
  xlab("Proportion of non-HH contacts traced/quarantined")+
  ggtitle("Effective Reproduction number")+
  labs(fill="Physical Distancing:\nMaximum number of\nnon-HH contacts") -> R_eff

# Proportion of outbreaks that go extinct
res %>% 
  group_by(intervention,prop.ascertain,lim.non.HH) %>% 
  subset(week==cap_weeks) %>% 
  summarise(extinct=sum(cumulative<control_max_cases & weekly_cases==0)/max(sim)) %>% 
  gather("measure","prop",-c(1,2,3)) %>% 
  ggplot(aes(x=factor(prop.ascertain),y=prop,fill=lim.non.HH))+
  geom_bar(stat="identity",position = "dodge",alpha=.5,color="black")+
  # geom_text(aes(label=round(prop,2)),vjust=-0.1,size=2,position=position_dodge(width=1))+
  theme_minimal() + 
  theme(legend.position = "none")+  
  ggtitle(paste0("Outbreaks that go\nextinct in 8 weeks"))+
  labs(fill="Reduction in\nnon-HH contacts")+
  xlab("Proportion of non-HH\ncontacts traced/quarantined")+
  ylab("Proportion")+ylim(0,1) -> p1

# Proportion of outbreaks that remain below 1000
control_max_cases <- 1000
res %>% 
  group_by(intervention,prop.ascertain,lim.non.HH) %>% 
  subset(week==cap_weeks) %>% 
  summarise(less=sum(cumulative<control_max_cases)/max(sim)) %>% 
  gather("measure","prop",-c(1,2,3)) %>% 
  ggplot(aes(x=factor(prop.ascertain),y=prop,fill=lim.non.HH))+
  geom_bar(stat="identity",position = "dodge",alpha=.5,color="black")+
  # geom_text(aes(label=round(prop,2)),vjust=-0.1,size=2,position=position_dodge(width=1))+
  ggtitle(paste0("Outbreaks with <",control_max_cases,"\ncases in 8 weeks"))+
  theme_minimal() + 
  theme(legend.position = "none")+
  labs(fill="Reduction in\nnon-HH contacts")+
  xlab("Proportion of non-HH\ncontacts traced/quarantined")+
  ylab("Proportion")+ylim(0,1) -> p2

png(paste0('figures/r0',r0.lev,'_disp',disp,'_relasym',rel.inf.asym*100,'_delay_',iso_delay,'_pre_',pre_symp,'_init5_',sims,'sims_',pop/1000,
           'pop/fullqu_r0',r0.lev,'_disp',disp,'_relasym',rel.inf.asym*100,'_delay_',iso_delay,'_pre_',pre_symp,'_init5_',sims,'sims_',pop/1000,'pop.png'),width = 12,height = 7,units = "in",res=300)
((p_cases/(p_hhcon|p_nhhcon))|(R_eff/(p1|p2)))+plot_layout(guides = 'collect')+
  plot_annotation(tag_levels = c('A'),tag_suffix = ')') & theme(legend.key.size = unit(1.0, 'cm')) &
  scale_fill_brewer(palette="Blues")
dev.off()
