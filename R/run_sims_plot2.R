require(magrittr)
require(patchwork)

## Run outbreak sims and do some plotting
source('R/contact_functions.R')
source('R/outbreak_functions2.R')

##Plot of parameters used-------------
n <- 10000
exposure <- 5
#Incubation periods
incfn <- dist_setup(dist_shape = 2.322737, dist_scale = 6.492272)
inc.period <- incfn(n)+exposure
mean.inc <- mean(inc.period)


#Serial intervals based on incubation samples
k = 0.7
inc.period.sample <- sample(inc.period,1)
# inc.period.sample <- max(inc.period)
serial.interval <- inf_fn(rep(inc.period.sample-exposure,n),k)+exposure
# serial.interval <- inf_fn(rep(inc.period.sample,n),k)
sum(serial.interval<=exposure)/n
mean.serial.interval <- mean(serial.interval)
median.serial.interval <- median(serial.interval)
perc.pre.symp <- sum(serial.interval<inc.period.sample)/n


#Delay to isolation long and short
delay_shape = c(1.651524, 2.305172)
delay_scale = c(4.287786, 9.483875)
delay_shape = c(2.5,1.651524, 2.305172)
delay_scale = c(5,4.287786, 9.483875)

iso.delay <- purrr::map2(delay_shape, 
                         delay_scale, function(x, y) {
                           delayfn <- dist_setup(x,y);
                           delayfn(n)
                         })
mean.iso.delay <- sapply(iso.delay,mean)
median.iso.delay <- sapply(iso.delay,median)
max.iso <- max(sapply(iso.delay,max))


##R0
r0community = 2.5
disp.com = 0.16
r0 <- rnbinom(n, size = disp.com, mu = r0community)
mean.r0 <- mean(r0)

params <- data.frame(n=1:n,
                     # r0,
                     iso.delay,inc.period,serial.interval)


#Plots
data.frame(n=1:n,days=inc.period) %>% 
  ggplot(aes(x=days))+
  geom_density(fill="red",alpha=.5)+
  # geom_vline(aes(xintercept=mean.inc))+
  geom_vline(aes(xintercept=inc.period.sample),color="red")+
  geom_vline(aes(xintercept=exposure),color="yellow")+
  xlim(0,max.iso)+
  ggtitle('Incubation Period of infector (time until onset of symptoms)') -> p1

data.frame(n=1:n,days=serial.interval) %>% 
  ggplot(aes(x=days))+
  geom_density(fill="blue",alpha=.5)+
  # geom_vline(aes(xintercept=mean.serial.interval))+
  # geom_vline(aes(xintercept=median.serial.interval),color="blue")+
  geom_vline(aes(xintercept=inc.period.sample),color="red")+
  geom_vline(aes(xintercept=exposure),color="yellow")+
  xlim(0,max.iso)+
  ggtitle('Exposure time of infected individual (Serial interval)')+
  annotate("text", 
           y= .2,label=paste0(round(perc.pre.symp*100),"% presymp"), x =max.iso,hjust=1)-> p2
data.frame(n=1:n,current=iso.delay[[1]],short=iso.delay[[2]],long=iso.delay[[3]]) %>% 
  gather('delay','days',-1) %>% 
  ggplot(aes(x=days,fill=delay,group=delay))+
  geom_density(alpha=.5)+
  # geom_vline(xintercept=mean.iso.delay)+
  # geom_vline(xintercept=median.iso.delay,color="blue")+
  geom_vline(xintercept=exposure,color="yellow")+
  geom_vline(xintercept=inc.period.sample,color="red")+
  xlim(0,max.iso)+theme(legend.position = c(.8,.8))+
  ggtitle('Delay to isolation of infector') -> p3

data.frame(n=1:n,sec.infs=r0) %>% 
  count(sec.infs) %>% mutate(density=n/sum(n)) %>% 
  ggplot(aes(x=(sec.infs),y=density))+
  geom_bar(stat="identity")+
  coord_cartesian(xlim=c(0,20))+
  # geom_histogram(fill="black",alpha=.5,aes(y = ..density..))+
  geom_vline(aes(xintercept=mean.r0))+ggtitle('R0')+xlab("") -> p4

png('plots/params2.png',width = 10,height = 5,units = "in",res=300)
# ((p1/p2/p3)) & theme_minimal()
((p1/p2/p3)|p4) & theme_minimal()
dev.off()

# 
# params %>% 
#   gather('param','value',-1) %>% 
#   group_by(param) %>% 
#   mutate(mean=mean(value)) %>% 
#   ggplot(aes(x=value,fill=param))+
#   geom_density()+
#   geom_vline(aes(xintercept=mean),color="red")+
#   facet_grid(param~.)
# 
#   



## Run outbreak simulations for different scenarios----------
##Setup-----
## Create synthetic population based on contact data
freq = c("Daily","Often","Regularly","Rarely","Never") #Can adjust this depending on which types of contacts to include
p.urban = .3
t <- 7
pop <- 300000
time.cd <- system.time(CD <- count.contacts3(contact_data,t=t,p.urban = p.urban,freq=freq,pop=pop))
setnames(CD[[1]],old = "caseid",new = "csid")
sims <- 50
cap_cases <- 2000
cap_weeks <- 8
cap_max_days <- (cap_weeks-1)*7
cap_contacts_per_day <- Inf
cap_max_hhcon <- cap_contacts_per_day*cap_max_days
cap_max_nhhcon <- cap_contacts_per_day*cap_max_days

num.initial.cases <- 5
r0community = 2.5;disp.com = 0.16
r0isolated = 0;disp.iso = 1
k = 0.7 ##DETERMINES PROPORTION OF PRE-SYMPTOMIC TRANSMISION, need to work out how exactly
delay_shape = c(2.5,1.651524, 2.305172)
delay_scale = c(5,4.287786, 9.483875)
delay_shape = 2.5;delay_scale = 5
delay_shape = 1.651524;delay_scale = 4.287786
prop.asym = 0.2

# No Isolation (assumes prop asymtomatic = 1)-----
quarantine <- FALSE
HH_trace <- FALSE
prop.ascertain <- 0

time.ext.par<-system.time(res.ext <- scenario_sim_ext_parallel(n.sim = sims,num.initial.cases = num.initial.cases, 
                                                               prop.ascertain = prop.ascertain, cap_max_days = cap_max_days, 
                                                               cap_cases = cap_cases, cap_max_hhcon=cap_max_hhcon, cap_max_nhhcon=cap_max_nhhcon, 
                                                               r0isolated = r0isolated, r0community = r0community, 
                                                               disp.iso = disp.iso, disp.com = disp.com, delay_shape = delay_shape, 
                                                               delay_scale = delay_scale, k = k, prop.asym = 1, 
                                                               quarantine = quarantine,HH_trace = HH_trace,
                                                               contact_data=CD,p.urban = p.urban,
                                                               t=t,freq = freq,
                                                               pop=pop))


time.ext.par

res.ext$res.week %>%
  # select(-intervention) %>%
  # mutate(hh.con=hh.con*HH_trace,
  #        nhh.con=nhh.con*prop.ascertain) %>%
  gather("measure","value",-c("week","sim")) %>%
  ggplot(aes(x=week,y=value,color=measure))+
  geom_boxplot(aes(group=week))+
  geom_line(aes(group=sim),alpha=.1,color="red")+
  stat_summary(fun.y = mean,geom="line")+
  facet_wrap(~measure,scales="free")

res.ext$res.gen %>%
  group_by(sim) %>% mutate(cum_cases_in_gen=cumsum(cases_in_gen)) %>%
  gather("measure","value",-c("gen","sim")) %>%
  ggplot(aes(x=gen,y=value,color=measure))+
  geom_boxplot(aes(group=gen))+
  geom_line(aes(group=sim),alpha=.1,color="red")+
  stat_summary(fun.y = mean,geom="line")+
  facet_wrap(~measure,scales="free")

res.no.int <- res.ext$res.week
res.no.int$intervention <- "No Interventions"

res.no.int.gen <- res.ext$res.gen
res.no.int.gen$intervention <- "No Interventions"

# Isolation only severe cases (i.e. 80% asymp)-----
quarantine <- FALSE
HH_trace <- FALSE
prop.ascertain <- 0

time.ext.par<-system.time(res.ext <- scenario_sim_ext_parallel(n.sim = sims,num.initial.cases = num.initial.cases, 
                                                               prop.ascertain = prop.ascertain, cap_max_days = cap_max_days, 
                                                               cap_cases = cap_cases, cap_max_hhcon=cap_max_hhcon, cap_max_nhhcon=cap_max_nhhcon, 
                                                               r0isolated = r0isolated, r0community = r0community, 
                                                               disp.iso = disp.iso, disp.com = disp.com, delay_shape = delay_shape, 
                                                               delay_scale = delay_scale, k = k, prop.asym = .8, 
                                                               quarantine = quarantine,HH_trace = HH_trace,
                                                               contact_data=CD,p.urban = p.urban,
                                                               t=t,freq = freq,
                                                               pop=pop))

time.ext.par

# res.ext$res.week %>% 
#   # select(-intervention) %>% 
#   # mutate(hh.con=hh.con*HH_trace,
#   #        nhh.con=nhh.con*prop.ascertain) %>% 
#   gather("measure","value",-c("week","sim")) %>% 
#   ggplot(aes(x=week,y=value,color=measure))+
#   geom_boxplot(aes(group=week))+
#   geom_line(aes(group=sim),alpha=.1,color="red")+
#   stat_summary(fun.y = mean,geom="line")+
#   facet_wrap(~measure,scales="free")
# 
# res.ext$res.gen %>% 
#   group_by(sim) %>% mutate(cum_cases_in_gen=cumsum(cases_in_gen)) %>% 
#   gather("measure","value",-c("gen","sim")) %>% 
#   ggplot(aes(x=gen,y=value,color=measure))+
#   geom_boxplot(aes(group=gen))+
#   geom_line(aes(group=sim),alpha=.1,color="red")+
#   stat_summary(fun.y = mean,geom="line")+
#   facet_wrap(~measure,scales="free")

res.no.tracing <- res.ext$res.week
res.no.tracing$intervention <- "Isolation only severe"

res.no.tracing.gen <- res.ext$res.gen
res.no.tracing.gen$intervention <- "Isolation only severe"

# Isolation only symptomatic (i.e. 20% asymp)-----
quarantine <- FALSE
HH_trace <- FALSE
prop.ascertain <- 0

time.ext.par<-system.time(res.ext <- scenario_sim_ext_parallel(n.sim = sims,num.initial.cases = num.initial.cases, 
                                                               prop.ascertain = prop.ascertain, cap_max_days = cap_max_days, 
                                                               cap_cases = cap_cases, cap_max_hhcon=cap_max_hhcon, cap_max_nhhcon=cap_max_nhhcon, 
                                                               r0isolated = r0isolated, r0community = r0community, 
                                                               disp.iso = disp.iso, disp.com = disp.com, delay_shape = delay_shape, 
                                                               delay_scale = delay_scale, k = k, prop.asym = prop.asym, 
                                                               quarantine = quarantine,HH_trace = HH_trace,
                                                               contact_data=CD,p.urban = p.urban,
                                                               t=t,freq = freq,
                                                               pop=pop))

time.ext.par

# res.ext$res.week %>% 
#   select(-intervention) %>% 
#   mutate(hh.con=hh.con*HH_trace,
#          nhh.con=nhh.con*prop.ascertain) %>% 
#   gather("measure","value",-c("week","sim")) %>% 
#   ggplot(aes(x=week,y=value,color=measure))+
#   geom_boxplot(aes(group=week))+
#   geom_line(aes(group=sim),alpha=.1,color="red")+
#   stat_summary(fun.y = mean,geom="line")+
#   facet_wrap(~measure,scales="free")
# 
# res.ext$res.gen %>% 
#   gather("measure","value",-c("gen","sim")) %>% 
#   ggplot(aes(x=gen,y=value,color=measure))+
#   geom_boxplot(aes(group=gen))+
#   geom_line(aes(group=sim),alpha=.1,color="red")+
#   stat_summary(fun.y = mean,geom="line")+
#   facet_wrap(~measure,scales="free")

res.no.tracing2 <- res.ext$res.week
res.no.tracing2$intervention <- "Isolation (all symptomatic)"

res.no.tracing.gen2 <- res.ext$res.gen
res.no.tracing.gen2$intervention <- "Isolation (all symptomatic)"

# Tracing HH-----
quarantine <- FALSE
HH_trace <- TRUE
prop.ascertain <- 0

time.ext.par<-system.time(res.ext <- scenario_sim_ext_parallel(n.sim = sims,num.initial.cases = num.initial.cases, 
                                                               prop.ascertain = prop.ascertain, cap_max_days = cap_max_days, 
                                                               cap_cases = cap_cases, cap_max_hhcon=cap_max_hhcon, cap_max_nhhcon=cap_max_nhhcon, 
                                                               r0isolated = r0isolated, r0community = r0community, 
                                                               disp.iso = disp.iso, disp.com = disp.com, delay_shape = delay_shape, 
                                                               delay_scale = delay_scale, k = k, prop.asym = prop.asym, 
                                                               quarantine = quarantine,HH_trace = HH_trace,
                                                               contact_data=CD,p.urban = p.urban,
                                                               t=t,freq = freq,
                                                               pop=pop))

time.ext.par
# res.ext$res.week %>%
#   mutate(hh.con=hh.con*HH_trace,
#          nhh.con=nhh.con*prop.ascertain) %>%
#   gather("measure","value",-c("week","sim")) %>%
#   ggplot(aes(x=week,y=value,color=measure))+
#   geom_boxplot(aes(group=week))+
#   geom_line(aes(group=sim),alpha=.1,color="red")+
#   stat_summary(fun.y = mean,geom="line")+
#   facet_wrap(~measure,scales="free")
# 
# res.ext$res.gen %>%
#   gather("measure","value",-c("gen","sim")) %>%
#   ggplot(aes(x=gen,y=value,color=measure))+
#   geom_boxplot(aes(group=gen))+
#   geom_line(aes(group=sim),alpha=.1,color="red")+
#   stat_summary(fun.y = mean,geom="line")+
#   facet_wrap(~measure,scales="free")

res.tracing.HH <- res.ext$res.week
res.tracing.HH$intervention <- "Tracing HH"

res.tracing.HH.gen <- res.ext$res.gen
res.tracing.HH.gen$intervention <- "Tracing HH"


# Tracing HH + 25% non-HH-----
quarantine <- FALSE
HH_trace <- TRUE
prop.ascertain <- 0.25

time.ext.par<-system.time(res.ext <- scenario_sim_ext_parallel(n.sim = sims,num.initial.cases = num.initial.cases, 
                                                               prop.ascertain = prop.ascertain, cap_max_days = cap_max_days, 
                                                               cap_cases = cap_cases, cap_max_hhcon=cap_max_hhcon, cap_max_nhhcon=cap_max_nhhcon, 
                                                               r0isolated = r0isolated, r0community = r0community, 
                                                               disp.iso = disp.iso, disp.com = disp.com, delay_shape = delay_shape, 
                                                               delay_scale = delay_scale, k = k, prop.asym = prop.asym, 
                                                               quarantine = quarantine,HH_trace = HH_trace,
                                                               contact_data=CD,p.urban = p.urban,
                                                               t=t,freq = freq,
                                                               pop=pop))

time.ext.par
# res.ext$res.week %>% 
#   mutate(hh.con=hh.con*HH_trace,
#     nhh.con=nhh.con*prop.ascertain) %>% 
#   gather("measure","value",-c("week","sim")) %>% 
#   ggplot(aes(x=week,y=value,color=measure))+
#   geom_boxplot(aes(group=week))+
#   geom_line(aes(group=sim),alpha=.1,color="red")+
#   stat_summary(fun.y = mean,geom="line")+
#   facet_wrap(~measure,scales="free")
# 
# res.ext$res.gen %>% 
#   gather("measure","value",-c("gen","sim")) %>% 
#   ggplot(aes(x=gen,y=value,color=measure))+
#   geom_boxplot(aes(group=gen))+
#   geom_line(aes(group=sim),alpha=.1,color="red")+
#   stat_summary(fun.y = mean,geom="line")+
#   facet_wrap(~measure,scales="free")

res.tracing.HH2 <- res.ext$res.week
res.tracing.HH2$intervention <- "Tracing HH + 25% non-HH"

res.tracing.HH2.gen <- res.ext$res.gen
res.tracing.HH2.gen$intervention <- "Tracing HH + 25% non-HH"

# Quaranting HH-----
quarantine <- TRUE
HH_trace <- TRUE
prop.ascertain <- 0

time.ext.par<-system.time(res.ext <- scenario_sim_ext_parallel(n.sim = sims,num.initial.cases = num.initial.cases, 
                                                               prop.ascertain = prop.ascertain, cap_max_days = cap_max_days, 
                                                               cap_cases = cap_cases, cap_max_hhcon=cap_max_hhcon, cap_max_nhhcon=cap_max_nhhcon, 
                                                               r0isolated = r0isolated, r0community = r0community, 
                                                               disp.iso = disp.iso, disp.com = disp.com, delay_shape = delay_shape, 
                                                               delay_scale = delay_scale, k = k, prop.asym = prop.asym, 
                                                               quarantine = quarantine,HH_trace = HH_trace,
                                                               contact_data=CD,p.urban = p.urban,
                                                               t=t,freq = freq,
                                                               pop=pop))

time.ext.par
# res.ext$res.week %>% 
#   mutate(hh.con=hh.con*HH_trace,
#          nhh.con=nhh.con*prop.ascertain) %>% 
#   gather("measure","value",-c("week","sim")) %>% 
#   ggplot(aes(x=week,y=value,color=measure))+
#   geom_boxplot(aes(group=week))+
#   geom_line(aes(group=sim),alpha=.1,color="red")+
#   stat_summary(fun.y = mean,geom="line")+
#   facet_wrap(~measure,scales="free")
# 
# res.ext$res.gen %>% 
#   gather("measure","value",-c("gen","sim")) %>% 
#   ggplot(aes(x=gen,y=value,color=measure))+
#   geom_boxplot(aes(group=gen))+
#   geom_line(aes(group=sim),alpha=.1,color="red")+
#   stat_summary(fun.y = mean,geom="line")+
#   facet_wrap(~measure,scales="free")

res.quarantine <- res.ext$res.week
res.quarantine$intervention <- "Quarantine HH"

res.quarantine.gen <- res.ext$res.gen
res.quarantine.gen$intervention <-  "Quarantine HH"


# Quaranting HH + 25% non-HH-----
quarantine <- TRUE
HH_trace <- TRUE
prop.ascertain <- 0.25

time.ext.par<-system.time(res.ext <- scenario_sim_ext_parallel(n.sim = sims,num.initial.cases = num.initial.cases, 
                                                               prop.ascertain = prop.ascertain, cap_max_days = cap_max_days, 
                                                               cap_cases = cap_cases, cap_max_hhcon=cap_max_hhcon, cap_max_nhhcon=cap_max_nhhcon, 
                                                               r0isolated = r0isolated, r0community = r0community, 
                                                               disp.iso = disp.iso, disp.com = disp.com, delay_shape = delay_shape, 
                                                               delay_scale = delay_scale, k = k, prop.asym = prop.asym, 
                                                               quarantine = quarantine,HH_trace = HH_trace,
                                                               contact_data=CD,p.urban = p.urban,
                                                               t=t,freq = freq,
                                                               pop=pop))

time.ext.par
# res.ext$res.week %>% 
#   mutate(hh.con=hh.con*HH_trace,
#          nhh.con=nhh.con*prop.ascertain) %>% 
#   gather("measure","value",-c("week","sim")) %>% 
#   ggplot(aes(x=week,y=value,color=measure))+
#   geom_boxplot(aes(group=week))+
#   geom_line(aes(group=sim),alpha=.1,color="red")+
#   stat_summary(fun.y = mean,geom="line")+
#   facet_wrap(~measure,scales="free")
# 
# res.ext$res.gen %>% 
#   gather("measure","value",-c("gen","sim")) %>% 
#   ggplot(aes(x=gen,y=value,color=measure))+
#   geom_boxplot(aes(group=gen))+
#   geom_line(aes(group=sim),alpha=.1,color="red")+
#   stat_summary(fun.y = mean,geom="line")+
#   facet_wrap(~measure,scales="free")

res.quarantine2 <- res.ext$res.week
res.quarantine2$intervention <- "Quarantine HH + 25% non-HH"

res.quarantine2.gen <- res.ext$res.gen
res.quarantine2.gen$intervention <-  "Quarantine HH + 25% non-HH"



##Scenarios to compare effec
# Tracing HH + p% non-HH-----
quarantine <- FALSE
HH_trace <- TRUE
prop.ascertain <- seq(0,1,length.out = 5)
res.ext <- list()
for (i in seq_along(prop.ascertain)) {
  time.ext.par<-system.time(res.ext[[i]] <- scenario_sim_ext_parallel(n.sim = sims,num.initial.cases = num.initial.cases, 
                                                                      prop.ascertain = prop.ascertain[i], cap_max_days = cap_max_days, 
                                                                      cap_cases = cap_cases, cap_max_hhcon=cap_max_hhcon, cap_max_nhhcon=cap_max_nhhcon, 
                                                                      r0isolated = r0isolated, r0community = r0community, 
                                                                      disp.iso = disp.iso, disp.com = disp.com, delay_shape = delay_shape, 
                                                                      delay_scale = delay_scale, k = k, prop.asym = prop.asym, 
                                                                      quarantine = quarantine,HH_trace = HH_trace,
                                                                      contact_data=contact_data,p.urban = p.urban,
                                                                      t=t,freq = freq,
                                                                      pop=pop))
  
  print(time.ext.par)
}

save(res.ext,file = 'R/traceHH_20asymp_2000cap_100sims.RData')

# Quaranting HH + p% non-HH-----
quarantine <- TRUE
HH_trace <- TRUE
prop.ascertain <- seq(0,1,length.out = 5)
res.ext <- list()
for (i in seq_along(prop.ascertain)) {
  time.ext.par<-system.time(res.ext[[i]] <- scenario_sim_ext_parallel(n.sim = sims,num.initial.cases = num.initial.cases, 
                                                                      prop.ascertain = prop.ascertain[i], cap_max_days = cap_max_days, 
                                                                      cap_cases = cap_cases, cap_max_hhcon=cap_max_hhcon, cap_max_nhhcon=cap_max_nhhcon, 
                                                                      r0isolated = r0isolated, r0community = r0community, 
                                                                      disp.iso = disp.iso, disp.com = disp.com, delay_shape = delay_shape, 
                                                                      delay_scale = delay_scale, k = k, prop.asym = prop.asym, 
                                                                      quarantine = quarantine,HH_trace = HH_trace,
                                                                      contact_data=contact_data,p.urban = p.urban,
                                                                      t=t,freq = freq,
                                                                      pop=pop))
  
  time.ext.par
}

save(res.ext,file = 'R/quarantineHH_20asymp_2000cap_100sims.RData')

## put the together and plot-----

res <- bind_rows(list(res.no.int,res.no.tracing,res.no.tracing2,res.tracing.HH,res.tracing.HH2,res.quarantine,res.quarantine2))
res.gen <- bind_rows(list(res.no.int.gen,res.no.tracing.gen,res.no.tracing.gen2,res.tracing.HH.gen,res.tracing.HH2.gen,res.quarantine.gen,res.quarantine2.gen))

res %<>% mutate(week=week+1) 

res$Intervention <- factor(res$intervention,levels = unique(res$intervention))
res.gen$Intervention <- factor(res.gen$intervention,levels = unique(res.gen$intervention))

figname <- 'plots/newCT_all_50sims_20asymp'

control_max_cases <- 1000
res %>% 
  group_by(Intervention) %>% 
  subset(week==cap_weeks) %>% 
  summarise(less=sum(cumulative<control_max_cases)/max(sim),
            extinct=sum(cumulative<control_max_cases & weekly_cases==0)/max(sim)) %>% 
  gather("measure","prop",-1) %>% 
  ggplot(aes(x=measure,y=prop,fill=Intervention))+
  geom_bar(stat="identity",position = "dodge")+
  ggtitle(paste0("Proportion of outbreaks that go extinct or with <",control_max_cases," cases in 8 weeks"))+
  theme(legend.position = "bottom",legend.direction = "vertical")+
  theme_minimal()+
  ylim(0,1) -> p

png(paste0(figname,'_Figure2.png'),width = 7,height = 5,units = "in",res=300)
p
dev.off()


res %>% 
  subset(week<=cap_weeks) %>%
  ggplot(aes(x=factor(week),y=weekly_cases,fill=Intervention,color=Intervention))+
  # geom_boxplot(alpha=.5)+
  geom_boxplot(alpha=.5,outlier.shape = NA)+
  # coord_cartesian(ylim = c(0,2000))+
  # stat_summary(aes(group=intervention),fun.y = mean,geom = "line")
  ggtitle("Weekly Incidence of Cases")+xlab("week")+ylab("Weekly Cases") -> p.cases

res %>% mutate(cumulative=ifelse(cumulative>cap_cases,cap_cases,cumulative))  %>% 
  subset(week<=cap_weeks) %>%
  ggplot(aes(x=factor(week),y=cumulative,fill=Intervention,color=Intervention))+
  # geom_boxplot(alpha=.5)+
  geom_boxplot(alpha=.5,outlier.shape = NA)+
  coord_cartesian(ylim = c(0,2000))+
  # stat_summary(aes(group=intervention),fun.y = mean,geom = "line")
  ggtitle("Cumulative Cases")+xlab("week")+ylab("Weekly Cases") -> p.cases.cum

res %>% 
  subset(week<=cap_weeks) %>%
  # subset(week<cap_weeks) %>% 
  mutate(hh.con=ifelse(intervention%in%c('No Interventions','Isolation only','Isolation only severe','Isolation (all symptomatic)'),0,hh.con)) %>% 
  # subset(intervention!="Isolation only") %>% 
  ggplot(aes(x=factor(week),y=hh.con,fill=Intervention,color=Intervention))+
  # geom_boxplot(alpha=.5)+
  geom_boxplot(alpha=.5,outlier.shape = NA)+
  # coord_cartesian(ylim = c(0,750))+
  ggtitle("HH contacts to trace")+xlab("week")+ylab("Weekly Contacts") -> p.hh.contacts

res %>% 
  subset(week<=cap_weeks) %>%
  # subset(week<cap_weeks) %>% 
  mutate(nhh.con=ifelse(intervention%in%c('Tracing HH + 25% non-HH','Quarantine HH + 25% non-HH'),nhh.con*.25,0)) %>% 
  # subset(intervention!="Isolation only") %>% 
  ggplot(aes(x=factor(week),y=nhh.con,fill=Intervention,color=Intervention))+
  # geom_boxplot(alpha=.5)+
  geom_boxplot(alpha=.5,outlier.shape = NA)+
  # coord_cartesian(ylim = c(0,4000))+
  ggtitle("Non-HH contacts to trace")+xlab("week")+ylab("Weekly Contacts")-> p.nhh.contacts

res %>% 
  subset(week<=cap_weeks) %>%
  # subset(week<cap_weeks) %>% 
  mutate(prop.hh=hh_cases/weekly_cases) %>% 
  # subset(intervention!="Isolation only") %>% 
  ggplot(aes(x=factor(week),y=prop.hh,fill=Intervention,color=Intervention))+
  # geom_boxplot(alpha=.5)+
  geom_boxplot(alpha=.5,outlier.shape = NA)+
  # coord_cartesian(ylim = c(0,4000))+
  ggtitle("Proportion of cases in the HH")+xlab("week")+ylab("Weekly Proportion of Cases")-> prop.HH.cases


png(paste0(figname,'_Figure1.png'),width = 12,height = 7,units = "in",res=300)
p.cases.cum+(p.hh.contacts/p.nhh.contacts) + plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = c('A'),tag_suffix = ')') &
  theme_minimal() & theme(legend.position="right",legend.direction = "vertical")
dev.off()


res.gen %>% gather("Measure","Value",-c('gen','sim','intervention','Intervention')) %>%
  group_by(gen,intervention,Measure) %>%
  summarise(mean=mean(Value),
            median=median(Value)) %>%
  ggplot(aes(x=(gen),color=intervention))+
  geom_line(aes(y=mean,linetype='mean'))+
  geom_line(aes(y=median,linetype='median'))+
  # geom_boxplot(alpha=.5)+
  # geom_boxplot(alpha=.5,outlier.shape = NA)+
  facet_wrap(~Measure,scales = "free")+
  theme_minimal()-> p

png(paste0(figname,'_Figure3.png'),width = 12,height = 7,units = "in",res=300)
p
dev.off()

