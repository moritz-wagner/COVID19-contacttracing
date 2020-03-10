require(magrittr)
require(patchwork)

## Run outbreak sims and do some plotting
source('R/contact_functions.R')
source('R/outbreak_functions.R')

##Plot of parameters used-------------
n <- 10000
#Incubation periods
incfn <- dist_setup(dist_shape = 2.322737, dist_scale = 6.492272)
inc.period <- incfn(n)
mean.inc <- mean(inc.period)


#Serial intervals based on incubation samples
k = 0.7
inc.period.sample <- sample(inc.period,1)
# inc.period.sample <- max(inc.period)
serial.interval <- inf_fn(rep(inc.period.sample,n),k)
mean.serial.interval <- mean(serial.interval)
perc.pre.symp <- sum(serial.interval<inc.period.sample)/n


#Delay to isolation
delay_shape = 2.5
delay_scale = 5
delayfn <- dist_setup(dist_shape = delay_shape,dist_scale = delay_scale)
iso.delay <- delayfn(n)+inc.period.sample
mean.iso.delay <- mean(iso.delay)


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
  geom_vline(aes(xintercept=mean.inc))+
  xlim(0,max(iso.delay))+ggtitle('Incubation Period') -> p1
data.frame(n=1:n,days=serial.interval) %>% 
  ggplot(aes(x=days))+
  geom_density(fill="blue",alpha=.5)+
  geom_vline(aes(xintercept=mean.serial.interval))+
  geom_vline(aes(xintercept=inc.period.sample),color="red")+
  xlim(0,max(iso.delay))+ggtitle('Serial interval')+
  annotate("text", 
           y= .2,label=paste0(round(perc.pre.symp*100),"% presymp"), x =max(iso.delay),hjust=1)-> p2
data.frame(n=1:n,days=iso.delay) %>% 
  ggplot(aes(x=days))+
  geom_density(fill="green",alpha=.5)+
  geom_vline(aes(xintercept=mean.iso.delay))+
  geom_vline(aes(xintercept=inc.period.sample),color="red")+
  xlim(0,max(iso.delay))+ggtitle('Delay to isolation') -> p3
data.frame(n=1:n,sec.infs=r0) %>% 
  count(sec.infs) %>% mutate(density=n/sum(n)) %>% 
  ggplot(aes(x=(sec.infs),y=density))+
  geom_bar(stat="identity")+
  coord_cartesian(xlim=c(0,20))+
  # geom_histogram(fill="black",alpha=.5,aes(y = ..density..))+
  geom_vline(aes(xintercept=mean.r0))+ggtitle('R0')+xlab("") -> p4

png('plots/params.png',width = 7,height = 5,units = "in",res=300)
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
time.CD <- system.time(CD <- count.contacts.full(contact_data,t=14))

sims <- 500
cap_cases <- 2000
cap_max_days <- 8*7

# Isolation only-----
quarantine <- FALSE
HH_trace <- FALSE
prop.ascertain <- 0

time.ext.par<-system.time(res.ext <- scenario_sim_ext_parallel(n.sim = sims,num.initial.cases = 5,cap_max_days = cap_max_days,
                                                               cap_cases = cap_cases,r0isolated = 0,r0community = 2.5,
                                                               disp.iso = 1,disp.com = 0.16,k = 0.7,delay_shape = 2.5,
                                                               delay_scale = 5,prop.asym = 0,prop.ascertain = prop.ascertain,
                                                               quarantine = quarantine,HH_trace = HH_trace,
                                                               contact_data = CD))

res.no.tracing <- res.ext$res.week
res.no.tracing$intervention <- "Isolation only"

# Tracing HH-----
quarantine <- FALSE
HH_trace <- TRUE
prop.ascertain <- 0


time.ext.par<-system.time(res.ext <- scenario_sim_ext_parallel(n.sim = sims,num.initial.cases = 5,cap_max_days = cap_max_days,
                                                               cap_cases = cap_cases,r0isolated = 0,r0community = 2.5,
                                                               disp.iso = 1,disp.com = 0.16,k = 0.7,delay_shape = 2.5,
                                                               delay_scale = 5,prop.asym = 0,prop.ascertain = prop.ascertain,
                                                               quarantine = quarantine,HH_trace = HH_trace,
                                                               contact_data = CD))
res.tracing.HH <- res.ext$res.week
res.tracing.HH$intervention <- "Tracing HH"

# Tracing HH + 25% non-HH-----
quarantine <- FALSE
HH_trace <- TRUE
prop.ascertain <- 0.25

time.ext.par<-system.time(res.ext <- scenario_sim_ext_parallel(n.sim = sims,num.initial.cases = 5,cap_max_days = cap_max_days,
                                                               cap_cases = cap_cases,r0isolated = 0,r0community = 2.5,
                                                               disp.iso = 1,disp.com = 0.16,k = 0.7,delay_shape = 2.5,
                                                               delay_scale = 5,prop.asym = 0,prop.ascertain = prop.ascertain,
                                                               quarantine = quarantine,HH_trace = HH_trace,
                                                               contact_data = CD))

res.tracing.HH2 <- res.ext$res.week
res.tracing.HH2$intervention <- "Tracing HH + 25% non-HH"

# Quaranting HH-----
quarantine <- TRUE
HH_trace <- TRUE
prop.ascertain <- 0

time.ext.par<-system.time(res.ext <- scenario_sim_ext_parallel(n.sim = sims,num.initial.cases = 5,cap_max_days = cap_max_days,
                                                               cap_cases = cap_cases,r0isolated = 0,r0community = 2.5,
                                                               disp.iso = 1,disp.com = 0.16,k = 0.7,delay_shape = 2.5,
                                                               delay_scale = 5,prop.asym = 0,prop.ascertain = prop.ascertain,
                                                               quarantine = quarantine,HH_trace = HH_trace,
                                                               contact_data = CD))

res.quarantine <- res.ext$res.week
res.quarantine$intervention <- "Quarantine HH"


## put the together and plot-----

res <- bind_rows(list(res.no.tracing,res.tracing.HH,res.tracing.HH2,res.quarantine))
res <- bind_rows(list(res.no.tracing,res.tracing.HH,res.tracing.HH2))

res %>% group_by(intervention) %>% 
  subset(week==max(week)) %>% 
  count(cumulative<cap_cases & weekly_cases==0) %>% 
  mutate(prop=n/sum(n))

res %>% 
  ggplot(aes(x=factor(week),y=cumulative,fill=intervention,color=intervention))+
  geom_boxplot(alpha=.5,outlier.shape = NA)+
  # coord_cartesian(ylim = c(0,2000))+
  # stat_summary(aes(group=intervention),fun.y = mean,geom = "line")
  ggtitle("Cumulative Cases")+xlab("week")+ylab("Cumulative Cases") -> p.cases

res %>% 
  mutate(hh.con=ifelse(intervention=='Isolation only',0,hh.con)) %>% 
  # subset(intervention!="Isolation only") %>% 
  ggplot(aes(x=factor(week),y=hh.con,fill=intervention,color=intervention))+
  geom_boxplot(alpha=.5,outlier.shape = NA)+
  # coord_cartesian(ylim = c(0,2000))+
  ggtitle("HH contacts to trace")+xlab("week")+ylab("Contacts") -> p.hh.contacts

res %>% 
  mutate(nhh.con=ifelse(intervention%in%c('Isolation only',"Tracing HH"),0,nhh.con*.25)) %>% 
  # subset(intervention!="Isolation only") %>% 
  ggplot(aes(x=factor(week),y=nhh.con,fill=intervention,color=intervention))+
  geom_boxplot(alpha=.5,outlier.shape = NA)+
  # coord_cartesian(ylim = c(0,4000))+
  ggtitle("Non-HH contacts to trace")+xlab("week")+ylab("Contacts")-> p.nhh.contacts


require(patchwork)
png('plots/example_plot.png',width = 10,height = 5,units = "in",res=300)
p.cases+(p.hh.contacts/p.nhh.contacts) + plot_layout(guides = 'collect') & theme_minimal()
dev.off()

