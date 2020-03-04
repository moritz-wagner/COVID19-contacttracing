# Quantifying the impact of and effort involved in contact tracing cases using contact data for Kenya

## Introduction

Multiple international cases of the novel coronavirus disease (COVID-19) have been observed outside of China since the start of the epidemic. Whilst the majority of these are exported cases with a direct link to the epidemic in China, some countries are seeing cases of sustained onwards transmissions. Outbreaks in other countries of the magnitude observed in China would pose a significant strain on public health resources, in particular in resource poor settings across Africa. To limit this, early detection and prevention of cases at the start of an outbreak is crucial.

Contact tracing forms one such prevention measure, where close contacts of an infected case are traced. The effectiveness of this depends heavily on the natural history of infection, in particular the proportion of pre- and asymptomatic transmission occurring, as has been [shown recently](https://www.thelancet.com/journals/langlo/article/PIIS2214-109X(20)30074-7/fulltext).
A further limitation is the amount of resources involved in tracing close contacts . Thus, informed decisions must be made of when and how to best implement contact tracing becomes particularly important.

This study makes use of diary-based contact data and includes it into a [previously developed model](https://www.thelancet.com/journals/langlo/article/PIIS2214-109X(20)30074-7/fulltext) of contact tracing. As such contact data typically represent contacts from a single day only, we resample contacts of each participant based on their recorded frequency to create a synthetic set of contacts across multiple days. The participants and their contacts then form our study population into which we introduce initial infected cases and then simulate the early stages of an outbreak. Using this model we explore the effectiveness of different strategies for contact tracing under different transmission scenarios. Additionally, we incorporate several other factors included in the contact data such as age, location (urban, rural), household membership, and frequency of contact.
There are several advantages of including contact data as these allow us to: 
1)	Create a more realistic model, e.g. allowing for household infection clusters
2)	Measure the effort involved in contact tracing, e.g. how many contacts need to be traced per infection at a given time of the outbreak
3)	Include targeted interventions, e.g. isolating household members only

As cases of COVID-19 rise internationally, we hope that this work can inform the Kenyan government and other African nations on how they may be able to effectively implement contact tracing strategies and at which point resources should be focused on other intervention strategies.


## Questions
* IMPACT: How effective is contact tracing, i.e.measured by reduction in Reff
	* What would make it more effective, i.e. delay to isolation, self isolation (of non symptomatic HH contacts), proportion of pre/asymptomatic transmission
* EFFORT: How many contacts need to be traced
* At what point does it become unfeasible in terms of impact and/or effort

## Methods 

### Contacts

#### Sampling from the contact data
For each outbreak simulation a bootstrapped dataset of the original contact data is created. First a set of participants is resampled, then for each of those a set of their contacts are resampled. In this step we also deal with some of the missing data:
* share_HH==TRUE & how_often=="Missing" -> how_often=="Daily"
* how_often=="Missing" & ever_met=="No" -> how_often=="Never"
* how_often=="Missing" that are left are sampled from the rest of the data with the same age of participant, age of contact and location of contact

#### Repeated vs unique contacts in a given time period
We have information on the frequency of a contact and the survey definitions of these, so we can assign a probability p of meeting that contact on a given day:
* Daily -> p=1
* Often 1-2 times per week -> p=(1.5/7)
* Regular 1-2 times per month -> p=(1.5/30)
* Rarely <1 times per month -> p=(0.5/30)
* Never -> p=0

We assume that, given a participant, the contacts for a single day are representative of any other day. So if we have a single contact of a certain frequency on day one, this contact is repeated for the remaining t-1 days. Then based on the frequencies above, this contact has probability p of being repeated. So for each day 1-p gives the probability that a contact is a new contact. If a new contact occurs, the probability of a further new contact on the following days is (1-p)^2, and so on. Random numbers between 0 and 1 are drawn each day to determine whether new contacts have occured. The current model accumulates contacts across a two week period (**sensitivity analysis?**). Weights are added to these contacts according to the expected frequency of contact during that period. Thus daily contacts will have a weight of 14, while new contacts have a weight of 1. Infections are later matched according to these weights during the infection step.

The final output is a full contact data set with participant and contact ids. This defines the susceptible population during the outbreak and for each id (contact and participant), we keep track of who becomes infected.


### Infection and isolation times
This is an extension of the [LSHTM model](https://github.com/epiforecasts/ringbp) to include realistic contact data. For each infector, a number of parameters are sampled:
* incubation period (time)
* delay to isolation (time): This is a measure of how long it takes from symptom onset for a case to be discovered (by seeking health care) and thus be isolated. Note that this does not represent isolation based on contact tracing.
* asymptomatic case (TRUE/FALSE): Each case has a certain probability of being asymptomatic/subclinical. If this is the case, they cannot be isolated. Transmission, however, is still possible
* R0 (integer): The expected number of secondary cases this individual may cause (this will be limited by the type of contacts they have and whether they have already been infected)
* serial interval (time): For each potential new infection (based on R0), this determines their time of exposure
* missed case (TRUE/FALSE): Each case has a certain probability to be missed by tracing in which case they will continue transmitting until symptom onset+delay to isolation

Based on these, for each infection step we determine, which individuals become infected and/or isolated and/or traced. Once a contact becomes isolated, their R0 reduces to 0 and they can no longer infect other individuals. The different scenarios are illustrated in [figure S8](https://www.medrxiv.org/content/medrxiv/suppl/2020/02/11/2020.02.08.20021162.DC1/2020.02.08.20021162-1.pdf) of the LSHTM study.

### Matching infections to contacts
Given an infector with a number of potential infecteds, these now need to be matched to the infector's contacts. We consider two scenarios:
1. The infector is a participant: We know all their contacts immediately.
2. The infector is a contact: Here we don't know the full set of contacts. We know their HH contacts, as these are linked to their participant. The nHH contacts, however, are sampled from participants of the same age group. Note that this may introduce some biases as their may be a relationship between the number of HH contacts and nHH by age.

 We then assign which of the contacts become infected based on their weight. If a contact has already been infected, they won't be infected again, which essentially reduces the reproduction number. We hope that most of these pre-existing infections will be HH contacts, which would create a HH herd immunity effect. If a big proportion of infections are prevented by already infected nHH members this would essentially represent community herd immunity effects, which would be undesirable as we are looking at the early stages of an outbreak with a fully susceptible population. (unless infections cluster in closed communities). **Will check!**

Another factor that might limit the reproduction number is the amount of contacts they have. There will be random draws of R0 that are larger than the total contacts. This is not a problem and if anything more realistic, but **will check how often this happens**.


### Scenarios

* **Isolation upon onset:** Here isolation occurs with a delay upon onset of symptoms, i.e. no tracing of contacts
* **Quarantine:** All traced contacts that are infected become isolated immediately when infector becomes isolated (apart from the ones that are missed)
* **Quarantine plus:** All traced contacts become isolated immediately when infector becomes isolated (apart from the ones that are missed)

For the quarantine scenarios, we distinguish between HH and nHH contacts and set the probability of a HH contact being missed to zero.

## Outputs

* Epidemic curves
* Number of contacts that need to be traced
* R effective
* Number of outbreaks that are controlled: [LSHTM](https://www.medrxiv.org/content/10.1101/2020.02.08.20021162v1.full.pdf) define this as outbreaks where transmission ended within 12 weeks or before 5000 cases in total. Perhaps a similar measure that also takes into account the amount of contacts that require tracing?
* How many contacts need to be traced to succesfully control an outbreak?


## Challenges/Extras
-Repeated contacts: Currently assuming that the contacts of a single day are repeated every day for a given time period. Is this assumption valid?

-Weighting by frequency of contacts: Daily contacts have a higher probability of becoming infected than regular ones. Reasonable?

-Currently drawing times when infection occurs. Should we limit contacts to that duration or consider contacts over a whole 2 week period. Weighting of contacts by frequency probably takes care of this to some extent. Issue: **Which contacts to count for effort of tracing?**

-Higher R0 for HH contacts? [Outside Hubei 70-80% of infections were within HH clusters](https://panopto.lshtm.ac.uk/Panopto/Pages/Viewer.aspx?id=83ba0783-b1ce-4053-aaa5-ab6600da76d8), although might be a by-product of social distancing? -> weighting might take care of this, as daily contacts are most likely to be HH contacts?

-


## Assumptions
* The types of contacts of a single day are representative of any other day and repeated
* The average frequency of a contact for a given period determines the weighting of infection, i.e. more frequent contacts are more likely to become infected
* Contacts over a two week period determine who gets traced and can get infected. Even if contacts become infected in a shorter period. Weighting of more frequent contacts, however, will take care of some potential infection biases. Likely to overestimte 
* Quarantine scenario assumes no delay to isolation, i.e. the time to find contacts is not considered


