# Impact of and effort involved in contact tracing cases using contact data

## Questions
* IMPACT: How effective is contact tracing, i.e.measured by reduction in Reff
	-What would make it more effective, i.e. delay to isolation, self isolation (of non symptomatic HH contacts), proportion of pre/asymptomatic transmission
-EFFORT: How many contacts need to be traced given n infecteds
-At what point does it become unfeasible in terms of impact and/or effort

## Methods 
### Contacts (Moritz)

#### Sampling from the contact data
For each outbreak simulation a bootstrapped dataset of the original contact data is created. First a set of participants is resampled, then for each of those a set of their contacts are resampled. In this step we also deal with some of the missing data:
* share_HH==TRUE & how_often=="Missing" -> how_often=="Daily"
* how_often=="Missing" & ever_met=="No" -> how_often=="Never"
* how_often=="Missing" that are left are sampled from the rest of the data based on age of participant and contact

#### Repeated vs unique contacts in a given time period
We have information on the frequency of a contact and the survey definitions of these, so we can assign a probability p of meeting that contact on a given day:
* Daily -> p=1
* Often 1-2 times per week -> p=(1.5/7)
* Regular 1-2 times per month -> p=(1.5/30)
* Rarely <1 times per month -> p=(0.5/30)
* Never -> p=0

We assume that, given a participant, the contacts for a single day are representative of any other day. So if we have a single contact of a certain frequency on day one, this contact is repeated for the remaining t-1 days. Then based on the frequencies above, this contact has probability p of being repeated. So for each day 1-p gives the probability that a contact is a new contact. If a new contact occurs, the probability of a further new contact on the following days is (1-p)^2, and so on. Random numbers between 0 and 1 are drawn each day to determine whether new contacts have occured. We also keep track of how often a contact is repeated in order to give an estimate for the contact intensity and thus the risk of transmission over a given period.


-For a given time period and a single infected we resample from the contact data to produce the number of unique new HH and non-HH contacts
-We keep track of HH infecteds to limit the susceptible pool of future infections


### Infection and isolation times (Ivy)
-Extending the LSHTM model to include contact data and the depletion of contacts due to HH clusters
-If for each infected we know who they contact, we can assign who becomes infected and who is successfully isolated based on sampling from probability distributions of R0, incubation periods, serial intervals, proportion of asymptomatic transmission, delay to and probability of isolation.

### Challenges
-Repeated contacts: Currently assuming that the contacts of a single day are repeated every day for a given time period. Is this assumption valid? Then using probabilities based on the recorded frequency of the contacts, the number of contacts that are new on a given day are calculated, e.g. Daily contacts will not be repeated, while regular contacts have a given probability of being repeated vs being new every day. (see figure)

-Weighting by frequency of contacts: Daily contacts have a higher probability of becoming infected than regular ones. Reasonable?

-Single generation step vs full outbreak: If we want to model a full outbreak, we have to keep track of who is in which HH. This is difficult and messy, as we don’t have HH ids. We have thought of ways around it, but it requires a number of assumptions. An alternative would be to model a single generation only, i.e. one infection step, and look at how effective isolation is in that scenario. It would give a lower bound to the impact of isolation, as any secondary transmissions are ignored.
