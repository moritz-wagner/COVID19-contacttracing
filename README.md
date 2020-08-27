# Using contact data to model the impact of contact tracing and physical distancing to control the SARS-CoV-2 outbreak in Kenya

Source code accompanying Wagner et al. "Using contact data to model the impact of contact tracing and physical distancing to control the SARS-CoV-2 outbreak in Kenya".

## Main files

* The relevant code to rerun the analyses can be found in the [R folder](https://github.com/moritz-wagner/COVID19-contacttracing/tree/master/R)
	* To run different simulations of the model, open COVID19-contacttracing.Rproj as a project in Rstudio and use [R/run_sims.R](https://github.com/moritz-wagner/COVID19-contacttracing/blob/master/R/run_sims.R)
	* To plot different simulations of the model that have been run, open COVID19-contacttracing.Rproj as a project in Rstudio and use [R/plotting.R](https://github.com/moritz-wagner/COVID19-contacttracing/blob/master/R/plotting.R)
* The contact data from [Kiti et al.](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0104786) can be found in the [data folder](https://github.com/moritz-wagner/COVID19-contacttracing/tree/master/data)
* The Supplementary Appendix accompanying the publication can be found under [submission/Supplementary Appendix.pdf](https://github.com/moritz-wagner/COVID19-contacttracing/blob/master/submission/Supplementary%20Appendix.pdf)


## Abstract

**Background** Across the African continent, other than South Africa, COVID-19 cases have remained relatively low. Nevertheless, in Kenya, despite early implementation of containment measures and restrictions, cases have consistently been increasing. Contact tracing forms one of the key strategies in Kenya, but may become infeasible as the caseload grows. Here we explore different contact tracing strategies by distinguishing between household and non-household contacts and how these may be combined with other non-pharmaceutical interventions.

**Methods** We extend a previously developed branching process model for contact tracing to include realistic contact data from Kenya. Using the contact data, we generate a synthetic population of individuals and their contacts categorised by age and household membership. We simulate the initial spread of SARS-CoV-2 through this population and look at the effectiveness of a number of non-pharmaceutical interventions with a particular focus on different contact tracing strategies and the potential effort involved in these.

**Results** General physical distancing and avoiding large group gatherings combined with contact tracing, where all contacts are isolated immediately, can be effective in slowing down the outbreak, but were, under our base assumptions, not enough to control it without implementing extreme stay at home policies. Under optimistic assumptions with a highly overdispersed R0 and a short delay from symptom onset to isolation, control was possible with less stringent physical distancing and by isolating household contacts only.

**Conclusions** Without strong physical distancing measures, controlling the spread of SARS-CoV-2 is difficult. With limited resources, physical distancing combined with the isolation of households of detected cases can form a moderately effective strategy, and control is possible under optimistic assumptions. More data are needed to understand transmission in Kenya, in particular by studying the settings that lead to larger transmission events, which may allow for more targeted responses, and collection of representative age-related contact data.