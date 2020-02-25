##Script to get new unique contacts given id and time t
library(readr)
library(tidyverse)
library(ggplot2)
library(data.table)

source('R/ct_functions.R')

##Get and tidy contact data
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


##Time t and id
t <- 14
cd <- contact_boot(contact_data)
id <- sample(cd$csid,1)
count.contacts(id,cd,t)
