# Load the required libraries 
library(Boruta)
library(DMwR2)
library(Hmisc)
library(mice)
library(caret)

# Load the required datasets


# Contacts dataset 
## Meaning of different columns in the dataset 
# cont_id : global identifying number of contact 
# part_id : global identfying number of the diary 
# cnt_age_exact : Exact age of the contact person reported 
# cnt_age_est_min : Left limit of age range of reported contact 
# cnt_age_est_max : Right limit of age range of reported contact 
# cnt_gender : Gender of the person (M/F)
# cnt_home : Contact at home 
# cnt_work : Contact at work 
# cnt_school : Contact at school 
# cnt_transport : Contact at transport 
# cnt_leisure : Contact during leisure activities 
# cnt_otherplace : Contact during other activities 
# frequency_multi : Frequency of contacting the reported contact 
# (1:daily, 2:weekly, 3:monthly, 4: a few times a year, 5: first time)
# phys_contact: 1 -> yes and 0 -> no 
# duration_multi: 1: < 5 mins, 2-> 5-15 mins, 3-> 15mins-1hr, 4-> 1-4hours, 5-> >4hours 

contact_common <- read.csv("2008_Mossong_POLYMOD_contact_common.csv")
head(contact_common)


# Participant Dataset 
# part_id: global identifying number of diary
# hh_id: household identifying number 
# part_age: age of the person whom the diary refers to (adult or child)
# part_gender: gender of the participant to whom the diary belongs to (F: Female, M: Male)
participant_common <- read.csv("2008_Mossong_POLYMOD_participant_common.csv")
head(participant_common)


# Participant Extra 
participant_extra <- read.csv("2008_Mossong_POLYMOD_participant_extra.csv")
head(participant_extra)

# Household common dataset 
hh_common <- read.csv("2008_Mossong_POLYMOD_hh_common.csv")


# Match household common with participant common 
household_lookup <- merge(participant_common[c("part_id","hh_id")], hh_common)

# Merge the tables 
df_participant_temp <- merge(participant_common,participant_extra,by = "part_id")
df <- merge(df_participant_temp,contact_common, by = "part_id")
df <- merge(df,household_lookup, by = c("part_id","hh_id"))

# Calculate estimated age by taking the average of cnt_age_est_min and cnt_age_est_max 
avg_est_age <- ifelse(!is.na(df[["cnt_age_est_min"]]) & !is.na(df[["cnt_age_est_max"]]), (df[["cnt_age_est_min"]] + df[["cnt_age_est_max"]]) / 2,
                      ifelse(is.na(df[["cnt_age_est_min"]]) & !is.na(df[["cnt_age_est_max"]]), df[["cnt_age_est_max"]],
                             ifelse(!is.na(df[["cnt_age_est_min"]]) & is.na(df[["cnt_age_est_max"]]), df[["cnt_age_est_min"]],
                                    NA)))


df[["contact_age"]] <- ifelse(is.na(df[["cnt_age_exact"]]),avg_est_age,df[["cnt_age_exact"]])
df[["contact_age"]][is.na(df[["contact_age"]])] <- mean(df[["contact_age"]])

# Subset out dataset of interest 
df <- df[c("part_id","part_gender","contact_age","part_age","country","hh_size","cnt_gender",
           "cnt_home","cnt_work","cnt_school","cnt_transport","cnt_leisure","cnt_otherplace")]


# Create separate binary columns for male and female 
df[["gender_female"]] <- ifelse(df[["part_gender"]] == "M",0,1)
df[["gender_male"]] <- ifelse(df[["part_gender"]] == "M",1,0)

# Remove the cnt_gender column 
df <- df[, -which(names(df) == "cnt_gender")]
df <- df[, -which(names(df) == "part_gender")]


# Create binary columns for countries 
df[["country_BE"]] <- ifelse(df[["country"]] == "BE",1,0)
df[["country_DE"]] <- ifelse(df[["country"]] == "DE",1,0)
df[["country_FI"]] <- ifelse(df[["country"]] == "FI",1,0)
df[["country_GB"]] <- ifelse(df[["country"]] == "GB",1,0)
df[["country_IT"]] <- ifelse(df[["country"]] == "IT",1,0)
df[["country_LU"]] <- ifelse(df[["country"]] == "LU",1,0)
df[["country_NL"]] <- ifelse(df[["country"]] == "NL",1,0)
df[["country_PL"]] <- ifelse(df[["country"]] == "PL",1,0)

df <- df[, -which(names(df) == "country")]


names(df) = c('participant_id',
        'contact_age',
        'participant_age',
        'household_size',
        'contact_home',
        'contact_work',
        'contact_school',
        'contact_transport',
        'contact_leisure',
        'contact_other',
        'gender_female',
        'gender_male',
        'country_BE',
        'country_DE',
        'country_FI',
        'country_GB',
        'country_IT',
        'country_LU',
        'country_NL',
        'country_PL')

df = na.omit(df)





