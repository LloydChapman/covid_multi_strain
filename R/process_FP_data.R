library(data.table)
library(ggplot2)
library(qs)
library(readxl)

source("R/vaccination.R")

get_min_age = function(x){
    as.numeric(sub("-.*","",sub("\\+|<","-",x)))  
}

get_max_age = function(x){
    as.numeric(sub(".*-","",sub("\\+|<","-",x)))    
}

# Age groups
age_groups <- c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70+")
min_ages <- get_min_age(age_groups)
n_age <- length(age_groups)

# Population
pop <- qread("~/OneDrive - London School of Hygiene and Tropical Medicine/LSHTM_RF/COVID/comix/data/unwpp_data.qs")
pop <- pop[country == "French Polynesia" & year == 2020]
pop[,age_group := cut(age,c(min_ages,Inf),labels = age_groups,right = F)]
agg_pop <- pop[,.(population = sum(total)),by = .(age_group)]
population <- agg_pop[,population]

# Contact matrix
contact_matrices <- fread("~/OneDrive - London School of Hygiene and Tropical Medicine/LSHTM_RF/COVID/FrenchPolynesia/synthetic_contacts_2020.csv")
setnames(contact_matrices,"age_cotactee","age_contactee")

# FOR NOW: use contact matrix for France (is Fiji an alternative as a Pacific island?)
contact <- contact_matrices[iso3c == "FRA" & setting == "overall" & location_contact == "all"]
age_cols <- names(contact)[grep("age",names(contact))]
contact[,(paste0("min_",age_cols)):=lapply(.SD,function(x) as.numeric(sub("\\+","",sub(" to.*","",x)))),.SDcols = age_cols] 

# Sum mean numbers of contacts over contact age groups being aggregated
contact[,age_contactee := cut(min_age_contactee,c(min_ages,Inf),labels = age_groups,right = F)]
contact1 <- contact[,.(contacts = sum(mean_number_of_contacts)),by = .(age_contactor,age_contactee,min_age_contactor)]

min_ages_contact <- contact[,unique(min_age_contactor)]
age_groups_contact <- contact[,unique(age_contactor)]
pop_contact <- copy(pop)
pop_contact[,age_group_contact := cut(age,c(min_ages_contact,Inf),labels = age_groups_contact,right = F)]
pop_contact <- pop_contact[,.(population = sum(population)),by = .(age_group_contact)]

contact1 <- merge(contact1,pop_contact,by.x = "age_contactor",by.y = "age_group_contact")
contact1[,age_contactor := cut(min_age_contactor,c(min_ages,Inf),labels = age_groups,right = F)]

# Take population-weighted average of mean number of contacts over "participant" age groups being aggregated
contact2 <- contact1[,.(contacts = sum(contacts * population)/sum(population)), by = .(age_contactor,age_contactee)]

# Plot
ggplot(contact2,aes(x = age_contactee,y = age_contactor,fill = contacts)) + geom_tile()

# Convert the contact matrix to the "transmission matrix" (the contact matrix weighted by the population in each age group)
contact_matrix <- matrix(contact2[,contacts],nrow = n_age,ncol = n_age, byrow = T) # N.B. individuals in rows, contacts in columns
transmission <- contact_matrix/rep(population, each = ncol(contact_matrix))

## Hospitalisations
# 1st wave
hosps_dt1 <- fread("~/OneDrive - London School of Hygiene and Tropical Medicine/LSHTM_RF/COVID/FrenchPolynesia/Hospi from  09-08-20 tot 22-06-2021 (wave 1).csv")
names(hosps_dt1) <- c("entry_date","age","sex","death","death_date")

hosps_dt1[,`:=`(sex = fifelse(sex == "Masculin",1,0),
               death = fifelse(death == "Oui",1,0))]


# 2nd wave
hosps_dt2 <- fread("~/OneDrive - London School of Hygiene and Tropical Medicine/LSHTM_RF/COVID/FrenchPolynesia/Hospi from 23-06-2021 to nov 21 (wave 2).csv")

# Drop empty columns
hosps_dt2[,(names(hosps_dt2)[ncol(hosps_dt2)-(0:2)]):=NULL]

names(hosps_dt2) <- c("case_num","age","sex","collection_date","vaccination_complete",
                      "vaccine_type","last_dose_date","hospital","death_date","history_of_covid")

hosps_dt2[,`:=`(age = as.integer(sub(" ans","",age)),
                sex = fifelse(sex == "M",1,0),
                collection_date = as.IDate(collection_date, format = "%d/%m/%Y"),
                vaccination_complete = fifelse(vaccination_complete == "Oui",1,0),
                last_dose_date = as.Date(last_dose_date, format = "%d/%m/%Y"),
                death_date = as.IDate(death_date, format = "%d/%m/%Y"),
                history_of_covid = fifelse(history_of_covid == "oui",1,0))]

# Bind data frames
hosps_dt <- rbind(hosps_dt1,hosps_dt2,fill=T)

# FOR NOW: Assume hospitalisation date and sample collection date are the same
hosps_dt[,date := fifelse(is.na(collection_date),entry_date,collection_date)]

# Add age groups
age_groups_hosp <- c(paste0(min_ages[c(1,5:(length(min_ages)-1))],"-",min_ages[5:length(min_ages)]-1), paste0(min_ages[length(min_ages)],"+"))
min_ages_hosp <- get_min_age(age_groups_hosp)
hosps_dt[,age_group := cut(age,c(min_ages_hosp,Inf), labels = age_groups_hosp, right = F)]
hosps_dt[,sum(is.na(age_group))]
# only 4 hospitalisations missing ages 
# FOR NOW: drop
hosps_dt <- hosps_dt[!is.na(age_group)]

# Aggregate over age groups
hosps <- hosps_dt[,.(hosps = .N),by = .(age_group,date)]

# Plot hospitalisations
ggplot(hosps,aes(x = date,y = hosps,group = age_group,color = age_group)) + geom_line()

## Deaths
deaths_dt <- hosps_dt[!is.na(death_date)]
deaths <- deaths_dt[,.(deaths=.N),by=.(age_group,date)]

# Plot deaths
ggplot(deaths,aes(x = date,y = deaths,group = age_group,color = age_group)) + geom_line()

## Seroprevalence
sero_pos_dt <- as.data.table(read_xlsx("~/OneDrive - London School of Hygiene and Tropical Medicine/LSHTM_RF/COVID/FrenchPolynesia/Serop_feb21.xlsx"))
setnames(sero_pos_dt,"age_group","age_group_sero")

sero_pos_dt[,date := as.Date("2021-02-14")] # FOR NOW: use middle of February 2021 as sample collection was done in February 2021

max_ages_sero <- sero_pos_dt[,get_max_age(age_group_sero)]
max_ages_sero[is.na(max_ages_sero)] <- Inf

sero_pos_dt[,age_group := cut(max_ages_sero,c(min_ages,Inf),labels = age_groups)]
sero_pos_dt <- sero_pos_dt[,lapply(.SD,function(x) sum(as.integer(round(x)))),.SDcols = c("n","seropos"),by = .(date,age_group)]

## Make data table of hospitalisations, deaths and seroprevalence for fitting 
strt_date <- hosps_dt[,min(date)] - 15
end_date <- as.Date("2021-11-21")
dates <- seq.Date(strt_date,end_date,by = 1)
base_dt <- CJ(date = dates,age_group = age_groups_hosp)

reformat_data <- function(x, base_dt, vrbl, fillna = F){
    x <- merge(base_dt,x,by = c("date","age_group"),all.x = T)
    if(fillna){
        setnafill(x,fill = 0,cols = vrbl)    
    }
    x_wide <- dcast(x, date ~ age_group, value.var = vrbl)
    idx <- 2:ncol(x_wide)
    names(x_wide)[idx] <- paste0(vrbl,"_",sub("-","_",sub("\\+","_plus",names(x_wide)[idx])))
    # x_wide[,day := as.integer(date - min(date) + 1L)]
    return(x_wide)
}

# Reformat data to required format for fitting
hosps_wide <- reformat_data(hosps,base_dt,"hosps",T)
deaths_wide <- reformat_data(deaths,base_dt,"deaths",T)

setnames(sero_pos_dt,c("n","seropos"),c("sero_tot_1","sero_pos_1"))
base_dt_sero <- CJ(date = dates,age_group = age_groups[3:length(age_groups)])
sero_pos_wide <- reformat_data(sero_pos_dt,base_dt_sero,"sero_pos_1")
sero_tot_wide <- reformat_data(sero_pos_dt,base_dt_sero,"sero_tot_1")

# Merge different data sources
data_raw <- Reduce(function(...) merge(...,all = T), list(hosps_wide, deaths_wide, sero_pos_wide, sero_tot_wide))

data_raw[,day := as.integer(date - min(date) + 1L)]
data_raw[,date := NULL]

# Add empty columns for total hospitalisations and deaths
data_raw$hosps <- NA
data_raw$deaths <- NA
data_raw$sero_pos_1 <- NA
data_raw$sero_tot_1 <- NA
data_raw$cases <- NA
data_raw$cases_non_variant <- NA

## Vaccinations
vax <- fread("~/OneDrive - London School of Hygiene and Tropical Medicine/LSHTM_RF/COVID/FrenchPolynesia/SPC_DF_COVID_VACCINATION_1.0_D.PF.COVIDVACAD1+COVIDVACAD2+COVIDVACADT.csv",)

# Remove bits before : in column names
names(vax) <- sub(".*: ","",names(vax))

cols <- names(vax)[!(names(vax) %in% c("Time","OBS_VALUE"))]
vax[,(cols) := lapply(.SD,function(x) sub(".*: ","",x)),.SDcols=cols]

# Drop columns
vax <- vax[,c("DATAFLOW","Frequency","Unit of measure","Unit multiplier","Observation Status","Data source","Comment"):=NULL]

# Rename columns
names(vax) <- c("date","country","dose","value")

# Drop total doses
vax <- vax[dose!="Total doses administered"]
# Recode dose variable
vax[,dose := fifelse(dose=="1st dose administered",1L,2L)]

# Plot
ggplot(vax,aes(x = date,y = value,group = factor(dose),color = factor(dose))) + geom_line()

# Make data table of daily vaccination doses from raw data
dates <- seq.Date(vax[,min(date)-7],vax[,max(date)],by = 1)

vax_dt <- CJ(dose = c(1L,2L), date = dates)
vax_dt <- merge(vax_dt,vax[,!"country"],by = c("dose","date"),all.x = T)

# Linearly interpolate cumulative number of doses over missing dates
vax_dt[date == min(date),value := 0]
vax_dt[,value_interp := approx(date,value,date)$y,by = .(dose)]

# Calculate approximate daily numbers of doses by differencing and rounding
vax_dt[,number := as.integer(round(diff(c(0,value_interp)))),by = .(dose)]

# Difference in rounded doses vs actual
vax_dt[,sum(number),by = .(dose)][,V1] - vax_dt[date == max(date),value,by = .(dose)][,value]
# small so ignore FOR NOW

vax_dt[,date := as.IDate(date)]

# Make vaccination schedule
pop_mat <- matrix(rep(population,1),nrow = length(population))

# Rough number of daily booster doses from eyeballing Fig. 5 plot in BEH health bulletin
booster_daily_doses <- c(rep(0,40*7),rep(5000/7,125))

# priority_population <- vaccine_priority_population(population, uptake = c(0.55,0.55,0.80,0.80,0.80,0.80,0.8,0.8))
schedule <- vaccine_schedule_future(as.integer(vax_dt[,min(date)] - hosps_dt[,min(date)]),
                                    vax_dt[,sum(number),by = .(date)][,V1],
                                    mean_days_between_doses = 28, # from eye-balling plot of cumulative 1st and 2nd doses
                                    pop_mat,
                                    booster_daily_doses_value = booster_daily_doses)

# Plot to check
doses_1 <- as.data.table(reshape2::melt(schedule$doses[,1, ],value.name = "number"))
doses_1[,`:=`(age_group = age_groups[Var1],date = dates[Var2])]
doses_1 <- merge(doses_1,agg_pop,by = "age_group")
doses_1[,prop := number/population]
doses_1[,cum_prop := cumsum(prop),by = .(age_group)]

ggplot(doses_1,aes(x = date,y = cum_prop,group = age_group, color = age_group)) + 
    geom_line()

doses_3 <- as.data.table(reshape2::melt(schedule$doses[,3, ],value.name = "number"))
doses_3[,`:=`(age_group = age_groups[Var1],date = dates[Var2])]
doses_3 <- merge(doses_3,agg_pop,by = "age_group")
doses_3[,prop := number/population]
doses_3[,cum_prop := cumsum(prop),by = .(age_group)]

ggplot(doses_3,aes(x = date,y = cum_prop,group = age_group, color = age_group)) + 
    geom_line()



