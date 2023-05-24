# - - - - - - - - - - - - - - - - - - - - - - - 
# Load and process FP data
# - - - - - - - - - - - - - - - - - - - - - - - 

# Set file path

if(Sys.info()["user"]=="akucharski") {
 filepath <- "~/Documents/COVID_data/ILM_data/OneDrive_surveillance/"
} else {
 filepath <- "~/OneDrive - London School of Hygiene and Tropical Medicine/LSHTM_RF/COVID/FrenchPolynesia/"
}

# Load packages
library(qs)
library(readxl)
library(data.table)
library(ggplot2)
library(lubridate)
library(ISOweek)

source("R/utils.R")
source("R/date.R")

# Load data files ---------------------------------------------------------

# Population
pop <- qread(paste0(filepath,"unwpp_data.qs")) 

# Weekly count data
weekly_dt <- as.data.table(read_xlsx(paste0(filepath,"Weekly data summary.xlsx")))

# Cases
cases_dt <- fread(paste0(filepath,"FP_processed_May_30/data cases.csv"))

# Hopsitalisations in 1st and 2nd wave
hosps_dt1 <- fread(paste0(filepath,"Hospi from  09-08-20 tot 22-06-2021 (wave 1).csv"))
hosps_dt2 <- fread(paste0(filepath,"Hospi from 23-06-2021 to nov 21 (wave 2).csv"))
hosps_dt3 <- fread(paste0(filepath,"FP_processed_May_30/data hospitalisations omicron.csv"))

# Serology Feb 21
sero_pos_dt1 <- as.data.table(read_xlsx(paste0(filepath,"Serop_feb21.xlsx")))

# Serology Nov 21
sero_pos_dt2 <- fread(paste0(filepath,"Serop_nov21.csv"))

# Vaccinations
# vax <- fread(paste0(filepath,"SPC_DF_COVID_VACCINATION_1.0_D.PF.COVIDVACAD1+COVIDVACAD2+COVIDVACBST+COVIDVACADT.csv"))
# dir_vax <- "FP_processed_May_30/"
dir_vax <- "2022_06_08_Vaccination_all_FP/"
files <- list.files(paste0(filepath,dir_vax))
files_vax <- files[grep("daily vaccination",files)]
age_groups_vax <- gsub("[a-z]+| |\\.","",files_vax)
vax_list <- vector("list", length(files_vax))
for (i in 1:length(files_vax)){
    vax_list1 <- as.data.table(read_xlsx(paste0(filepath,dir_vax,files_vax[i])),
                               sheet = 1)
    vax_list1[,`Type d'injection V2` := "Primo injection"]
    vax_list2 <- as.data.table(read_xlsx(paste0(filepath,dir_vax,files_vax[i]),
                                         sheet = 2))
    vax_list[[i]] <- rbind(vax_list1,vax_list2)
}
vax <- rbindlist(vax_list,idcol = "file")
vax[,age_group := age_groups_vax[file]]
vax[,file := NULL]

# Process data ---------------------------------------------------------
# Age groups
age_groups <- c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70+")
min_ages <- get_min_age(age_groups)

# Population
pop <- pop[country == "French Polynesia" & year == 2020]
pop[,age_group := cut(age,c(min_ages,Inf),labels = age_groups,right = F)]
write.csv(pop,"data/population.csv",row.names = F)

## Weekly data
weekly_dt[,Date := sub(".*-","",Date)]
weekly_dt[,Year := as.integer(paste0("20",sub("[0-9]+/[0-9]+/","",Date)))]
weekly_dt[Date == "03/01",Year := 2021L]
weekly_dt[Date == "02/01",Year := 2022L]
setnafill(weekly_dt,type = "locf",cols = "Year")
weekly_dt[,Date := paste0(sub("([0-9]+/[0-9]+).*","\\1",Date),"/",Year)]
weekly_dt[,Date := dmy(Date) - 6]

## Cases
# Drop empty columns
cases_dt[,Column1 := NULL]

names(cases_dt) <- c("case_number","age_group","date","test_result_date",
                     "vaccination_complete","vaccine_type","last_dose_date",
                     "hospital","death_date","history_of_covid","archipelago")

# Convert Excel date numbers to dates
cases_dt[,death_date := as.IDate(death_date,origin = "1899-12-30")]

# Correct typos in entry date
cases_dt[case_number == 6729, date := as.IDate("2020-10-21")] # date incorrectly given as 2020-01-21
cases_dt[case_number == 73214, date := as.IDate("2022-03-10")] # date incorrectly given as 2020-03-10

# Clean age group variable
cases_dt[,age_group := sub("_","-",age_group)]
cases_dt[,age_group := sub("#NUM!|#VALUE!","",age_group)]

# Change age groups to model age groups
cases_dt[,min_age := get_min_age(age_group)]
cases_dt[,age_group := cut(min_age,c(min_ages,Inf),labels = age_groups,right = F)]

# Plot total cases against weekly data to check
# Total cases
total_cases_dt <- cases_dt[!is.na(date),.(cases = .N),by = .(date)]
total_cases_dt[,iso_week := ISOweek(date)]
total_cases_dt <- total_cases_dt[,.(cases = sum(cases,na.rm = T)),by = .(iso_week)]
total_cases_dt[,date := ISOweek2date(paste0(iso_week,"-1"))]

ggplot() + 
    geom_line(aes(x = date,y = cases),total_cases_dt) +
    geom_point(aes(x = Date,y = `Total nombre de nouveaux cas confirmés locaux`),weekly_dt,color = "red")
ggsave("output/total_cases.pdf",width = 5,height = 4)

# Impute missing confirmation dates with dates of nearest cases
cases_dt1 <- copy(cases_dt) 
setnafill(cases_dt1, type = "locf", cols = "date")
total_cases_dt1 <- cases_dt1[,.(cases = .N),by = .(date)]
total_cases_dt1[,iso_week := ISOweek(date)]
total_cases_dt1 <- total_cases_dt1[,.(cases = sum(cases,na.rm = T)),by = .(iso_week)]
total_cases_dt1[,date := ISOweek2date(paste0(iso_week,"-1"))]

ggplot() + 
    geom_line(aes(x = date,y = cases),total_cases_dt1) +
    geom_point(aes(x = Date,y = `Total nombre de nouveaux cas confirmés locaux`),weekly_dt,color = "red")
ggsave("output/total_cases_imputed_missing_dates.pdf",width = 5,height = 4)
# Most cases with missing dates are early in first wave, so use data with imputed missing dates

# Aggregate cases by age group and date
cases <- cases_dt[!is.na(date) & !is.na(age_group),.(cases = .N),by = .(age_group,date)]

# Plot cases
ggplot(cases,aes(x = date,y = cases,group = age_group,color = age_group)) + 
    geom_line() #+ 
# facet_wrap(~age_group)

# Aggregate cases with imputed missing dates by age group and date
cases1 <- cases_dt1[!is.na(age_group),.(cases = .N),by = .(age_group,date)]

# Plot cases
ggplot(cases1,aes(x = date,y = cases,group = age_group,color = age_group)) + 
    geom_line()

## Hospitalisations
# 1st wave 
names(hosps_dt1) <- c("date","age","sex","death","death_date")
hosps_dt1[,`:=`(sex = fifelse(sex == "Masculin",1,0),
               death = fifelse(death == "Oui",1,0))]

# 2nd wave
# Drop empty columns
hosps_dt2[,(names(hosps_dt2)[ncol(hosps_dt2)-(0:2)]):=NULL]

names(hosps_dt2) <- c("case_num","age","sex","collection_date","vaccination_complete",
                      "vaccine_type","last_dose_date","hospital","death_date","history_of_covid")

hosps_dt2[,`:=`(age = as.integer(sub(" ans","",age)),
                sex = fifelse(sex == "M",1,0),
                collection_date = as.IDate(collection_date, format = "%d/%m/%Y"),
                vaccination_complete = fifelse(vaccination_complete == "Oui",1,0),
                last_dose_date = as.IDate(last_dose_date, format = "%d/%m/%Y"),
                death_date = as.IDate(death_date, format = "%d/%m/%Y"),
                history_of_covid = fifelse(history_of_covid == "oui",1,0))]

# Omicron hospitalisations
# Drop empty columns
hosps_dt3[,c("Column1","Column2","Column3") := NULL]

names(hosps_dt3) <- c("date","age","vaccination_complete","vaccine_type",
                      "hospital","death_date","history_of_covid")

hosps_dt3[,`:=`(age = as.integer(sub("_.*","",age)),
                vaccination_complete = fifelse(vaccination_complete == "Oui",1,0))]

# Hospitalisations Nov-Dec 2021
hosps_dt_Nov21 <- cases_dt[
    hospital != "" & between(date,
                             hosps_dt2[,max(collection_date)],
                             hosps_dt3[,min(date)],incbounds = F),
    !c("test_result_date","archipelago")]
hosps_dt_Nov21[,age := as.integer(sub("_.*","",age_group))]


# Bind data frames
hosps_dt <- rbind(hosps_dt1,hosps_dt2,hosps_dt_Nov21,hosps_dt3,fill=T)

# FOR NOW: Assume hospitalisation date and sample collection date are the same
hosps_dt[,date := fifelse(is.na(collection_date),date,collection_date)]

# Plot total hospitalisations against weekly data to check
# Total hospitalisations by hospital
total_hosps_by_hosp_dt <- hosps_dt[,.(hosps = .N),by = .(date,hospital)]
total_hosps_by_hosp_dt[,iso_week := ISOweek(date)]
total_hosps_by_hosp_dt <- total_hosps_by_hosp_dt[,.(hosps = sum(hosps)),by = .(iso_week,hospital)]
total_hosps_by_hosp_dt[,date := ISOweek2date(paste0(iso_week,"-1"))]

# CHPF (main hospital on Tahiti)
ggplot() + 
    geom_line(aes(x = date,y = hosps),total_hosps_by_hosp_dt[is.na(hospital) | hospital == "CHPF"]) +
    geom_point(aes(x = Date,y = `Nombre de nouvelles hospitalisations Covid CHPf Tahiti`),weekly_dt,color = "red") + 
    labs(title = "CHPF")
ggsave("output/total_hosps_CHPF.pdf",width = 5,height = 4)

# hospitals <- total_hosps_by_hosp_dt[,unique(hospital)]
# hospitals <- hospitals[!is.na(hospitals) & hospitals != "CHPF"]
# for (i in seq_along(hospitals)){
#     y <- sym(paste0("Nombre de nouvelles hospitalisations Covid ",hospitals[i]))
#     print(ggplot() + 
#               geom_line(aes(x = date,y = hosps),total_hosps_by_hosp_dt[hospital == hospitals[i]]) + 
#               geom_point(aes(x = Date,y = !!y),weekly_dt,color = "red") + 
#               labs(title = hospitals[i]))    
# }

# Total hospitalisations in all hospitals
total_hosps_dt <- hosps_dt[,.(hosps = .N),by = .(date)]
total_hosps_dt[,iso_week := ISOweek(date)]
total_hosps_dt <- total_hosps_dt[,.(hosps = sum(hosps)),by = .(iso_week)]
total_hosps_dt[,date := ISOweek2date(paste0(iso_week,"-1"))]

ggplot() + 
    geom_line(aes(x = date,y = hosps),total_hosps_dt) +
    geom_point(aes(x = Date,y = `Nombre total de nouvelles hospitalisations tous hôpitaux`),weekly_dt,color = "red")
ggsave("output/total_hosps.pdf",width = 5,height = 4)

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

# Correct typo in date of death
deaths_dt[date == as.IDate("2021-01-14") & death_date == as.IDate("2021-01-06"),
          death_date := as.IDate("2021-01-16")]

# Aggregate over age groups
deaths <- deaths_dt[,.(deaths=.N),by=.(age_group,date = death_date)]

# Plot deaths
ggplot(deaths,aes(x = date,y = deaths,group = age_group,color = age_group)) + geom_line()

# Plot total deaths against weekly data to check
total_deaths_by_hosp_dt <- deaths_dt[,.(deaths = .N),by = .(date = death_date,hospital)]
total_deaths_by_hosp_dt[,iso_week := ISOweek(date)]
total_deaths_by_hosp_dt <- total_deaths_by_hosp_dt[,.(deaths = sum(deaths)),by = .(iso_week,hospital)]
total_deaths_by_hosp_dt[,date := ISOweek2date(paste0(iso_week,"-1"))]

ggplot() + 
    geom_line(aes(x = date,y = deaths),total_deaths_by_hosp_dt[is.na(hospital) | hospital == "CHPF"]) +
    geom_point(aes(x = Date,y = `Nombre de décès CHPf`),weekly_dt,color = "red") + 
    labs(title = "CHPF")
ggsave("output/total_deaths_CHPF.pdf",width = 5,height = 4)

total_deaths_dt <- deaths_dt[,.(deaths = .N),by = .(date = death_date)]
total_deaths_dt[,iso_week := ISOweek(date)]
total_deaths_dt <- total_deaths_dt[,.(deaths = sum(deaths)),by = .(iso_week)]
total_deaths_dt[,date := ISOweek2date(paste0(iso_week,"-1"))]

ggplot() + 
    geom_line(aes(x = date,y = deaths),total_deaths_dt) + 
    geom_point(aes(x = Date,y = `Nombre total de décès`),weekly_dt,color = "blue") + 
    geom_point(aes(x = Date,y = `Nombre total de décès hospitaliers`),weekly_dt,color = "red")
ggsave("output/total_deaths.pdf",width = 5, height = 4)

## Seroprevalence
process_sero_data <- function(sero_dt,sample_date,age_groups,min_ages){
    x <- copy(sero_dt)
    
    setnames(x,"age_group","age_group_sero")
    
    x[,date := sample_date]
    
    max_ages <- x[,get_max_age(age_group_sero)]
    max_ages[is.na(max_ages)] <- Inf
    
    x[,age_group := cut(max_ages,c(min_ages,Inf),labels = age_groups)]
    x <- x[,lapply(.SD, function(y) sum(as.integer(round(y)))),.SDcols = c("n","seropos"),by = .(date,age_group)]
    x[,seroprev := mapply(
        function(y,z) paste0(round(100*y/z,1)," (", 
                             round(100*binom.test(y,z)$conf.int[1],1),"--",
                             round(100*binom.test(y,z)$conf.int[2],1),")"),
        seropos,n
    )]
    
    return(x)
}

sero_pos_dt1 <- process_sero_data(sero_pos_dt1, as.Date("2021-02-14"), age_groups, min_ages)
write.csv(sero_pos_dt1,"data/seroprev_feb21.csv",row.names = F)
print(binom.test(sero_pos_dt1[,sum(seropos)],sero_pos_dt1[,sum(n)]))
sero_pos_dt2 <- process_sero_data(sero_pos_dt2, as.Date("2021-11-30"), age_groups, min_ages)
write.csv(sero_pos_dt2,"data/seroprev_nov21.csv",row.names = F)
print(binom.test(sero_pos_dt2[,sum(seropos)],sero_pos_dt2[,sum(n)]))

sero_pos_dt <- rbind(sero_pos_dt1,sero_pos_dt2)

## Variant sequencing
variant_dt <- as.data.table(read_xlsx(paste0(filepath,"FP_processed_May_30/variant_screening.xlsx"),skip = 1))
cols_to_keep <- c("Date","Week","ALPHA...3","DELTA...4","BA1...5","BA2...6","GAMMA...7","MU...8")
variant_dt <- variant_dt[!is.na(Date),..cols_to_keep]
setnames(variant_dt,cols_to_keep,tolower(sub("...[0-9]","",cols_to_keep)))
variant_dt[,date := as.Date(date)]

# Exclude values before 2022
variant_dt <- variant_dt[date >= as.Date("2022-01-01")]
variant_dt <- variant_dt[, `:=`(strain_tot = ba1 + ba2, strain_non_variant = ba1)]


## Make data table of hospitalisations, deaths, cases and seroprevalence for fitting 
strt_date <- hosps_dt[,min(date,na.rm = T)] - 20 # 2020-07-20
end_date <- as.Date("2022-05-06") # last death date in data files
dates <- seq.Date(strt_date,end_date,by = 1)
base_dt <- CJ(date = dates,age_group = age_groups_hosp)

reformat_data <- function(x, base_dt, vrbl, fillna = F){
    max_date <- x[,max(date)]
    x <- merge(base_dt,x,by = c("date","age_group"),all.x = T)
    if(fillna){
        x[date <= max_date, (vrbl) := nafill(get(vrbl), fill = 0)]
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

base_dt_case <- CJ(date = dates, age_group = age_groups)
cases_wide <- reformat_data(cases1,base_dt_case,"cases",T)

setnames(sero_pos_dt,c("n","seropos"),c("sero_tot_1","sero_pos_1"))
base_dt_sero <- CJ(date = dates,age_group = age_groups[3:length(age_groups)])
sero_pos_wide <- reformat_data(sero_pos_dt,base_dt_sero,"sero_pos_1")
sero_tot_wide <- reformat_data(sero_pos_dt,base_dt_sero,"sero_tot_1")

base_dt_vrnt <- data.table(date = dates)
vrnt <- merge(base_dt_vrnt,variant_dt[,.(date,strain_tot,strain_non_variant)],by = "date",all.x = T)

# Merge different data sources
data_raw <- Reduce(function(...) merge(...,all = T), list(hosps_wide, deaths_wide, cases_wide, sero_pos_wide, sero_tot_wide, vrnt))

data_raw[,day := covid_multi_strain_date(date)]
data_raw[,date := NULL]

# Add empty columns for total hospitalisations and deaths
# data_raw$hosps <- NA
data_raw[,hosps := hosps_0_39 + hosps_40_49 + hosps_50_59 + hosps_60_69 + hosps_70_plus]
# data_raw$deaths <- NA
data_raw[,deaths := deaths_0_39 + deaths_40_49 + deaths_50_59 + deaths_60_69 + deaths_70_plus]
data_raw[,cases := cases_0_9 + cases_10_19 + cases_20_29 + cases_30_39 + cases_40_49 + cases_50_59 + cases_60_69 + cases_70_plus]
data_raw[,sero_pos_1 := sero_pos_1_20_29 + sero_pos_1_30_39 + sero_pos_1_40_49 + sero_pos_1_50_59 + sero_pos_1_60_69 + sero_pos_1_70_plus]
data_raw[,sero_tot_1 := sero_tot_1_20_29 + sero_tot_1_30_39 + sero_tot_1_40_49 + sero_tot_1_50_59 + sero_tot_1_60_69 + sero_tot_1_70_plus]
# data_raw[,strain_tot := NA]
# data_raw[,strain_non_variant := NA]
write.csv(data_raw,"data/data_cases_hosps_deaths_serology.csv",row.names = F)

## Vaccinations

# Change names
setnames(vax,c("Date","Nombre de doses injectées","Type d'injection V2"),c("date","number","dose"))

# Convert date column to Date type
vax[,date := as.Date(date)]

# Recode dose variable
vax[,dose := fcase(dose == "Primo injection","dose1",
                   dose == "Schéma vaccinal complet","dose2",
                   dose == "Rappel","dose3")]

# Aggregate doses in the same age group on the same day
vax <- vax[,.(number = sum(number)),by = .(date,dose,age_group)]
write.csv(vax,"data/data_vaccination.csv",row.names = F)

# Probability of death in the community given severe disease
prob_death_community <- weekly_dt[
    ,sum(`Nombre de décès à domicile`,na.rm = T)/
        sum(`Nombre total de nouvelles hospitalisations tous hôpitaux`,na.rm = T)]
write.csv(prob_death_community,"data/prob_death_community.csv",row.names = F)
