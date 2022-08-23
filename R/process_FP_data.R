# - - - - - - - - - - - - - - - - - - - - - - - 
# Load and process FP data
# - - - - - - - - - - - - - - - - - - - - - - - 

# Set file path

if(Sys.info()["user"]=="akucharski") {
 filepath <- "~/Documents/COVID_data/ILM_data/OneDrive_surveillance/"
} else {
 filepath <- "~/OneDrive - London School of Hygiene and Tropical Medicine/LSHTM_RF/COVID/FrenchPolynesia/"
}

# Load data files ---------------------------------------------------------

# Population
pop <- qread(paste0(filepath,"unwpp_data.qs")) 

# Contact matrix
contact_matrices <- fread(paste0(filepath,"synthetic_contacts_2020.csv")) 
setnames(contact_matrices,"age_cotactee","age_contactee")

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
pop <- pop[country == "French Polynesia" & year == 2020]
pop[,age_group := cut(age,c(min_ages,Inf),labels = age_groups,right = F)]
agg_pop <- pop[,.(population = sum(total)),by = .(age_group)]
population <- agg_pop[,population]

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
pop_contact <- pop_contact[,.(population = sum(total)),by = .(age_group_contact)]

contact1 <- merge(contact1,pop_contact,by.x = "age_contactor",by.y = "age_group_contact")
contact1[,age_contactor := cut(min_age_contactor,c(min_ages,Inf),labels = age_groups,right = F)]

# Take population-weighted average of mean number of contacts over "participant" age groups being aggregated
contact2 <- contact1[,.(contacts = sum(contacts * population)/sum(population)), by = .(age_contactor,age_contactee)]

# Plot
ggplot(contact2,aes(x = age_contactee,y = age_contactor,fill = contacts)) + geom_tile()

# Convert the contact matrix to the "transmission matrix" (the contact matrix weighted by the population in each age group)
contact_matrix <- matrix(contact2[,contacts],nrow = n_age,ncol = n_age, byrow = T) # N.B. individuals in rows, contacts in columns
transmission <- contact_matrix/rep(population, each = ncol(contact_matrix))

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

# Aggregate cases by age group and date
cases <- cases_dt[!is.na(date) & age_group != "",.(cases = .N),by = .(age_group,date)]

# Plot cases
ggplot(cases,aes(x = date,y = cases,group = age_group,color = age_group)) + 
    geom_line() #+ 
    # facet_wrap(~age_group)

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
# ggsave("output/total_hosps_CHPF.pdf",width = 5,height = 4)

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
# ggsave("output/total_hosps.pdf",width = 5,height = 4)

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
# ggsave("output/total_deaths_CHPF.pdf",width = 5,height = 4)

total_deaths_dt <- deaths_dt[,.(deaths = .N),by = .(date = death_date)]
total_deaths_dt[,iso_week := ISOweek(date)]
total_deaths_dt <- total_deaths_dt[,.(deaths = sum(deaths)),by = .(iso_week)]
total_deaths_dt[,date := ISOweek2date(paste0(iso_week,"-1"))]

ggplot() + 
    geom_line(aes(x = date,y = deaths),total_deaths_dt) + 
    geom_point(aes(x = Date,y = `Nombre total de décès`),weekly_dt,color = "blue") + 
    geom_point(aes(x = Date,y = `Nombre total de décès hospitaliers`),weekly_dt,color = "red")
# ggsave("output/total_deaths.pdf",width = 5, height = 4)

## Seroprevalence
process_sero_data <- function(sero_dt,sample_date,age_groups,min_ages){
    x <- copy(sero_dt)
    
    setnames(x,"age_group","age_group_sero")
    
    x[,date := sample_date]
    
    max_ages <- x[,get_max_age(age_group_sero)]
    max_ages[is.na(max_ages)] <- Inf
    
    x[,age_group := cut(max_ages,c(min_ages,Inf),labels = age_groups)]
    x <- x[,lapply(.SD, function(y) sum(as.integer(round(y)))),.SDcols = c("n","seropos"),by = .(date,age_group)]
    
    return(x)
}

sero_pos_dt1 <- process_sero_data(sero_pos_dt1, as.Date("2021-02-14"), age_groups, min_ages)
sero_pos_dt2 <- process_sero_data(sero_pos_dt2, as.Date("2021-11-30"), age_groups, min_ages)

sero_pos_dt <- rbind(sero_pos_dt1,sero_pos_dt2)

## Make data table of hospitalisations, deaths, cases and seroprevalence for fitting 
strt_date <- hosps_dt[,min(date,na.rm = T)] - 20 # 2020-07-20
end_date <- as.Date("2022-05-23") # last date in vaccination data files
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

base_dt_case <- CJ(date = dates, age_group = age_groups)
cases_wide <- reformat_data(cases,base_dt_case,"cases",T)

setnames(sero_pos_dt,c("n","seropos"),c("sero_tot_1","sero_pos_1"))
base_dt_sero <- CJ(date = dates,age_group = age_groups[3:length(age_groups)])
sero_pos_wide <- reformat_data(sero_pos_dt,base_dt_sero,"sero_pos_1")
sero_tot_wide <- reformat_data(sero_pos_dt,base_dt_sero,"sero_tot_1")

# Merge different data sources
data_raw <- Reduce(function(...) merge(...,all = T), list(hosps_wide, deaths_wide, cases_wide, sero_pos_wide, sero_tot_wide))

data_raw[,day := as.integer(date - min(date) + 1L)]
data_raw[,date := NULL]

# Add empty columns for total hospitalisations and deaths
# data_raw$hosps <- NA
data_raw[,hosps := hosps_0_39 + hosps_40_49 + hosps_50_59 + hosps_60_69 + hosps_70_plus]
# data_raw$deaths <- NA
data_raw[,deaths := deaths_0_39 + deaths_40_49 + deaths_50_59 + deaths_60_69 + deaths_70_plus]
data_raw[,cases := cases_0_9 + cases_10_19 + cases_20_29 + cases_30_39 + cases_40_49 + cases_50_59 + cases_60_69 + cases_70_plus]
data_raw[,sero_pos_1 := sero_pos_1_20_29 + sero_pos_1_30_39 + sero_pos_1_40_49 + sero_pos_1_50_59 + sero_pos_1_60_69 + sero_pos_1_70_plus]
data_raw[,sero_tot_1 := sero_tot_1_20_29 + sero_tot_1_30_39 + sero_tot_1_40_49 + sero_tot_1_50_59 + sero_tot_1_60_69 + sero_tot_1_70_plus]
data_raw[,strain_tot := NA]
data_raw[,strain_non_variant := NA]

## Vaccinations

# # Remove bits before : in column names
# names(vax) <- sub(".*: ","",names(vax))
# 
# cols <- names(vax)[!(names(vax) %in% c("Time","OBS_VALUE"))]
# vax[,(cols) := lapply(.SD,function(x) sub(".*: ","",x)),.SDcols = cols]
# 
# # Drop columns
# vax <- vax[,c("DATAFLOW","Frequency","Unit of measure","Unit multiplier","Observation Status","Data source","Comment"):=NULL]
# 
# # Rename columns
# names(vax) <- c("date","country","dose","value")
# 
# # Drop total doses and booster doses
# vax <- vax[!(dose %in% c("Total doses administered","Booster doses administered"))]
# # Recode dose variable
# vax[,dose := fifelse(dose=="1st dose administered",1L,2L)]
# 
# # Plot
# ggplot(vax,aes(x = date,y = value,group = factor(dose),color = factor(dose))) + geom_line()
# 
# # Make data table of daily vaccination doses from raw data
# dates_vax <- seq.Date(vax[,min(date)-7],vax[,max(date)],by = 1)
# 
# vax_dt <- CJ(dose = c(1L,2L), date = dates_vax)
# vax_dt <- merge(vax_dt,vax[,!"country"],by = c("dose","date"),all.x = T)
# 
# # Linearly interpolate cumulative number of doses over missing dates
# vax_dt[date == min(date),value := 0]
# vax_dt[,value_interp := approx(date,value,date)$y,by = .(dose)]
# 
# # Calculate approximate daily numbers of doses by differencing and rounding
# vax_dt[,number := as.integer(round(diff(c(0,value_interp)))),by = .(dose)]
# 
# # Difference in rounded doses vs actual
# vax_dt[,sum(number),by = .(dose)][,V1] - vax_dt[date == max(date),value,by = .(dose)][,value]
# # small so ignore FOR NOW
# 
# vax_dt[,date := as.IDate(date)]
# vax_dt <- vax_dt[date <= end_date]
# 
# doses <- vax[,unique(dose)]

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

# # # Plot to check
# # ggplot(vax,aes(x = date,y = number,group = age_group,color = age_group)) +
# #     geom_line() +
# #     facet_wrap(~dose)
# # ggplot(vax[,.(date,number = cumsum(number)),by = .(dose,age_group)],aes(x = date,y = number,group = age_group,color = age_group)) +
# #     geom_line() +
# #     facet_wrap(~dose)
# doses_by_age_and_dose <- dcast(vax[,.(number = sum(number)),by=.(age_group,dose)],age_group ~ dose)
# doses_by_age_and_dose <- rbind(doses_by_age_and_dose,cbind(data.table(age_group="Total"),doses_by_age_and_dose[,lapply(.SD,sum),.SDcols = c("dose1","dose2","dose3")]))
# write.csv(doses_by_age_and_dose,"output/doses_by_age_and_dose.csv",row.names = F)

# Add different delays for immune response to different doses
delay_dose1 <- 28
delay_dose2 <- 14
vax[dose == "dose1", date := date + delay_dose1]
vax[dose == "dose2", date := date + delay_dose2]

# Reaggregate vaccine doses by model age groups
dates_vax <- seq.Date(vax[,min(date)],vax[,max(date)],by = 1)
doses <- vax[,unique(dose)]
base_vax_dt <- CJ(date = dates_vax, dose = doses, age = pop[,unique(age)])
base_vax_dt <- merge(base_vax_dt,pop[,.(age,total)],by = "age")
age_groups_vax <- sort(age_groups_vax)
min_ages_vax <- get_min_age(age_groups_vax)
base_vax_dt[, age_group := cut(age,c(min_ages_vax,Inf),labels = age_groups_vax,right = F)]
# Merge with vaccinations data table
# N.B. This duplicates doses across age groups
vax_dt <- merge(base_vax_dt,vax,by = c("date","dose","age_group"),all.x = T)
# Split vaccine doses by population proportion
vax_dt[,number := number*total/sum(total),by = .(date,dose,age_group)]
# Change age groups
vax_dt[,age_group := cut(age,c(min_ages,Inf),labels = age_groups,right = F)]
# Fill missing values with 0s (i.e. assume all doses were recorded)
setnafill(vax_dt,fill = 0,cols = "number")
# Sum doses over age groups
vax_dt <- vax_dt[,.(number = sum(number)),by = .(date,dose,age_group)]
# Limit vaccine schedule to end date
vax_dt <- vax_dt[date <= end_date]

# Plot to check
ggplot(vax_dt[age_group!="0-9"],aes(x = date,y = number,group = age_group,color = age_group)) +
    geom_line() +
    facet_wrap(~dose)

# # Check totals are the same
# print(vax[,sum(number)]) # 465247
# print(vax_dt[,sum(number)]) # 465247

# Cast to wide format
vax_dt_wide <- dcast(vax_dt,date + age_group ~ dose,value.var = "number")
vax_dt_wide[,age_band_min := get_min_age(age_group)]
vax_dt_wide[,age_group := NULL]

# # Make vaccination schedule
# pop_mat <- matrix(rep(population,1),nrow = length(population))
# 
# # Rough number of daily booster doses from eyeballing Fig. 5 plot in BEH health bulletin No. 84
# booster_daily_doses <- c(rep(0,40*7),rep(5000/7,as.integer(as.IDate(end_date) - vax_dt[,min(date)] - 40*7 + 1)))
# 
# # priority_population <- vaccine_priority_population(population, uptake = c(0.55,0.55,0.80,0.80,0.80,0.80,0.8,0.8))
# mean_days_between_doses <- 28 # from eye-balling plot of cumulative 1st and 2nd doses
# schedule <- vaccine_schedule_future(as.integer(vax_dt[,min(date)] - strt_date + 1),
#                                     vax_dt[,sum(number),by = .(date)][,V1],
#                                     mean_days_between_doses = mean_days_between_doses,
#                                     pop_mat,
#                                     booster_daily_doses_value = booster_daily_doses)

# Matrix of uptake rates (age group x dose)
uptake <- matrix(1,nrow = length(age_groups),ncol = length(doses))
schedule <- vaccine_schedule_from_data(as.data.frame(vax_dt_wide),min_ages,population,uptake)

# Plot to check
doses_dt <- as.data.table(schedule$doses,value.name = "number")
doses_dt[,`:=`(age_group = age_groups[V1],dose = doses[V2], date = dates_vax[V3])]
doses_dt <- merge(doses_dt,agg_pop,by = "age_group")
doses_dt[,prop := number/population]
doses_dt[,cum_prop := cumsum(prop),by = .(age_group,dose)]

ggplot(doses_dt,aes(x = date,y = cum_prop,group = age_group, color = age_group)) + 
    geom_line() + 
    xlim(as.Date("2021-01-01"),NA_Date_) + 
    labs(x = "Date", y = "Proportion vaccinated") + 
    scale_color_discrete(name = "Age group") + 
    facet_wrap(~dose, 
               labeller = labeller(
                   dose = c("dose1" = "Dose 1", 
                            "dose2" = "Dose 2", 
                            "dose3" = "Dose 3")))
ggsave("output/vax_cov_by_dose.pdf",width = 9,height = 3.5)
