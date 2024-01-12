
# TOPIC: FACTORS ASSOCIATED WITH THE RISK OF DEATH FROM ASTHMA AND IF 
# THEY HAVE CHANGED AS A RESULT OF INFECTION BY COVID, 
# A PROSPECTIVE COHORT STUDY OF ASTHMA PATIENTS (2019 To 2021)
#-------------------------------------------------------------------------------
# Install and run the packages below.
install.packages("dplyr")
install.packages("survival")
install.packages("magrittr")
install.packages("tidyverse")
install.packages("caret")


library(dplyr)
library(magrittr)
library(survival)
library(tidyverse)
library(lubridate)
library(tcltk)
library(RODBC)
library(ROCR)
#-------------------------------------------------------------------------------
## Set working directory and to assess SAIL Databank  use the query below

setwd(".../2210245_202 Assignment/202 Assignment/2210245_Assignment A2/
      
      210245_PMIM202 Assignment Annotated Script A1")

source("login.R")
#------------------------------------------------------------------------------
# A. Run Query to import Asthma and Covid-19 patients and demographics
#    from SAIL1281V project Databank using PEDW_SINGLE_DIAG_20220228 and 
#    PEDW_SINGLE_DIAG_20220228 .
#-------------------------------------------------------------------------------
## Query  for patient with Asthma ICD-10 codes and Read codes from SAIL and join

a <- sqlQuery(channel, "SELECT ALF_PE, min(START_DATE) as start_date FROM 

                        SAIL1281V.PEDW_SINGLE_DIAG_20220228 WHERE 
                            
                        DIAG_CD_1234 LIKE 'J45%' OR DIAG_CD_1234 LIKE 'H33%'
                            
                        GROUP BY ALF_PE")

b <- sqlQuery(channel, "SELECT DISTINCT (ALF_PE ), min(EVENT_DT ) as 
                        
                        start_date_b FROM 
                        
                        SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 WHERE 
                        
                        EVENT_CD LIKE 'H33%' OR EVENT_CD LIKE '14BA%' 
              
                        GROUP BY ALF_PE")

asthma_patients <- a %>% full_join(b, by = 'ALF_PE')

asthma_patients <- asthma_patients %>%
  
                  mutate(start_date = pmin(START_DATE, 
                  
                  START_DATE_B, na.rm = TRUE)) %>% 
    
                  select(ALF_PE, start_date) %>% na.omit()

#------------------------------------------------------------------------------
# B. Table for COVID-19 Patients using PEDW_SINGLE_DIAG_20220228 and
#    WLGP_GP_EVENT_CLEANSED_20220201
#------------------------------------------------------------------------------
## Query  for patient with COVID ICD-10 and Read codes from SAIL and join data

c <- sqlQuery(channel, "SELECT ALF_PE, 
                           
                        min(START_DATE) as start_date  FROM 
                            
                        SAIL1281V.PEDW_SINGLE_DIAG_20220228 WHERE 
                          
                        DIAG_CD_1234 LIKE 'U071%' OR DIAG_CD_1234 LIKE 
                          
                        'U072%' OR DIAG_CD_1234 LIKE 'R060%'GROUP BY ALF_PE")


d <- sqlQuery(channel, "SELECT DISTINCT (ALF_PE ), 
          
                       min(EVENT_DT ) as start_date_b FROM 
    
                       SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 WHERE 
                       
                       EVENT_CD LIKE '4J3R1%' OR EVENT_CD LIKE '1JX1.%' 
              
                       GROUP BY ALF_PE")


covid_patients <- c %>% full_join(d, by = 'ALF_PE')

covid_patients <- covid_patients %>%
  
                 mutate(covid_st_date = pmin(START_DATE, 
                                             
                 START_DATE_B, na.rm = TRUE)) %>% 
  
                 select(ALF_PE, covid_st_date) %>% na.omit()



#------------------------------------------------------------------------------
# C. Table for patients demographics from WLGP_PATIENT_ALF_CLEANSED_20220201
#------------------------------------------------------------------------------

# Query  for patient demographics from SAIL and convert GNDR_CD
# to gender.

patient_demographic <- sqlQuery(channel, "SELECT ALF_PE, WOB, GNDR_CD FROM 
                                
                                SAIL1281V.WLGP_PATIENT_ALF_CLEANSED_20220201")

patient_demographic <- patient_demographic %>% group_by(ALF_PE) %>% 
                            
                                              filter(row_number() == 1)

patient_demographic <- patient_demographic %>%  
        
                      mutate(GNDR_CD=ifelse(GNDR_CD==1, 0, 1)) %>% 
  
                      rename(gender=GNDR_CD)

#------------------------------------------------------------------------------
## D. Create joins for tables imported, restrict study period from 2019 to 2022
#     and filter for analysis.
#-------------------------------------------------------------------------------

# Create a Join for asthma and Covid tables and rename columns.

asthma_covid_combined <- inner_join(asthma_patients, covid_patients,
                                    
                                                           by = 'ALF_PE')

colnames(asthma_covid_combined)[2] <- 'asthma_onset'

colnames(asthma_covid_combined)[3] <- 'covid_onset'


## Join the created asthma and covid table to patient demographics.

asthma_covid_combined1 <- left_join(asthma_covid_combined, patient_demographic, 
                                      
                                                                 by = 'ALF_PE')


## Derive the age at diagnosis using WOB and date of asthma diagnosis.

asthma_covid_combined1 <- asthma_covid_combined1 %>% mutate(YOB = year(WOB)) %>% 
                   
                           mutate(asthma_onset_year = year(asthma_onset)) %>% 
                   
                           mutate(age_at_diag = asthma_onset_year - YOB) %>%  
                    
                           select(-c(YOB, asthma_onset_year))


## Remove NA values from the combined table and view the summary.

a_c_c <- asthma_covid_combined1 %>% filter(!is.na(WOB))

a_c_c <- na.omit(a_c_c)

summary(a_c_c)

# Create a variable to remove ages with negative values from combined table

negative_age <- a_c_c$age_at_diag < 0

a_c_c <- a_c_c[!negative_age, ]


# Restrict study period for asthma and covid patients to start from 2019-01-01 
# for patients diagnosed with COVID

a_c_c_pop <- a_c_c %>% filter(covid_onset>"2018-12-31") %>% 
  
                                            arrange(ALF_PE, covid_onset)


#------------------------------------------------------------------------------
# E. Create a table for Control Cases for Asthma table through anti-join, rename
#    column to asthma-onset for first diagnosis and derive age at diagnosis.
#-----------------------------------------------------------------------------

# Create an anti-join of asthma patients and combined asthma and covid patients

control_cases <- asthma_patients %>% anti_join(asthma_covid_combined1, 
                                           
                                                         by = 'ALF_PE')

colnames(control_cases)[2] <- 'asthma_onset'

control_cases1 <- left_join(control_cases, patient_demographic,
                            
                                                          by = 'ALF_PE')


control_cases1 <- control_cases1 %>% mutate(YOB = year(WOB)) %>% 
  
                  mutate(asthma_onset_year = year(asthma_onset)) %>% 
  
                  mutate(age_at_diag = asthma_onset_year - YOB) %>%  
   
                  select(-c(YOB, asthma_onset_year))


# Remove NAs from the derived control population

f_control_cases <- control_cases1 %>% filter(!is.na(WOB))

f_control_cases <- na.omit(f_control_cases)

summary(f_control_cases)

# Create a variable to remove ages with negative values from control population 

neg_age <- f_control_cases$age_at_diag < 0

f_control_cases <- f_control_cases[!neg_age, ]


# Select the the sample size for control population and exclude population
# diagnosed with asthma before start date of the study

f_control_cases_pop <- f_control_cases 

#-------------------------------------------------------------------------------
## F. Important Risk factors for the hypothesis from PEDW_SINGLE_DIAG_20220228,
#     EDDS_EDDS_20220301 and, WLGP_GP_EVENT_CLEANSED_20220201 and 
#     merge with control and study population for hypothesis analysis.
#-------------------------------------------------------------------------------

# RISK FACTORS FOR HYPOTHESIS
# Select ALFE_PE with risk factors and create the risk factor variable

obs <- sqlQuery(channel, "SELECT ALF_PE FROM(SELECT ALF_PE FROM

                    SAIL1281V.PEDW_SINGLE_DIAG_20220228 WHERE DIAG_CD_1234 
                   
                    LIKE 'E66%' UNION  SELECT ALF_PE FROM 
                   
                    SAIL1281V.EDDS_EDDS_20220301 WHERE DIAG_CD_1 LIKE 'E66%'
                   
                    UNION SELECT ALF_PE FROM 
                   
                    SAIL1281V.PEDW_SINGLE_DIAG_20220228 WHERE DIAG_CD_1234 =
                   
                    'C380' UNION  SELECT ALF_PE FROM 
                   
                    SAIL1281V.EDDS_EDDS_20220301 WHERE DIAG_CD_1 = 'C380') 
                   
                    AS RESULT GROUP BY ALF_PE HAVING COUNT(*) = 1")

obs1 <-  sqlQuery(channel, "SELECT DISTINCT (ALF_PE ) FROM

                    SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 WHERE 
                    
                    EVENT_CD LIKE '22A%' OR EVENT_CD LIKE '22K%'
                  
                    GROUP BY ALF_PE")

obesity <- obs %>% full_join(obs1, by = 'ALF_PE')

obesity$obesity <- 1



sap <- sqlQuery(channel, "SELECT ALF_PE FROM(SELECT ALF_PE FROM 
                    
                    SAIL1281V.PEDW_SINGLE_DIAG_20220228 WHERE DIAG_CD_1234 
                    
                    LIKE 'G47%' UNION  SELECT ALF_PE FROM 
                    
                    SAIL1281V.EDDS_EDDS_20220301 WHERE DIAG_CD_1 LIKE 'G47%'
                    
                    UNION SELECT ALF_PE FROM 
                    
                    SAIL1281V.PEDW_SINGLE_DIAG_20220228 WHERE DIAG_CD_1234 =
                    
                    'R068' UNION  SELECT ALF_PE FROM 
                    
                    SAIL1281V.EDDS_EDDS_20220301 WHERE DIAG_CD_1 = 'R068') 
                    
                    AS RESULT GROUP BY ALF_PE HAVING COUNT(*) = 1")


sap1 <- sqlQuery(channel, "SELECT DISTINCT (ALF_PE ) FROM 
                    
                    SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 WHERE EVENT_CD 
                    
                    LIKE 'Fy03%' OR EVENT_CD LIKE 'H5B%' GROUP BY ALF_PE")

sleep_apnea <- sap %>% full_join(sap1, by = 'ALF_PE')

sleep_apnea$sleep_apnea <- 1


smk <- sqlQuery(channel, "SELECT ALF_PE FROM(SELECT ALF_PE FROM 
                    
                    SAIL1281V.PEDW_SINGLE_DIAG_20220228 WHERE DIAG_CD_1234 
                    
                    IN ('F172', 'Z878', '1375') UNION  SELECT ALF_PE FROM 
                    
                    SAIL1281V.EDDS_EDDS_20220301 WHERE DIAG_CD_1  IN 
                    
                    ('F172', 'Z878', '1375') UNION SELECT ALF_PE FROM 
                    
                    SAIL1281V.PEDW_SINGLE_DIAG_20220228 WHERE DIAG_CD_1234 IN
                    
                    ('Z270', 'T652', '1373') UNION  SELECT ALF_PE FROM 
                    
                    SAIL1281V.EDDS_EDDS_20220301 WHERE DIAG_CD_1 IN
                    
                    ('Z270', 'T652', '1373')) 
                    
                    AS RESULT GROUP BY ALF_PE HAVING COUNT(*) = 1")

smk1 <- sqlQuery(channel, "SELECT DISTINCT (ALF_PE ) FROM 
               
                    SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 WHERE EVENT_CD
                    
                    LIKE '137P%' OR EVENT_CD LIKE 'E251%' GROUP BY ALF_PE")

smoking <- smk %>% full_join(smk1, by = 'ALF_PE')

smoking$smoking <- 1


alco <- sqlQuery(channel, "SELECT ALF_PE FROM(SELECT ALF_PE FROM 
                    
                    SAIL1281V.PEDW_SINGLE_DIAG_20220228 WHERE DIAG_CD_1234 
                    
                    IN ('Z721', 'F100', 'Z502', 'Z714') UNION  SELECT ALF_PE 
                    
                    FROM SAIL1281V.EDDS_EDDS_20220301 WHERE DIAG_CD_1  IN 
                    
                    ('Z721', 'F100', 'Z502', 'Z714') UNION SELECT ALF_PE FROM 
                    
                    SAIL1281V.PEDW_SINGLE_DIAG_20220228 WHERE DIAG_CD_1234 IN
                    
                    ('1282', '1366', '1462') UNION  SELECT ALF_PE FROM 
                    
                    SAIL1281V.EDDS_EDDS_20220301 WHERE DIAG_CD_1 IN
                    
                    ('1282', '1366', '1462')) 
                    
                    AS RESULT GROUP BY ALF_PE HAVING COUNT(*) = 1")


alco1 <- sqlQuery(channel, "SELECT DISTINCT (ALF_PE ) FROM 
     
                    SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 WHERE EVENT_CD 
                    
                    LIKE '136..%' GROUP BY ALF_PE")


alcohol <- alco %>% full_join(alco1, by = 'ALF_PE')

alcohol$alcohol <- 1



fam <- sqlQuery(channel, "SELECT ALF_PE FROM(SELECT ALF_PE FROM 
                    
                    SAIL1281V.PEDW_SINGLE_DIAG_20220228 WHERE DIAG_CD_1234 = 
                    
                    'Z825' UNION  SELECT ALF_PE FROM 
                    
                    SAIL1281V.EDDS_EDDS_20220301 WHERE DIAG_CD_1 = 'Z825' 
                    
                    UNION SELECT ALF_PE FROM 
                    
                    SAIL1281V.PEDW_SINGLE_DIAG_20220228 WHERE DIAG_CD_1234 =
                    
                    'Z836' UNION  SELECT ALF_PE FROM 
                    
                    SAIL1281V.EDDS_EDDS_20220301 WHERE DIAG_CD_1  =
                    
                    'Z836') AS RESULT GROUP BY ALF_PE HAVING COUNT(*) = 1")

fam1 <- sqlQuery(channel, "SELECT DISTINCT (ALF_PE ) FROM 
    
                   SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 WHERE EVENT_CD 
                   
                   LIKE '12D2%' GROUP BY ALF_PE")


fam_hist_asthma <- fam %>% full_join(fam1, by = 'ALF_PE')

fam_hist_asthma$fam_hist_asthma <- 1

#------------------------------------------------------------------------------
# G. Join  hypothesis risk factors to the Asthma_COVid table and create table
#------------------------------------------------------------------------------
# Join risk factors to the study population table and rename

a_c_c1 <- a_c_c_pop %>% left_join(obesity, by = 'ALF_PE') %>% 
         
         left_join(sleep_apnea, by = 'ALF_PE') %>% 
         
         left_join(smoking, by = 'ALF_PE') %>% 
         
         left_join(alcohol, by = 'ALF_PE') %>%
         
         left_join(fam_hist_asthma, by = 'ALF_PE') %>%
         
         mutate(obesity = if_else(!is.na(obesity), 1, 0)) %>% 
         
         mutate(sleep_apnea = if_else(!is.na(sleep_apnea), 1, 0)) %>%
         
         mutate(smoking = if_else(!is.na(smoking), 1, 0)) %>%
         
         mutate(alcohol = if_else(!is.na(alcohol), 1, 0)) %>%
         
         mutate(fam_hist_asthma = if_else(!is.na(fam_hist_asthma), 1, 0))
  

## Join hypothesis risk factors to control cases table and create a  new table

f_control_cases1 <- f_control_cases_pop %>% left_join(obesity, by = 'ALF_PE') %>% 
          
          left_join(sleep_apnea, by = 'ALF_PE') %>% 
          
          left_join(smoking, by = 'ALF_PE') %>% 
          
          left_join( alcohol, by = 'ALF_PE') %>%
          
          left_join(fam_hist_asthma, by = 'ALF_PE') %>%
          
          mutate(obesity = if_else(!is.na(obesity), 1, 0)) %>% 
          
          mutate(sleep_apnea = if_else(!is.na(sleep_apnea), 1, 0)) %>%
          
          mutate(smoking = if_else(!is.na(smoking), 1, 0)) %>%
          
          mutate(alcohol = if_else(!is.na(alcohol), 1, 0)) %>%
          
          mutate(fam_hist_asthma = if_else(!is.na(fam_hist_asthma), 1, 0))

##-----------------------------------------------------------------------------
## H. RISK FACTORS FOR MACHINE LERANING
#     From PEDW_SINGLE_DIAG_20220228, EDDS_EDDS_20220301 and
#     WLGP_GP_EVENT_CLEANSED_20220201
#-----------------------------------------------------------------------------

# Select ALFE_PE with risk factors for machine learning 

dep <- sqlQuery(channel, "SELECT ALF_PE FROM(SELECT ALF_PE FROM 

                    SAIL1281V.PEDW_SINGLE_DIAG_20220228 WHERE DIAG_CD_1234 
                    
                    LIKE 'F32%' UNION  SELECT ALF_PE FROM 
                    
                    SAIL1281V.EDDS_EDDS_20220301 WHERE DIAG_CD_1  LIKE 'F32%'
                    
                    UNION SELECT ALF_PE FROM 
                    
                    SAIL1281V.PEDW_SINGLE_DIAG_20220228 WHERE DIAG_CD_1234 IN
                    
                    ('Y492', 'F32A', 'F331', 'F330') UNION  SELECT ALF_PE FROM
                    
                    SAIL1281V.EDDS_EDDS_20220301 WHERE DIAG_CD_1 IN
                    
                    ('Y492', 'F32A', 'F331', 'F330')) 
                    
                    AS RESULT GROUP BY ALF_PE HAVING COUNT(*) = 1")

dep1 <- sqlQuery(channel, "SELECT DISTINCT (ALF_PE ) FROM 
    
                   SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 WHERE EVENT_CD 
                   
                   LIKE 'E112.%' OR EVENT_CD LIKE 'E2B.%' GROUP BY ALF_PE")


depression <- dep %>% full_join(dep1, by = 'ALF_PE')

depression$depression <- 1


anx <- sqlQuery(channel, "SELECT ALF_PE FROM(SELECT ALF_PE FROM 

                    SAIL1281V.PEDW_SINGLE_DIAG_20220228 WHERE DIAG_CD_1234 
                    
                    LIKE 'F41%' UNION  SELECT ALF_PE FROM 
                    
                    SAIL1281V.EDDS_EDDS_20220301 WHERE DIAG_CD_1  LIKE 'F41%'
                    
                    UNION SELECT ALF_PE FROM 
                    
                    SAIL1281V.PEDW_SINGLE_DIAG_20220228 WHERE DIAG_CD_1234 IN
                    
                    ('E20', 'E200') UNION  SELECT ALF_PE FROM 
                    
                    SAIL1281V.EDDS_EDDS_20220301 WHERE DIAG_CD_1 IN
                    
                    ('E20', 'E200')) 
                    
                    AS RESULT GROUP BY ALF_PE HAVING COUNT(*) = 1")

anx1 <- sqlQuery(channel, "SELECT DISTINCT (ALF_PE ) FROM 

                   SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 WHERE EVENT_CD 
                   
                   LIKE 'E20%'  GROUP BY ALF_PE")

anxiety <- anx %>% full_join(anx1, by = 'ALF_PE')

anxiety$anxiety <- 1


alle <- sqlQuery(channel, "SELECT ALF_PE FROM(SELECT ALF_PE FROM

                    SAIL1281V.PEDW_SINGLE_DIAG_20220228 WHERE DIAG_CD_1234 
                    
                    IN ('Z888', 'Z889', 'T784') UNION  SELECT ALF_PE FROM 
                    
                    SAIL1281V.EDDS_EDDS_20220301 WHERE DIAG_CD_1 IN 
                    
                    ('Z888', 'Z889', 'T784') UNION SELECT ALF_PE FROM 
                    
                    SAIL1281V.PEDW_SINGLE_DIAG_20220228 WHERE DIAG_CD_1234 LIKE
                    
                    'J30%' UNION  SELECT ALF_PE FROM 
                    
                    SAIL1281V.EDDS_EDDS_20220301 WHERE DIAG_CD_1 LIKE 'J30%') 
                    
                    AS RESULT GROUP BY ALF_PE HAVING COUNT(*) = 1")


alle1 <- sqlQuery(channel, "SELECT DISTINCT (ALF_PE ) FROM

                     SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 WHERE EVENT_CD 
                     
                     LIKE 'H17%' GROUP BY ALF_PE")


allergies <- alle %>% full_join(alle1, by = 'ALF_PE')

allergies$allergies <- 1


GER <- sqlQuery(channel, "SELECT ALF_PE FROM(SELECT ALF_PE FROM

                    SAIL1281V.PEDW_SINGLE_DIAG_20220228 WHERE DIAG_CD_1234 
                    
                    IN ('K219', 'K210') UNION  SELECT ALF_PE FROM
                    
                    SAIL1281V.EDDS_EDDS_20220301 WHERE DIAG_CD_1 IN 
                    
                    ('K219', 'K210') UNION SELECT ALF_PE FROM 
                    
                    SAIL1281V.PEDW_SINGLE_DIAG_20220228 WHERE DIAG_CD_1234 =
                    
                    'T855' UNION  SELECT ALF_PE FROM 
                    
                    SAIL1281V.EDDS_EDDS_20220301 WHERE DIAG_CD_1 = 'T855') 
                    
                    AS RESULT GROUP BY ALF_PE HAVING COUNT(*) = 1")

GER1 <- sqlQuery(channel, "SELECT DISTINCT (ALF_PE ) FROM 

                      SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 WHERE EVENT_CD 
                      
                      LIKE 'K21%' GROUP BY ALF_PE")


GERD <- GER %>% full_join(GER1, by = 'ALF_PE')

GERD$GERD <- 1


sinu <- sqlQuery(channel, "SELECT ALF_PE FROM(SELECT ALF_PE FROM 

                    SAIL1281V.PEDW_SINGLE_DIAG_20220228 WHERE DIAG_CD_1234 
                    
                    LIKE 'J32%' UNION  SELECT ALF_PE FROM 
                    
                    SAIL1281V.EDDS_EDDS_20220301 WHERE DIAG_CD_1 LIKE 'J32%' 
                    
                    UNION SELECT ALF_PE FROM 
                    
                    SAIL1281V.PEDW_SINGLE_DIAG_20220228 WHERE DIAG_CD_1234 =
                    
                    'H135' UNION  SELECT ALF_PE FROM 
                    
                    SAIL1281V.EDDS_EDDS_20220301 WHERE DIAG_CD_1 = 'H135')
                    
                    AS RESULT GROUP BY ALF_PE HAVING COUNT(*) = 1")


sinu1 <- sqlQuery(channel, "SELECT DISTINCT (ALF_PE ) FROM 
        
                    SAIL1281V.WLGP_GP_EVENT_CLEANSED_20220201 WHERE EVENT_CD 
                    
                    LIKE 'H01.%' GROUP BY ALF_PE")

sinusitis <- sinu %>% full_join(sinu1, by = 'ALF_PE')

sinusitis$sinusitis <- 1

#------------------------------------------------------------------------------

# I. Table for outcome of interest(Dead from Asthma or respiratory infections)
#   from ADDE_DEATHS_20220301
#------------------------------------------------------------------------------
outcome <- sqlQuery(channel, "SELECT DISTINCT (ALF_PE ), 
                      
                      DEATH_DT FROM SAIL1281V.ADDE_DEATHS_20220301 WHERE 
                      
                      DEATHCAUSE_DIAG_1_CD LIKE 'J45%' OR DEATHCAUSE_DIAG_1_CD 
                      
                      LIKE 'J44%' OR DEATHCAUSE_DIAG_1_CD LIKE 'J46%' OR 
                      
                      DEATHCAUSE_DIAG_1_CD LIKE 'H5B0%' OR DEATHCAUSE_DIAG_1_CD
                      
                      LIKE 'J99%' OR DEATHCAUSE_DIAG_1_CD LIKE '663%' OR
 
                      DEATHCAUSE_DIAG_1_CD LIKE 'J98%' OR DEATHCAUSE_DIAG_1_CD 
                      
                      LIKE 'J96%' OR DEATHCAUSE_DIAG_1_CD LIKE 'J22%' OR
                      
                      DEATHCAUSE_DIAG_1_CD LIKE 'J80%' OR DEATHCAUSE_DIAG_1_CD
                      
                      LIKE 'I46%' OR DEATHCAUSE_DIAG_2_CD LIKE 'J45%' OR 
  
                      DEATHCAUSE_DIAG_2_CD LIKE 'J46%' OR DEATHCAUSE_DIAG_2_CD 
                      
                      LIKE 'J99%' OR DEATHCAUSE_DIAG_2_CD LIKE 'J98%' OR
                      
                      DEATHCAUSE_DIAG_2_CD LIKE 'J96%' OR DEATHCAUSE_DIAG_2_CD 
                      
                      LIKE 'J22%' OR DEATHCAUSE_DIAG_2_CD LIKE 'J80%' OR
                      
                      DEATHCAUSE_DIAG_2_CD LIKE 'I46%' OR
                      
                      DEATHCAUSE_DIAG_3_CD LIKE 'J45%' OR
                      
                      DEATHCAUSE_DIAG_3_CD LIKE 'J46%' OR
                      
                      DEATHCAUSE_DIAG_3_CD LIKE 'J99%' OR
                      
                      DEATHCAUSE_DIAG_3_CD LIKE 'J98%' OR
                      
                      DEATHCAUSE_DIAG_3_CD LIKE 'J96%' OR
                      
                      DEATHCAUSE_DIAG_3_CD LIKE 'J22%' OR
                      
                      DEATHCAUSE_DIAG_3_CD LIKE 'J80%' OR
                      
                      DEATHCAUSE_DIAG_3_CD LIKE 'I46%' OR
                      
                      DEATHCAUSE_DIAG_4_CD LIKE 'J45%' OR
                      
                      DEATHCAUSE_DIAG_4_CD LIKE 'J46%' OR 
                      
                      DEATHCAUSE_DIAG_4_CD LIKE 'J99%' OR
                      
                      DEATHCAUSE_DIAG_4_CD LIKE 'J98%' OR
                      
                      DEATHCAUSE_DIAG_4_CD LIKE 'J96%' OR
                      
                      DEATHCAUSE_DIAG_4_CD LIKE 'J22%' OR
                      
                      DEATHCAUSE_DIAG_4_CD LIKE 'J80%' OR
                      
                      DEATHCAUSE_DIAG_4_CD LIKE 'I46%' OR
                      
                      DEATHCAUSE_DIAG_5_CD LIKE 'J45%' OR
                      
                      DEATHCAUSE_DIAG_5_CD LIKE 'J46%' OR 
                      
                      DEATHCAUSE_DIAG_5_CD LIKE 'J99%' OR
                      
                      DEATHCAUSE_DIAG_5_CD LIKE 'J98%' OR
                      
                      DEATHCAUSE_DIAG_5_CD LIKE 'J96%' OR
                      
                      DEATHCAUSE_DIAG_5_CD LIKE 'J22%' OR
                      
                      DEATHCAUSE_DIAG_5_CD LIKE 'J80%' OR
                      
                      DEATHCAUSE_DIAG_5_CD LIKE 'I46%' OR
                      
                      DEATHCAUSE_DIAG_6_CD LIKE 'J45%' OR
                      
                      DEATHCAUSE_DIAG_6_CD LIKE 'J46%' OR
                      
                      DEATHCAUSE_DIAG_6_CD LIKE 'J99%' OR
                      
                      DEATHCAUSE_DIAG_6_CD LIKE 'J98%' OR
                      
                      DEATHCAUSE_DIAG_6_CD LIKE 'J96%' OR
                      
                      DEATHCAUSE_DIAG_6_CD LIKE 'J22%' OR
                      
                      DEATHCAUSE_DIAG_6_CD LIKE 'J80%' OR
                      
                      DEATHCAUSE_DIAG_6_CD LIKE 'I46%'")


colnames(outcome)[2]<- 'SAA_EVENT'

outcome$outcome_var <- 1

#------------------------------------------------------------------------------
# J. Table of all recorded death from ADDE_DEATHS_20220301
#------------------------------------------------------------------------------

dod <- sqlQuery(channel, "SELECT DISTINCT (ALF_PE), DEATH_DT FROM 
                   
                SAIL1281V.ADDE_DEATHS_20220301 ad")

dod <- na.omit(dod)

colnames(dod)[2] <- 'DEATH_DATE'

dod$dead <- 1

#-----------------------------------------------------------------------------
# K. Attach outcome of interest and dead to Study and Control population and
#    restrict population to age >= 5 years and those who died before the start
#    of the study period
#------------------------------------------------------------------------------
# Attach outcome of interest and Death of death to our Tables
#  and selecting period of study for Hypothesis

asthma_covid_hypo <- a_c_c1 %>% left_join(outcome, by = 'ALF_PE') %>% 
            
                    mutate(outcome_var = if_else(!is.na(outcome_var), 
                                                       
                    1, 0)) %>% left_join(dod, by = 'ALF_PE') %>% mutate(dead =
                                                                          
                    if_else(!is.na(dead), 1, 0))


## Restrict Study to exclude population with death date > study period

asthma_covid_hypo <- asthma_covid_hypo %>% filter(is.na(DEATH_DATE) | !is.na(
  
                      DEATH_DATE) & DEATH_DATE>"2018-12-31") %>% 
  
                      filter(is.na(SAA_EVENT) | !is.na(SAA_EVENT) & 
                               
                      SAA_EVENT>"2018-12-31") %>% filter(age_at_diag >=5) 


#-----------------------------------------------------------------------------
## Control Population for the hypothesis

control_pop_hypo <- f_control_cases1 %>% left_join(outcome, by = 'ALF_PE') %>% 
           
                     mutate(outcome_var = if_else(!is.na(outcome_var),
                                                  
                     1, 0)) %>% left_join(dod, by ='ALF_PE') %>% mutate(dead =
                                                                          
                     if_else(!is.na(dead), 1, 0))


## Restrict Study to exclude population with death date > study period

control_pop_hypo <- control_pop_hypo  %>% filter(is.na(DEATH_DATE) | !is.na(
  
                    DEATH_DATE) & DEATH_DATE>"2018-12-31") %>% 
  
                    filter(is.na(SAA_EVENT) | !is.na(SAA_EVENT) & 
                             
                    SAA_EVENT>"2018-12-31") %>% filter(age_at_diag >=5)


# Randomly select control population equals to study population

# control_pop_hypo1 <- control_pop_hypo[sample(nrow(control_pop_hypo), "40321"), ]

#-----------------------------------------------------------------------------
## Send Study Population and Control Population to CSV files

# write.csv(control_pop_hypo1, file =  "control_pop_hypo1.csv")

# write.csv(asthma_covid_hypo, file = "asthma_covid_hypo.csv")


#----------------------------------------------------------------------------
## L. Attach machine learning risk factors to the created study and control
#   Population.
#-----------------------------------------------------------------------------
# Attach the machine learning risk factors to the Asthma_COVid table

asthma_covid_ML <- asthma_covid_hypo %>% left_join(depression,
                                                        
                                                          by = 'ALF_PE') %>%
  
                 left_join(anxiety, by = 'ALF_PE') %>%
  
                 left_join(allergies, by = 'ALF_PE') %>%
  
                 left_join(GERD, by = 'ALF_PE') %>%
  
                 left_join(sinusitis, by = 'ALF_PE') %>%
  
                 mutate(depression = if_else(!is.na(depression), 1, 0)) %>%
  
                 mutate(anxiety = if_else(!is.na(anxiety), 1, 0)) %>%
  
                 mutate(allergies = if_else(!is.na(allergies), 1, 0)) %>%
  
                 mutate(GERD = if_else(!is.na(GERD), 1, 0)) %>%
  
                 mutate(sinusitis = if_else(!is.na(sinusitis), 1, 0))


## Attach the machine learning risk factors to the control population table

control_pop_ML <- control_pop_hypo1 %>% left_join(depression, by = 'ALF_PE') %>%
  
                     left_join(anxiety, by = 'ALF_PE') %>% 
  
                     left_join(allergies, by = 'ALF_PE') %>%
  
                     left_join(GERD, by = 'ALF_PE') %>%
  
                     left_join(sinusitis, by = 'ALF_PE') %>%
  
                     mutate(depression = if_else(!is.na(depression), 1, 0)) %>%
  
                     mutate(anxiety = if_else(!is.na(anxiety), 1, 0)) %>%
  
                     mutate(allergies = if_else(!is.na(allergies), 1, 0)) %>%
  
                     mutate(GERD = if_else(!is.na(GERD), 1, 0)) %>%
  
                     mutate(sinusitis = if_else(!is.na(sinusitis), 1, 0))

#-----------------------------------------------------------------------------
## M. Analysis for Logistics Regression
#      CONDUCT LOGISTICS REGRESSION FOR THE STUDY POPULATION
#------------------------------------------------------------------------------

library(rpart)                 # For Decision trees
library(caret)                 # Machine learning
install.packages('broom')
library(broom)                 # For odds ratio

## Remove variables not important to GLM
asthma_covid_glm <- asthma_covid_hypo %>% select(-c(ALF_PE, asthma_onset,
                                                    
                                covid_onset, SAA_EVENT, WOB, DEATH_DATE, dead))


par(mfrow=c(2,2), xpd=NA)

# show pairs relationship

pairs(asthma_covid_glm, panel = panel.smooth)

# Fitting the maximal model

model1 <- glm(outcome_var ~ gender + age_at_diag + obesity + 
                
               sleep_apnea + smoking + alcohol + fam_hist_asthma, 
              
               data = asthma_covid_glm, family = "binomial")

summary(model1)

# odds ratio

model1 %>% tidy(conf.int = TRUE, exp = TRUE)

# check and update model

model2 <- step(model1)

aba <- summary(model2)

abr <- aba$coefficients

plot(model2)


# Save significant risk factors table

write.csv(abr, file = "asthma_covid_glm_coeff.csv")


## CONDUCT LOGISTICS REGRESSION FOR THE CONTROL POPULATION

# Remove variables not important to GLM

control_pop_glm <- control_pop_hypo1 %>% select(-c(ALF_PE, asthma_onset,
                                   
                                            SAA_EVENT, WOB, DEATH_DATE, dead))

# show pairs relationship

pairs(control_pop_glm, panel = panel.smooth)

# Fitting the maximal model

modelc1 <- glm(outcome_var ~ gender + age_at_diag + obesity + 
                 
                sleep_apnea + smoking + alcohol + fam_hist_asthma, 
               
                data = control_pop_glm, family = 'binomial')

summary(modelc1)

# odds ratio

modelc1 %>% tidy(conf.int = TRUE, exp = TRUE)

# check and update model

modelc2 <- step(modelc1)

bb <-summary(modelc2)

abk <- bb$coefficients

par(mfrow=c(2,2), xpd=NA)

plot(modelc2)

# Save significant coeff for Control population

write.csv(abk, file = "control_pop_glm_coeff.csv")


#-----------------------------------------------------------------------------
# N.  Machine learning Analysis for Study Population and Control population
#     using decision tree
#------------------------------------------------------------------------------

# Machine learning for study Population 
# Remove variables not important to Machine learning
asthma_covid_dt <- asthma_covid_ML %>% select(-c(ALF_PE, asthma_onset,
                                                 
                                 covid_onset, SAA_EVENT, WOB, DEATH_DATE, dead))

set.seed(1792)

# Partition data for ASTHMA and COVID  Population
df <- createDataPartition(asthma_covid_dt$outcome_var, p=0.9, list = F)

train <- asthma_covid_dt[df, ]

test <- asthma_covid_dt[-df, ]

par(mfrow=c(1,1), xpd=NA)

fit <- rpart(outcome_var~., data =train, control =rpart.control(cp=0.0001))

plot(fit)

text(fit, cex=0.5, xpd=TRUE)

printcp(fit)

fit_pruned <- prune(fit, cp=0.00084)

fit_pruned$frame$yval <- round(fit_pruned$frame$yval, digits = 3)

plot(fit_pruned, cex=0.6, xpd=TRUE)

text(fit_pruned, cex=0.6, xpd=TRUE)

pred <- predict(fit_pruned, test)

pred <- data.frame(pred)

pred <- pred %>% mutate(pred = ifelse(pred >= 0.5, 1, 0))

pred <- unlist(pred)


#ROC Curve and AUC for asthma patients with COVID-19

modelr < rpart(data = train, family = binomial, outcome_var~gender + 
                 
                 age_at_diag + obesity + sleep_apnea + smoking + alcohol + 
                 
                 fam_hist_asthma + allergies + GERD + sinusitis + depression)

predroc <- predict(modelr, test, type = "response")

preditionroc <- prediction(predictions = predroc, labels = test$outcome_var)

perfroc <- performance(preditionroc, "tpr", "fpr")

plot(perfroc, lwd=3, colorize=TRUE)

auc <- performance(preditionroc, "auc")@y.values[[1]]


##  Machine Learning for the Control population 

# Remove variables not important to GLM
control_pop_dt <- control_pop_ML %>% select(-c(ALF_PE, asthma_onset,
                                               
                                   SAA_EVENT, WOB, DEATH_DATE, dead))
set.seed(2333)
     
# Partition data for ASTHMA and COVID  Population

df1 <- createDataPartition(control_pop_dt$outcome_var, p=0.9, list = F)

train1 <- control_pop_dt[df1, ]

test1 <- control_pop_dt[-df1, ]
     
par(mfrow=c(1,1), xpd=NA)

fit1 <- rpart(outcome_var~., data =train1, control =rpart.control(cp=0.0001))

plot(fit1)

text(fit1, cex=0.5, xpd=TRUE)

printcp(fit1)

fit_pruned1 <- prune(fit1, cp=0.001)

fit_pruned1$frame$yval <- round(fit_pruned1$frame$yval, digits = 3)

plot(fit_pruned1, cex=0.6, xpd=TRUE)

text(fit_pruned1, cex=0.6, xpd=TRUE)

pred1 <- predict(fit_pruned1, test1)

pred1 <- data.frame(pred1)

pred1 <- pred1 %>% mutate(pred1 = ifelse(pred >= 0.5, 1, 0))

pred1 <- unlist(pred1)

table(pred1)

# ROC Curve and AUC for asthma patients 

modelz < rpart(data = train1, family = binomial, outcome_var~gender + 
                 
                 age_at_diag + obesity + sleep_apnea + smoking + alcohol + 
                 
                 fam_hist_asthma + allergies + GERD + sinusitis + depression)

predroc1 <- predict(modelz, test1, type = "response")

preditionroc1 <- prediction(predictions = predroc1, labels = test1$outcome_var)

perfroc1 <- performance(preditionroc1, "tpr", "fpr")

plot(perfroc1, lwd=3, colorize=TRUE)

auc1 <- performance(preditionroc1, "auc")@y.values[[1]]


# Plot Combined ROCS and AUC

plot(perfroc, col = "blue", lwd=3, print.auc=TRUE, main="ROC Curves for Study
     
     and Control Population")

text(0.5, 0.4, paste("AUC", round(auc, 2)), col="blue", cex=1)

plot(perfroc1, col="red", lwd=3, add=TRUE)

text(0.6, 0.5, paste("AUC", round(auc1, 2)), col="red", cex=1)

abline(a=0, b=1, lwd=2, lty=3, col="grey")

legend('bottomright', c('Asthma Patients', 'Asthma Patients with COVID'),
       
       lwd = 2, col = c("blue", "red"))



#------------------------------------------------------------------------------
# O. Time_to_ event and COX proportional hazard Regression analysis for both 
# control and study population.
#-----------------------------------------------------------------------------
install.packages('survminer')
install.packages('psych')


library(survival)       # The survival object.
library(survminer)      # survival plots. 
library(psych)         # describe()



#------------------------------------------------------------------------------
# Load, tidy and explore Data
#------------------------------------------------------------------------------
#Assign cohort to the two study population and change date format
asthma_covid_ML$cohort <- 1

asthma_covid_ML$DEATH_DATE <- as.Date(asthma_covid_ML$DEATH_DATE, 
                                                        origin = "1970-01-01")

control_pop_ML$covid_onset <- NA

control_pop_ML$cohort <- 2

str(asthma_covid_ML)

asthma_covid_ML$SAA_EVENT <- as.Date(asthma_covid_ML$SAA_EVENT)

asthma_covid_ML <- asthma_covid_ML %>% mutate(begin = 
  ifelse(outcome_var == 1 & covid_onset > "2019-01-01",
         
  covid_onset, as.Date("2019-01-01", origin = "1970-01-01"))) %>%
  
  mutate(end = ifelse(dead == 1 & outcome_var == 0 & DEATH_DATE < "2022-12-31",
                      
  DEATH_DATE, as.Date("2022-12-31", origin = "1970-01-01"))) %>%
  
  mutate(end = ifelse(outcome_var==1 & SAA_EVENT < "2022-12-31",
                      
  SAA_EVENT, as.Date("2022-12-31", origin = "1970-01-01"))) %>%
  
  mutate(begin = as.Date(begin, origin = "1970-01-01")) %>% 
  
  mutate(end = as.Date(end, origin = "1970-01-01")) %>%
  
  mutate(time_to_event = round(end - begin, 0))


table(asthma_covid_ML$time_to_event)


str(control_pop_ML)

control_pop_ML$SAA_EVENT <- as.Date(control_pop_ML$SAA_EVENT)

control_pop_ML$DEATH_DATE <- as.Date(control_pop_ML$DEATH_DATE)

control_pop_ML <- control_pop_ML %>% mutate(begin = 
                 
  ifelse(outcome_var == 1 & asthma_onset > "2019-01-01",
                 
  asthma_onset, as.Date("2019-01-01", origin = "1970-01-01"))) %>%
  
  mutate(end = ifelse(dead == 1 & outcome_var == 0 & DEATH_DATE < "2022-12-31",
                      
  DEATH_DATE, as.Date("2022-12-31", origin = "1970-01-01"))) %>%
  
  mutate(end = ifelse(outcome_var==1 & SAA_EVENT < "2022-12-31",
                      
  SAA_EVENT, as.Date("2022-12-31", origin = "1970-01-01"))) %>%
  
  mutate(begin = as.Date(begin, origin = "1970-01-01")) %>% 
  
  mutate(end = as.Date(end, origin = "1970-01-01")) %>%
  
  mutate(time_to_event = round(end - begin, 0))
  


# rbind cohort and arrange by ALF_PE for analysis   
asthma_covid_surv <- rbind(asthma_covid_ML, control_pop_ML) %>% 
  
                                                           arrange(ALF_PE)

 
## Time_to_event Analysis
surv_model <- survfit(Surv(ceiling(time_to_event/365), outcome_var)~cohort, 
                                                    
                                                    data = asthma_covid_surv)


plot(surv_model, main = "Time to death for study cohort in years ",
     
     xlab = "Time(Years)",
     
     ylab = "Survival Probobility", col =c("blue", "red"), lwd = 2,
     
     ylim = c(.98, 1))

legend('bottomleft', c('Asthma Patients', 'Covid and Asthma Patients'), 
    
        lwd = 2, col =c("blue", "red"))


## Smooth the curve plot by creation of strata

# Create the data frames

s1 <- surv_model$strata[1]

s2 <- surv_model$strata[2]

# Combine the data frames

pls1 <- data.frame(time = surv_model$time[1:s1],
                   
                        chz = 1 - surv_model$cumhaz[1:s1])

pls2 <- data.frame(time = surv_model$time[(s1+1):(s1+s2)],
                   
                         chz1 = 1 - surv_model$cumhaz[(s1+1):(s1+s2)])

pl <- pls1 %>% full_join(pls2, by = 'time') %>% filter(time >= 0)

# Create the probability plot using ggplot

ggplot(data = pl, aes(x = time), conf.int = FALSE) +
       
       geom_smooth(aes(y = chz, color = 'COVID and Asthma Patients'),
                   
                   method = 'loess') +
       
       geom_smooth(aes(y = chz1, color = 'Asthma Patients'), method = 
                     
                     'loess') +
       
       ggtitle(' COVID and Asthma Patients Versus Asthma Patients Survival
               
                Probability in months') + labs(y = 'survival probability') +
       
        theme(legend.position = "bottom")



#Test for survival  difference between the cohorts

survdiff(Surv(time_to_event, outcome_var) ~ cohort, data = asthma_covid_surv)



# Cox proportional hazard regression analysis

survtime <- Surv(asthma_covid_surv$time_to_event, asthma_covid_surv$outcome_var)

cox_model <- coxph(survtime~ gender + age_at_diag + obesity + anxiety + 
                     
                     sleep_apnea + smoking + alcohol + fam_hist_asthma +
                     
                     allergies + GERD + sinusitis + depression, 
                   
                     data = asthma_covid_surv)



mina <- summary(cox_model)

mina1 <- step(cox_model)

mina2 <- summary(mina1)

bnn <- mina2$coefficients

bnx <- mina2$conf.int


# Plot the hazard ratio for factors with hazard for outcome variable

cox_h <- coxph(survtime~ age_at_diag + obesity + sleep_apnea + smoking + 
                 
              fam_hist_asthma + GERD + sinusitis + depression, 
              
              data = asthma_covid_surv)

ggforest(cox_h, data = asthma_covid_surv)



## Save coefficients of cox hazard analysis

write.csv(bnn, file = "cox_hazard_coeff.csv")

write.csv(bnx, file = "cox_hazard_ratio.csv")

## Schoenfeld Test

ggcoxzph(cox.zph(cox_model))

#------------------------------------------------------------------------------
#P.  Descriptive Statistics
-------------------------------------------------------------------------------
#  Age distribution for studies population


ggplot(data = asthma_covid_surv, aes(x = age_at_diag)) + 
  
      geom_histogram(binwidth=15,
          
      fill = "grey", alpha=0.9) + theme_classic() +
      
      ggtitle("Distribution of studies Population by Age") + 
  
      xlab("Age") + ylab("Frequency")


## Table of events for Study population

install.packages('gtsummary')

library(gtsummary)                # for table summary

# Remove variables not important 

com1 <- asthma_covid_ML %>% select(-c(ALF_PE, asthma_onset, covid_onset, WOB,
                                      
                                      SAA_EVENT, DEATH_DATE, dead, cohort,
                                      
                                      begin, end, time_to_event))



# Convert variables to character

asthm33 <- asthma_covid_surv %>% mutate(outcome_var = ifelse(
  
                  outcome_var == 1, "Died", "Alive"),
                                        
                  gender = ifelse(gender == 1, "Female", "Male"),
                                        
                  obesity = ifelse(obesity ==1, "Yes", "No"),
                            
                  sleep_apnea = ifelse(sleep_apnea ==1, "Yes", "No"),
                                        
                  smoking = ifelse(smoking == 1, "Yes", "No"),
                                        
                  alcohol = ifelse(alcohol == 1, "Yes", "No"),
                                        
                  fam_hist_asthma = ifelse(fam_hist_asthma == 1, "Yes", "No"),
                                        
                  depression = ifelse(depression == 1, "Yes", "No"),
                          
                  anxiety = ifelse(anxiety == 1, "Yes", "No"),
                          
                  GERD = ifelse(GERD == 1, "Yes", "No"),
                                        
                  allergies = ifelse(allergies == 1, "Yes", "No"),
                              
                  sinusitis = ifelse(sinusitis == 1, "Yes", "No"),
                  
                  cohort = ifelse(cohort == 1, "Asthma Patients with COVID-19
                                  
                  Infection", "Asthma Patients"))


## Generate Characteristics of the Research Population by Cohort

asth <- asthm33 %>% select(-c(ALF_PE, asthma_onset, covid_onset, WOB, 
                              
                              SAA_EVENT, DEATH_DATE, dead, begin, end,
                              
                              time_to_event))

# Characteristics of cohorts by outcome variable

asth %>% tbl_strata(strata = cohort, ~.x) %>% tbl_summary(by = outcome_var)
    
# Characteristics of the whole research population 

asth %>% tbl_summary(by = cohort,
                     
                     statistic = list(
                       
                     age_at_diag ~ "{mean} ({sd})")) %>% add_p() %>% 
  
                     add_n() %>% add_ci() %>% add_overall() %>% 
  
                     add_stat_label(label = all_continous() ~ "Median (IQR)")

#-------------------------------------------------------------------------------
#-===================THE END===================================================

                     
