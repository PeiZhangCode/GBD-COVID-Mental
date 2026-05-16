######################################################################################################
# this file is for constructing policy indices
# Pei Zhang
# 2026-04-30
######################################################################################################
### packages
packages <- c('foreign','tsModel','lmtest','Epi','splines','vcd','scales','dplyr','ggplot2','ggsci',
              'RColorBrewer','ggpubr','tidyr','Epi','ggtext','lubridate','countrycode','sf','RColorBrewer',
              'patchwork','gtable','grid','broom','performance','forcats','readr')
lapply(packages, require, character.only=T) 

dat0 <- read.csv("OxCGRT_compact_national_v1.csv")

names(dat0) <- c(
  "CountryName",
  "CountryCode",
  "RegionName",
  "RegionCode",
  "Jurisdiction",
  "Date",
  "C1",
  "C1M_Flag",
  "C2",
  "C2M_Flag",
  "C3",
  "C3M_Flag",
  "C4",
  "C4M_Flag",
  "C5",
  "C5M_Flag",
  "C6",
  "C6M_Flag",
  "C7",
  "C7M_Flag",
  "C8",
  "E1",
  "E1_Flag",
  "E2",
  "E3",
  "E4",
  "H1",
  "H1_Flag",
  "H2",
  "H3",
  "H4",
  "H5",
  "H6",
  "H6M_Flag",
  "H7",
  "H7_Flag",
  "H8",
  "H8M_Flag",
  "V1",
  "V2",
  "V2B_Vaccine.age.eligibility.availability.age.floor..general.population.summary.",
  "V2C_Vaccine.age.eligibility.availability.age.floor..at.risk.summary.",
  "V2D_Medically..clinically.vulnerable..Non.elderly.",
  "V2E_Education",
  "V2F_Frontline.workers...non.healthcare.",
  "V2G_Frontline.workers...healthcare.",
  "V3",
  "V4",
  "ConfirmedCases",
  "ConfirmedDeaths",
  "MajorityVaccinated",
  "PopulationVaccinated",
  "StringencyIndex_Average",
  "GovernmentResponseIndex_Average",
  "ContainmentHealthIndex_Average",
  "EconomicSupportIndex")

vars_to_binary <- c("E3", "E4", "H4", "H5")
dat0[vars_to_binary] <- lapply(dat0[vars_to_binary], function(x) ifelse(x > 0, 1, 0))

dat0 <- dat0 %>%
  mutate(
    C1 = case_when(
      C1 == 0 ~ 0,
      !is.na(C1M_Flag) ~ 100 * (C1 + C1M_Flag) / (3 + 1),
      TRUE ~ 100 * C1 / 3),
    
    C2 = case_when(
      C2 == 0 ~ 0,
      !is.na(C2M_Flag) ~ 100 * (C2 + C2M_Flag) / (3 + 1),
      TRUE ~ 100 * C2 / 3),
    
    C3 = case_when(
      C3 == 0 ~ 0,
      !is.na(C3M_Flag) ~ 100 * (C3 + C3M_Flag) / (2 + 1),
      TRUE ~ 100 * C3 / 2),
    
    C4 = case_when(
      C4 == 0 ~ 0,
      !is.na(C4M_Flag) ~ 100 * (C4 + C4M_Flag) / (4 + 1),
      TRUE ~ 100 * C4 / 4),
    
    C5 = case_when(
      C5 == 0 ~ 0,
      !is.na(C5M_Flag) ~ 100 * (C5 + C5M_Flag) / (2 + 1),
      TRUE ~ 100 * C5 / 2),
    
    C6 = case_when(
      C6 == 0 ~ 0,
      !is.na(C6M_Flag) ~ 100 * (C6 + C6M_Flag) / (3 + 1),
      TRUE ~ 100 * C6 / 3),
    
    C7 = case_when(
      C7 == 0 ~ 0,
      !is.na(C7M_Flag) ~ 100 * (C7 + C7M_Flag) / (2 + 1),
      TRUE ~ 100 * C7 / 2),
    
    C8 = case_when(
      C8 == 0 ~ 0,
      TRUE ~ 100 * C8 / 4))


dat0 <- dat0 %>%
  mutate(
    E1 = case_when(
      E1 == 0 ~ 0,
      !is.na(E1_Flag) ~ 100 * (E1 + E1_Flag) / (2 + 1),
      TRUE ~ 100 * E1 / 2),
    
    E2 = case_when(
      E2 == 0 ~ 0,
      TRUE ~ 100 * E2 / 2),
    
    E3 = case_when(
      E3 == 0 ~ 0,
      TRUE ~ 100 * as.numeric(E3) / 1),    
    
    E4 = case_when(
      E4 == 0 ~ 0,
      TRUE ~ 100 * as.numeric(E4) / 1))

dat0 <- dat0 %>%
  mutate(
    H1 = case_when(
      H1 == 0 ~ 0,
      !is.na(H1_Flag) ~ 100 * (H1 + H1_Flag) / (2 + 1),
      TRUE ~ 100 * H1 / 2),
    
    H2 = case_when(
      H2 == 0 ~ 0,
      TRUE ~ 100 * H2 / 3),
    
    H3 = case_when(
      H3 == 0 ~ 0,
      TRUE ~ 100 * H3 / 2),
    
    H4 = case_when(
      H4 == 0 ~ 0,
      TRUE ~ 100 * as.numeric(H4) / 1),    
    
    H5 = case_when(
      H5 == 0 ~ 0,
      TRUE ~ 100 * as.numeric(H5) / 1),   
    
    H6 = case_when(
      H6 == 0 ~ 0,
      !is.na(H6M_Flag) ~ 100 * (H6 + H6M_Flag) / (4 + 1),
      TRUE ~ 100 * H6 / 4),
    
    H7 = case_when(
      H7 == 0 ~ 0,
      !is.na(H7_Flag) ~ 100 * (H7 + H7_Flag) / (5 + 1),
      TRUE ~ 100 * H7 / 5),
    
    H8 = case_when(
      H8 == 0 ~ 0,
      !is.na(H8M_Flag) ~ 100 * (H8 + H8M_Flag) / (3 + 1),
      TRUE ~ 100 * H8 / 3))

dat0 <- dat0 %>%
  mutate(
    V1 = case_when(
      V1 == 0 ~ 0,
      TRUE ~ 100 * as.numeric(V1) / 2),
    
    V2 = case_when(
      V2 == 0 ~ 0,
      TRUE ~ 100 * as.numeric(V2) / 3),
    
    V3 = case_when(
      V3 == 0 ~ 0,
      TRUE ~ 100 * as.numeric(V3) / 5),
    
    V4 = case_when(
      V4 == 0 ~ 0,
      TRUE ~ 100 * as.numeric(V4) / 1))  

###### policy indices construction ######

dat1 <- dat0 %>% mutate(ContainmentIndex = rowMeans(select(., C1,C2,C3,C4,C5,C6,C7,C8), na.rm = TRUE),
                        EconomicSupportIndex = rowMeans(select(., E1,E2), na.rm = TRUE),
                        HealthSystemIndex = rowMeans(select(., H1,H2,H3,H6,H7,H8), na.rm = TRUE),
                        VaccinationIndex = rowMeans(select(., V1,V2,V3,V4) %>% mutate_all(as.numeric), na.rm = TRUE)) %>%
  select(location_name=CountryName,RegionName,Date,ContainmentIndex,EconomicSupportIndex,
         HealthSystemIndex,VaccinationIndex)

write.csv(dat1, file = "policy_indices_four.csv")

dat1 <- dat0 %>% mutate(ContainmentIndex = rowMeans(select(., C1,C2,C3,C4,C5,C6,C7), na.rm = TRUE),
                        EconomicSupportIndex = rowMeans(select(., E1,E2), na.rm = TRUE),
                        HealthSystemIndex = rowMeans(select(., H1,H2,H3,H7,H8), na.rm = TRUE),
                        VaccinationIndex = rowMeans(select(., V1,V2,V3,V4) %>% mutate_all(as.numeric), na.rm = TRUE),
                        BorderIndex = C8,
                        FacemaskIndex = H6) %>%
  select(location_name=CountryName,RegionName,Date,ContainmentIndex,EconomicSupportIndex,
         HealthSystemIndex,VaccinationIndex,BorderIndex,FacemaskIndex)

write.csv(dat1, file = "policy_index_six.csv")

###### alternative policy indices construction for sensitivity analysis ######
# dat1 <- dat0 %>% mutate(ContainmentIndex = rowMeans(select(., C1,C2,C3,C4,C5,C6,C7,C8), na.rm = TRUE),
#                         EconomicSupportIndex = rowMeans(select(., E1,E2,E3,E4), na.rm = TRUE),
#                         HealthSystemIndex = rowMeans(select(., H1,H2,H3,H4,H5,H6,H7,H8), na.rm = TRUE),
#                         VaccinationIndex = rowMeans(select(., V1,V2,V3,V4) %>% mutate_all(as.numeric), na.rm = TRUE)) %>%
#   select(location_name=CountryName,RegionName,Date,ContainmentIndex,EconomicSupportIndex,
#          HealthSystemIndex,VaccinationIndex)
# 
# write.csv(dat1, file = "policy_indices_four_all.csv")

# dat1 <- dat0 %>% mutate(ContainmentIndex = rowMeans(select(., C1,C2,C3,C4,C5,C6,C7), na.rm = TRUE),
#                         EconomicSupportIndex = rowMeans(select(., E1,E2,E3,E4), na.rm = TRUE),
#                         HealthSystemIndex = rowMeans(select(., H1,H2,H3,H4,H5,H7,H8), na.rm = TRUE),
#                         VaccinationIndex = rowMeans(select(., V1,V2,V3,V4) %>% mutate_all(as.numeric), na.rm = TRUE),
#                         BorderIndex = C8,
#                         FacemaskIndex = H6) %>%
#   select(location_name=CountryName,RegionName,Date,ContainmentIndex,EconomicSupportIndex,
#          HealthSystemIndex,VaccinationIndex,BorderIndex,FacemaskIndex)
# 
# write.csv(dat1, file = "policy_index_six_all.csv")
