######################################################################################################
# this file is for calculating excess mental disorder burden using ITS
# ARIMA + ARIMA-LSTM was used in main analysis, while ITS was for sensitivity analysis
# Pei Zhang
# 2026-04-30
######################################################################################################
### packages
packages <- c('foreign','tsModel','lmtest','Epi','splines','vcd','scales','dplyr','ggplot2','ggsci',
              'RColorBrewer','ggpubr','tidyr','Epi','ggtext','lubridate','countrycode','sf','RColorBrewer')
lapply(packages, require, character.only=T) 

###### 1. calculate excess for global ######
dat <- read.csv("IHME-GBD_2023_DATA-a33adbeb-1.csv") %>%
  select(location_id, location_name, measure_name, sex_name, age_name, cause_name, year, val) %>%
  filter(cause_name == "Mental disorders", year>2015)

names(dat) <- c("location_id","location","measure","sex","age","cause","year","val")

dat <- dat %>% arrange(year)
dat$time <- dat$year - 2015
dat$covid <- as.integer(dat$year >= 2020)
dat$time_after_covid <- pmax(0, dat$year - 2019)

model2 <- lm(val ~ time + covid + time_after_covid, data = dat)

### counterfactual
cf <- data.frame(
  year = 2020:2023,
  time = 2020:2023 - 2015,
  covid = 0,
  time_after_covid = 0)
pred <- predict(model2, newdata = cf, se.fit = TRUE)

mental_post <- dat %>%
  filter(year >= 2020) %>%
  mutate(
    predicted = pred$fit,                         
    se = pred$se.fit,                             
    ci_lower = predicted - 1.96 * se,
    ci_upper = predicted + 1.96 * se,
    extra_dalys = val - predicted,
    extra_low = val - ci_upper,                 
    extra_high = val - ci_lower) %>%
  select(year, predicted,
         extra_dalys, 
         extra_low, 
         extra_high)

### global excess mental health burden due to covid-19 during 2020-2023 (consider population change)
pop <- read.csv("IHME-GBD_2023_DATA-0154aaec-1.csv") %>% 
  select(year, popn=val) %>%
  mutate(popn=popn/100000)

mental_post2 <- mental_post %>%
  left_join(pop, by = "year")

summary_dalys <- mental_post2 %>%
  filter(year %in% 2020:2023) %>%
  summarise(
    extra = sum(extra_dalys * popn, na.rm = TRUE),
    extral = sum(extra_low   * popn, na.rm = TRUE),
    extrah = sum(extra_high  * popn, na.rm = TRUE),
    predicted = sum(predicted   * popn, na.rm = TRUE),
    prop_extra = round(extra / predicted * 100, 1),
    prop_extral = round(extral / predicted * 100, 1),
    prop_extrah = round(extrah / predicted * 100, 1))

summary_dalys

###### 2. calculate excess for countries #######
dat0 <- read.csv("IHME-GBD_2023_DATA-2b3e5b00-1.csv") %>% filter(year>2015)

country_ext <- function(country, case) {
  dat <- dat0 %>%
    filter(cause_name == case,
           metric_name == "Rate",
           location_id == country) %>%
    select(location_id, location_name, cause_name, year, val) %>%
    arrange(year)

  dat <- dat %>%
    mutate(
      covid = ifelse(year >= 2020, 1, 0),
      time = row_number())

  first_post_row <- which(dat$covid == 1)[1]
  dat <- dat %>%
    mutate(time_after_covid = ifelse(covid == 1, time - first_post_row + 1, 0))

  model2 <- lm(val ~ time + covid + time_after_covid, data = dat)

  post_years <- dat %>% filter(year >= 2020) %>% pull(year)
  if (length(post_years) == 0) return(NULL)  

  cf <- tibble(
    year = post_years,
    time = match(post_years, dat$year),
    covid = 0,
    time_after_covid = 0)

  pred <- predict(model2, newdata = cf, se.fit = TRUE)

  mental_post <- dat %>%
    filter(year %in% post_years) %>%
    mutate(
      predicted = pred$fit,
      se = pred$se.fit,
      ci_lower = predicted - 1.96 * se,
      ci_upper = predicted + 1.96 * se,
      ave_ext = val - predicted,
      ave_extl = val - ci_upper,
      ave_exth = val - ci_lower) %>%
    select(location_id, location_name, cause_name, year,
           ave_ext, ave_extl, ave_exth)
  return(mental_post)
}

country_list <- unique(dat0$location_id)
results_year <- purrr::map_dfr(country_list, ~country_ext(.x, "Mental disorders"))

names(results_year) <- c("location_id","location","cause","year","excess","excess_lower","excess_upper")

# write.csv(results_year, file = "ITS1623_excess_country.csv")

