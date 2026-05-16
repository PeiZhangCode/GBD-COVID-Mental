######################################################################################################
# this file is for interaction analysis between policy indices and excess mental disorder burden 
# Pei Zhang
# 2026-04-30
######################################################################################################
### packages
packages <- c('foreign','tsModel','lmtest','Epi','splines','vcd','scales','dplyr','ggplot2','ggsci',
              'RColorBrewer','ggpubr','tidyr','Epi','ggtext','lubridate','countrycode','sf','RColorBrewer',
              'patchwork','gtable','grid','broom','performance','forcats','readr','purrr','pwr')
lapply(packages, require, character.only = TRUE)

country_name_map <- c(
  "Taiwan" = "Taiwan (Province of China)",
  "South Korea" = "Republic of Korea",
  "United States" = "United States of America",
  "Vietnam" = "Viet Nam",
  "Cote d'Ivoire" = "Côte d'Ivoire",
  "Tanzania" = "United Republic of Tanzania",
  "Czech Republic" = "Czechia",
  "Laos" = "Lao People's Democratic Republic",
  "Moldova" = "Republic of Moldova",
  "Russia" = "Russian Federation",
  "Kyrgyz Republic" = "Kyrgyzstan",
  "Venezuela" = "Venezuela (Bolivarian Republic of)",
  "Slovak Republic" = "Slovakia",
  "Syria" = "Syrian Arab Republic",
  "Cape Verde" = "Cabo Verde",
  "Turkey" = "Türkiye",
  "Iran" = "Iran (Islamic Republic of)",
  "Brunei" = "Brunei Darussalam",
  "Bolivia" = "Bolivia (Plurinational State of)",
  "Democratic Republic of Congo" = "Democratic Republic of the Congo")

interaction_pairs_tbl <- tibble::tribble(
  ~x,                     ~z,
  "mental",               "ContainmentIndex",
  "mental",               "EconomicSupportIndex",
  "mental",               "HealthSystemIndex",
  "mental",               "VaccinationIndex",
  "mental",               "BorderIndex",
  "mental",               "FacemaskIndex",
  "ContainmentIndex",     "EconomicSupportIndex",
  "ContainmentIndex",     "HealthSystemIndex",
  "ContainmentIndex",     "VaccinationIndex",
  "ContainmentIndex",     "BorderIndex",
  "ContainmentIndex",     "FacemaskIndex",
  "EconomicSupportIndex", "HealthSystemIndex",
  "EconomicSupportIndex", "VaccinationIndex",
  "EconomicSupportIndex", "BorderIndex",
  "EconomicSupportIndex", "FacemaskIndex",
  "HealthSystemIndex",    "VaccinationIndex",
  "HealthSystemIndex",    "BorderIndex",
  "HealthSystemIndex",    "FacemaskIndex",
  "VaccinationIndex",     "BorderIndex",
  "VaccinationIndex",     "FacemaskIndex",
  "BorderIndex",          "FacemaskIndex")

joint_vars_to_scale <- c(
  "ave_ext",
  "mental",
  "ContainmentIndex",
  "EconomicSupportIndex",
  "HealthSystemIndex",
  "VaccinationIndex",
  "BorderIndex",
  "FacemaskIndex")

standardize_location <- function(x, name_map = country_name_map) {
  x <- as.character(x)
  ifelse(x %in% names(name_map), unname(name_map[x]), x)}

prepare_interaction_data <- function(target_year,
                                     name_map = country_name_map,
                                     policy_file = "policy_index_six.csv",
                                     excess_file = "F_country_prediction_with_excess_select_noshuffle_seed2026.csv",
                                     country_file = "gbd_country_id.csv",
                                     mental_file = "IHME-GBD_2023_DATA-2b3e5b00-1.csv",
                                     covid_file = "IHME-GBD_2023_DATA-2dadd859-1.csv") {
  
  # policy data
  policy_df <- read.csv(policy_file) %>%
    mutate(
      Date = as.Date(as.character(Date), format = "%Y%m%d"),
      year = year(Date)) %>%
    filter(year %in% target_year) %>%
    group_by(location_name) %>%
    summarise(
      ContainmentIndex = mean(ContainmentIndex, na.rm = TRUE),
      EconomicSupportIndex = mean(EconomicSupportIndex, na.rm = TRUE),
      HealthSystemIndex = mean(HealthSystemIndex, na.rm = TRUE),
      VaccinationIndex = mean(VaccinationIndex, na.rm = TRUE),
      BorderIndex = mean(BorderIndex, na.rm = TRUE),
      FacemaskIndex = mean(FacemaskIndex, na.rm = TRUE),
      .groups = "drop") %>%
    mutate(location_name = standardize_location(location_name, name_map))
  
  # excess burden
  country_id <- read.csv(country_file)
  excess_df <- read.csv(excess_file) %>%
    select(
      location_name = location,
      year,
      ave_ext = excess) %>%
    mutate(
      location_name = ifelse(location_name == "Taiwan",
                             "Taiwan (Province of China)",
                             location_name)) %>%
    left_join(country_id, by = "location_name") %>%
    filter(year %in% target_year) %>%
    group_by(location_name) %>%
    summarise(
      ave_ext = sum(ave_ext, na.rm = TRUE),
      .groups = "drop")
  
  # baseline mental burden
  mental_df <- read.csv(mental_file) %>%
    filter(year < 2020 & year > 2015, cause_name == "Mental disorders") %>%
    select(location_name, val) %>%
    mutate(
      location_name = ifelse(location_name == "Taiwan",
                             "Taiwan (Province of China)",
                             location_name)) %>%
    group_by(location_name) %>%
    summarise(
      mental = mean(val, na.rm = TRUE),
      .groups = "drop")
  
  # covid deaths
  covid_df <- read.csv(covid_file) %>%
    filter(year %in% target_year, measure_name == "Deaths") %>%
    select(location_name, val) %>%
    mutate(
      location_name = ifelse(location_name == "Taiwan",
                             "Taiwan (Province of China)",
                             location_name)) %>%
    group_by(location_name) %>%
    summarise(
      covid = mean(val, na.rm = TRUE),
      .groups = "drop") %>%
    mutate(location_name = standardize_location(location_name, name_map))
  
  mental_df %>%
    left_join(excess_df, by = "location_name") %>%
    left_join(covid_df, by = "location_name") %>%
    left_join(policy_df, by = "location_name") %>%
    na.omit()
}

fit_interaction_models <- function(dat3,
                                   target_year,
                                   interaction_pairs = interaction_pairs_tbl,
                                   vars_to_scale = joint_vars_to_scale) {
  
  dat3_std <- dat3 %>%
    mutate(across(all_of(vars_to_scale), ~ as.numeric(scale(.x))))
  
  base_rhs <- c(
    "mental",
    "ContainmentIndex",
    "EconomicSupportIndex",
    "HealthSystemIndex",
    "VaccinationIndex",
    "BorderIndex",
    "FacemaskIndex")
  
  base_formula <- as.formula(
    paste("ave_ext ~", paste(base_rhs, collapse = " + ")))
  
  reduced_model <- lm(base_formula, data = dat3_std)
  r2_reduced <- summary(reduced_model)$r.squared
  
  year_label <- if (length(target_year) == 1) {
    as.character(target_year)
  } else {
    paste0(min(target_year), "-", max(target_year))}
  
  classify_f2 <- function(f2) {
    case_when(
      is.na(f2)      ~ NA_character_,
      f2 < 0.02      ~ "very small",
      f2 < 0.15      ~ "small",
      f2 < 0.35      ~ "medium",
      TRUE           ~ "large")}
  
  # test interaction one by one
  fit_one_interaction <- function(x, z) {
    fml <- as.formula(
      paste(
        "ave_ext ~",
        paste(base_rhs, collapse = " + "),
        "+",
        paste0(x, ":", z)))
    
    full_model <- lm(fml, data = dat3_std)
    r2_full <- summary(full_model)$r.squared
    
    partial_R2 <- (r2_full - r2_reduced) / (1 - r2_reduced)
    f2 <- partial_R2 / (1 - partial_R2)
    
    list(
      model = full_model,
      partial_R2 = partial_R2,
      f2 = f2,
      effect_size = classify_f2(f2))}
  
  models_int <- pmap(interaction_pairs, fit_one_interaction)
  model_names <- paste0(interaction_pairs$x, " × ", interaction_pairs$z)
  names(models_int) <- model_names
  
  extract_coef <- function(mod_obj, x, z, label) {
    mod <- mod_obj$model
    partial_R2 <- mod_obj$partial_R2
    f2 <- mod_obj$f2
    effect_size <- mod_obj$effect_size
    
    term1 <- paste0(x, ":", z)
    term2 <- paste0(z, ":", x)
    
    out <- broom::tidy(mod, conf.int = TRUE) %>%
      filter(term %in% c(term1, term2)) %>%
      transmute(
        year = year_label,
        interaction = label,
        term,
        beta = estimate,
        std.error,
        p.value,
        conf.low,
        conf.high,
        partial_R2 = partial_R2,
        f2 = f2,
        effect_size = effect_size)
    
    if (nrow(out) == 0) {
      out <- tibble(
        year = year_label,
        interaction = label,
        term = term1,
        beta = NA_real_,
        std.error = NA_real_,
        p.value = NA_real_,
        conf.low = NA_real_,
        conf.high = NA_real_,
        partial_R2 = partial_R2,
        f2 = f2,
        effect_size = effect_size)}
    
    out %>%
      mutate(
        p = ifelse(
          is.na(p.value), NA_character_,
          ifelse(p.value < 0.001, "<0.001", formatC(p.value, digits = 3, format = "f")))) %>%
      select(year, interaction, term, beta, std.error, p.value, p,
             conf.low, conf.high, partial_R2, f2, effect_size)}
  
  coef_table <- pmap_dfr(
    list(models_int, interaction_pairs$x, interaction_pairs$z, model_names),
    extract_coef)
  
  get_vif_table <- function(mod_obj, label) {
    mod <- mod_obj$model
    check_collinearity(mod) %>%
      as.data.frame() %>%
      transmute(
        year = year_label,
        interaction = label,
        term = Term,
        VIF)}
  
  vif_table <- map2_dfr(models_int, model_names, get_vif_table)
  
  partial_r2_table <- tibble(
    year = year_label,
    interaction = model_names,
    partial_R2 = map_dbl(models_int, "partial_R2"),
    f2 = map_dbl(models_int, "f2"),
    effect_size = map_chr(models_int, "effect_size"))
  
  list(
    data = dat3,
    data_std = dat3_std,
    reduced_model = reduced_model,
    coef_table = coef_table,
    vif_table = vif_table,
    partial_r2_table = partial_r2_table,
    models = models_int)
}

run_interaction_analysis <- function(target_year,
                                     name_map = country_name_map,
                                     interaction_pairs = interaction_pairs_tbl,
                                     vars_to_scale = joint_vars_to_scale) {
  
  dat3 <- prepare_interaction_data(
    target_year = target_year,
    name_map = name_map)
  
  fit_interaction_models(
    dat3 = dat3,
    target_year = target_year,
    interaction_pairs = interaction_pairs,
    vars_to_scale = vars_to_scale)
}

year_sets <- list(
  2020,
  2021,
  2022)

all_int_results <- lapply(year_sets, run_interaction_analysis)

coef_int_all <- bind_rows(lapply(all_int_results, \(x) x$coef_table))
vif_int_all  <- bind_rows(lapply(all_int_results, \(x) x$vif_table))
partialR2_int_all <- bind_rows(lapply(all_int_results, \(x) x$partial_r2_table))

coef_int_all <- coef_int_all %>%
  mutate(
    across(c(beta, std.error, conf.low, conf.high, partial_R2, f2), ~ round(.x, 3)),
    CI = sprintf("(%.3f, %.3f)", conf.low, conf.high)) %>%
  select(interaction, year, beta, CI, p.value, p, partial_R2, f2, effect_size)

coef_plot <- coef_int_all %>%
  filter(!is.na(p.value) & p.value < 0.05) %>%
  arrange(interaction, year)

coef_plot

###### a brief power analysis ######
n <- 179
k_full <- 8   # 7 main effects + 1 interaction
v <- n - k_full - 1   # error df

mdes <- pwr.f2.test(
  u = 1,          # one interaction term
  v = v,
  sig.level = 0.05,
  power = 0.80)

mdes

f2_min <- mdes$f2
partial_R2_min <- f2_min / (1 + f2_min)

f2_min
partial_R2_min


