######################################################################################################
# this file is for independent associations (mutually adjusted) between policy indices and excess mental disorder burden 
# Pei Zhang
# 2026-04-30
######################################################################################################
### packages
packages <- c('foreign','tsModel','lmtest','Epi','splines','vcd','scales','dplyr','ggplot2','ggsci',
              'RColorBrewer','ggpubr','tidyr','Epi','ggtext','lubridate','countrycode','sf','RColorBrewer',
              'patchwork','gtable','grid','broom','performance','forcats','readr')
lapply(packages, require, character.only=T)

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

joint_vars_to_scale <- c(
  "ave_ext",
  "mental",
  "covid",
  "ContainmentIndex",
  "EconomicSupportIndex",
  "HealthSystemIndex",
  "VaccinationIndex",
  "BorderIndex",
  "FacemaskIndex")

standardize_location <- function(x, name_map) {
  x <- as.character(x)
  ifelse(x %in% names(name_map), unname(name_map[x]), x)}

prepare_joint_policy_data <- function(target_year,
                                      name_map,
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
  results_year <- read.csv(excess_file) %>%
    select(
      location_name = location,
      cause_name = cause,
      year,
      ave_ext = excess,
      ave_extl = excess_lower,
      ave_exth = excess_upper) %>%
    mutate(
      location_name = ifelse(location_name == "Taiwan",
                             "Taiwan (Province of China)",
                             location_name)) %>%
    left_join(country_id, by = "location_name")

  excess_df <- results_year %>%
    filter(year %in% target_year) %>%
    group_by(location_name) %>%
    summarise(
      ave_ext = sum(ave_ext, na.rm = TRUE),
      .groups = "drop")

  # baseline mental burden
  mental_df <- read.csv(mental_file) %>%
    filter(year < 2020 & year > 2015, cause_name == "Mental disorders") %>%
    select(location_name, val) %>%
    group_by(location_name) %>%
    summarise(
      mental = mean(val, na.rm = TRUE),
      .groups = "drop")%>%
    mutate(location_name = standardize_location(location_name, name_map))

  # covid deaths
  covid_df <- read.csv(covid_file) %>%
    filter(year %in% target_year, measure_name == "Deaths") %>%
    select(location_name, val) %>%
    group_by(location_name) %>%
    summarise(
      covid = mean(val, na.rm = TRUE),
      .groups = "drop") %>%
    mutate(location_name = standardize_location(location_name, name_map))

  # final merged data
  dat3 <- mental_df %>%
    left_join(excess_df, by = "location_name") %>%
    left_join(covid_df, by = "location_name") %>%
    left_join(policy_df, by = "location_name") %>%
    na.omit()

  dat3
}

fit_joint_policy_model <- function(dat3,
                                   target_year,
                                   vars_to_scale = joint_vars_to_scale) {

  dat3_std <- dat3 %>%
    mutate(across(all_of(vars_to_scale), ~ as.numeric(scale(.x))))

  model_formula <- ave_ext ~ mental +
    ContainmentIndex +
    EconomicSupportIndex +
    HealthSystemIndex +
    VaccinationIndex +
    BorderIndex +
    FacemaskIndex

  m_joint <- lm(model_formula, data = dat3_std)

  year_label <- if (length(target_year) == 1) {
    as.character(target_year)
  } else {
    paste0(min(target_year), "-", max(target_year))}

  vif_table <- check_collinearity(m_joint) %>%
    as.data.frame() %>%
    transmute(
      year = year_label,
      term = Term,
      VIF)

  coef_table <- broom::tidy(m_joint, conf.int = TRUE) %>%
    filter(term != "(Intercept)") %>%
    transmute(
      year = year_label,
      term,
      beta = estimate,
      std.error,
      p.value,
      conf.low,
      conf.high) %>%
    mutate(
      p = ifelse(p.value < 0.001, "<0.001", formatC(p.value, digits = 3, format = "f"))) %>%
    select(year, term, beta, std.error, p, conf.low, conf.high)

  coef_table_formatted <- coef_table %>%
    mutate(
      beta = formatC(as.numeric(beta), digits = 2, format = "f"),
      CI = sprintf("(%s, %s)",
                     formatC(conf.low, digits = 2, format = "f"),
                     formatC(conf.high, digits = 2, format = "f")),
      p = ifelse(
        p %in% c("<0.001", "NA"),
        p,
        formatC(as.numeric(p), digits = 3, format = "f"))) %>%
    select(year, term, beta, CI, p)

  list(
    data = dat3,
    data_std = dat3_std,
    model = m_joint,
    vif_table = vif_table,
    coef_table = coef_table,
    coef_table_formatted = coef_table_formatted)
}

run_joint_policy_analysis <- function(target_year,
                                      name_map = country_name_map,
                                      vars_to_scale = joint_vars_to_scale) {

  dat3 <- prepare_joint_policy_data(
    target_year = target_year,
    name_map = name_map)

  fit_joint_policy_model(
    dat3 = dat3,
    target_year = target_year,
    vars_to_scale = vars_to_scale)
}

years <- list(2020, 2021, 2022, 2020:2022)

all_results <- lapply(years, run_joint_policy_analysis)

vif_all <- bind_rows(lapply(all_results, \(x) x$vif_table))
coef_all <- bind_rows(lapply(all_results, \(x) x$coef_table))
coef_all_fmt <- bind_rows(lapply(all_results, \(x) x$coef_table_formatted))

vif_all

### prepare plotting data
plot_df <- coef_all %>%
  filter(
    year %in% c("2020", "2021", "2022"),
    !term %in% c("mental")) %>% 
  mutate(
    beta = as.numeric(beta),
    SizeL = as.numeric(conf.low),
    SizeH = as.numeric(conf.high),
    p_num = parse_number(as.character(p)),
    sig = if_else(!is.na(p_num) & p_num <= 0.05, 1L, 0L),
    year = factor(year, levels = c("2020", "2021", "2022")),
    term = factor(
      term,
      levels = c(
        "ContainmentIndex",
        "BorderIndex",
        "EconomicSupportIndex",
        "HealthSystemIndex",
        "FacemaskIndex",
        "VaccinationIndex"),
      labels = c(
        "Redefined Containment\nIndex",
        "Border control Index",
        "Economic Support\nIndex",
        "Redefined Health-System\nIndex",
        "Facial Covering\nIndex",
        "Vaccination Index")))

mytheme <- theme_bw() +
  theme(
    panel.border = element_blank(),
    panel.background = element_rect(fill = NA, colour = "grey90"),
    strip.background = element_rect(fill = NA, colour = NA),
    strip.text.x = element_text(size = 14, color = "gray10", face = "bold", family = "serif"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(linewidth = 0.5, linetype = "solid", colour = "grey90"),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.key.width = unit(1.1, "cm"),
    legend.text = element_text(size = 12, color = "black", family = "serif"),
    legend.spacing.y = unit(0.1, "cm"),
    legend.background = element_rect(color = NA),
    legend.box.margin = margin(0, 0, 0, 0, "cm"),
    legend.margin = margin(0, 0, 0, 0, "cm"),
    plot.margin = margin(0.5, 0.5, 0.1, 0.5, "cm"),
    axis.title.x = element_text(size = 16, family = "serif"),
    axis.title.y = element_text(size = 16, family = "serif"),
    axis.text.x = element_text(color = "black", size = 16, family = "serif",
                               angle = 45, hjust = 1, vjust = 1),
    axis.text.y = element_text(color = "black", size = 16, family = "serif"),
    plot.title = element_text(size = 20, hjust = 0.5, family = "serif"))

### plot
pd <- position_dodge(width = 0.6)

p_joint <- ggplot(plot_df, aes(x = term, group = interaction(term, year))) +
  geom_errorbar(
    aes(ymin = SizeL, ymax = SizeH),
    color = "gray",
    position = pd,
    width = 0,
    linewidth = 1,
    alpha = 0.8) +
  geom_point(
    aes(y = beta, fill = beta, shape = year),
    color = "gray",
    size = 4,
    stroke = 0.8,
    position = pd) +
  geom_hline(yintercept = 0, linetype = 2, color = "gray") +
  scale_y_continuous(
    limits = c(-0.5, 0.5),
    breaks = c(-0.4, 0, 0.4)) +
  scale_fill_gradientn(
    name = "Effect size",
    colors = rev(brewer.pal(11, "RdBu")), 
    limits = c(-0.4, 0.4),
    breaks = c(-0.4, 0, 0.4),
    labels = c("-0.4", "0", "0.4")) +
  scale_shape_manual(
    name = "Year",
    values = c(
      "2020" = 21,
      "2021" = 22,
      "2022" = 23,
      "2020-2022" = 24)) +
  labs(
    x = "",
    y = "Effect estimates") +
  mytheme +
  guides(
    shape = guide_legend(
      label.position = "bottom",
      title.position = "top",
      nrow = 1)) +
  theme(
    panel.border = element_blank(),
    axis.line.y = element_blank(),
    legend.key = element_blank(),
    legend.background = element_blank())

p_joint

