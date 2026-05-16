######################################################################################################
# this file is for analyzing the associations between policy indicators and excess mental disorder burden 
# Pei Zhang
# 2026-04-30
######################################################################################################
### packages
packages <- c('foreign','tsModel','lmtest','Epi','splines','vcd','scales','dplyr','ggplot2','ggsci',
              'RColorBrewer','ggpubr','tidyr','Epi','ggtext','lubridate','countrycode','sf','RColorBrewer',
              'patchwork','gtable','grid','broom','performance','forcats','readr')
lapply(packages, require, character.only=T) 

country_id <- read.csv("gbd_country_id.csv") 

results_year <- read.csv("F_country_prediction_with_excess_select_noshuffle_seed2026.csv") %>% 
  select(location_name=location, cause_name=cause,
         year,ave_ext=excess, ave_extl=excess_lower,
         ave_exth=excess_upper) %>%
  left_join(country_id, by="location_name")

results_year$location_name[results_year$location_name == "Taiwan"] <- "Taiwan (Province of China)"

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

dat1 <- dat0 %>% select(location_name=CountryName,RegionName,Date,C1,C2,C3,C4,C5,C6,C7,C8,E1,E2,E3,E4,
                        H1,H2,H3,H4,H5,H6,H7,H8,V1,V2,V3,V4)
dat1$Date <- as.Date(as.character(dat1$Date), format="%Y%m%d")

policy_index_vars <- c(
  "C1","C2","C3","C4","C5","C6","C7","C8",
  "E1","E2","E3","E4",
  "H1","H2","H3","H4","H5","H6","H7","H8",
  "V1","V2","V3","V4")

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

prepare_policy_data <- function(dat1, year_values, index_vars, name_map) {
  dat1 %>%
    filter(year(Date) %in% year_values) %>%
    group_by(location_name) %>%
    summarise(
      across(all_of(index_vars), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop") %>%
    mutate(
      across(all_of(index_vars), ~ ifelse(is.infinite(.), NA, .)),
      location_name = as.character(location_name),
      location_name = ifelse(
        location_name %in% names(name_map),
        unname(name_map[location_name]),
        location_name))}

prepare_burden_data <- function(results_year,
                                covid_file = "IHME-GBD_2023_DATA-2dadd859-1.csv",
                                mental_file = "IHME-GBD_2023_DATA-2b3e5b00-1.csv",
                                year_values) {
  
  df_ext <- results_year %>%
    filter(year %in% year_values) %>%
    group_by(location_name) %>%
    summarise(
      ave_ext = sum(ave_ext, na.rm = TRUE),
      .groups = "drop")
  
  df_covid <- read.csv(covid_file) %>%
    filter(year %in% year_values, measure_name == "Deaths") %>%
    select(location_name, val) %>%
    group_by(location_name) %>%
    summarise(
      covid = mean(val, na.rm = TRUE),
      .groups = "drop")
  
  df_mental <- read.csv(mental_file) %>%
    filter(year < 2020 & year > 2015, cause_name == "Mental disorders") %>%
    select(location_name, val) %>%
    group_by(location_name) %>%
    summarise(
      mental = mean(val, na.rm = TRUE),
      .groups = "drop")
  
  df_mental %>%
    left_join(df_covid, by = "location_name") %>%
    left_join(df_ext, by = "location_name")
}

prepare_partial_data <- function(dat1,
                                 results_year,
                                 year_values,
                                 index_vars,
                                 name_map) {
  policy_df <- prepare_policy_data(
    dat1 = dat1,
    year_values = year_values,
    index_vars = index_vars,
    name_map = name_map)
  
  burden_df <- prepare_burden_data(
    results_year = results_year,
    year_values = year_values)
  
  dat3 <- burden_df %>%
    left_join(policy_df, by = "location_name")
  
  df_long <- dat3 %>%
    select(all_of(index_vars), ave_ext, covid, mental) %>%
    pivot_longer(
      cols = all_of(index_vars),
      names_to = "IndexType",
      values_to = "IndexValue")
  
  list(
    wide_data = dat3,
    long_data = df_long)
}

get_partial_table <- function(df, 
                              xvar = "IndexValue", 
                              yvar = "ave_ext",
                              facet = "IndexType",
                              controls = c("mental")) {
  
  stopifnot(all(c(xvar, yvar, facet, controls) %in% colnames(df)))
  
  z <- paste(controls, collapse = " + ")
  f_full <- as.formula(paste(yvar, "~", xvar, "+", z))
  
  stdz <- function(v) {
    if (!is.numeric(v)) return(v)
    s <- stats::sd(v, na.rm = TRUE)
    if (is.na(s) || s == 0) return(rep(NA_real_, length(v)))
    as.numeric(scale(v))}
  
  cor_table <- df %>%
    group_by(.data[[facet]]) %>%
    group_modify(~{
      dat <- .x %>%
        select(all_of(c(xvar, yvar, controls))) %>%
        drop_na()
      
      if (nrow(dat) < (length(controls) + 3) ||
          all(is.na(dat[[xvar]])) ||
          all(is.na(dat[[yvar]]))) {
        return(tibble(
          term = xvar,
          estimate = NA_real_,
          std.error = NA_real_,
          statistic = NA_real_,
          p.value = NA_real_,
          adj_r2 = NA_real_,
          n = nrow(dat)))}
      
      dat_std <- dat %>%
        mutate(across(all_of(c(xvar, yvar, controls)), stdz))
      
      if (all(is.na(dat_std[[xvar]])) || all(is.na(dat_std[[yvar]]))) {
        return(tibble(
          term = xvar,
          estimate = NA_real_,
          std.error = NA_real_,
          statistic = NA_real_,
          p.value = NA_real_,
          adj_r2 = NA_real_,
          n = nrow(dat)))}
      
      m <- lm(f_full, data = dat_std)
      k <- broom::tidy(m)
      kx <- k[k$term == xvar, , drop = FALSE]
      
      tibble(
        term = xvar,
        estimate = kx$estimate,
        std.error = kx$std.error,
        statistic = kx$statistic,
        p.value = kx$p.value,
        adj_r2 = summary(m)$adj.r.squared,
        n = stats::nobs(m))
    }) %>%
    ungroup() %>%
    rename(!!facet := 1) %>%
    mutate(
      beta = round(estimate, 3),
      se = round(std.error, 3),
      t = round(statistic, 3),
      p = round(p.value, 3),
      adj_r2 = round(adj_r2, 3)) %>%
    select(all_of(facet), beta, se, t, p, adj_r2, n)
  
  cor_table
}

run_partial_analysis_table <- function(dat1,
                                       results_year,
                                       year_values,
                                       controls = c("mental"),
                                       index_vars = policy_index_vars,
                                       name_map = country_name_map) {
  
  prepared <- prepare_partial_data(
    dat1 = dat1,
    results_year = results_year,
    year_values = year_values,
    index_vars = index_vars,
    name_map = name_map)
  
  year_label <- if (length(year_values) == 1) {
    as.character(year_values)
  } else {
    paste0(min(year_values), "-", max(year_values))}
  
  get_partial_table(
    df = prepared$long_data,
    controls = controls) %>%
    mutate(year = year_label)}

cor_table_2020 <- run_partial_analysis_table(dat1, results_year, 2020, controls = c("mental"))
cor_table_2021 <- run_partial_analysis_table(dat1, results_year, 2021, controls = c("mental"))
cor_table_2022 <- run_partial_analysis_table(dat1, results_year, 2022, controls = c("mental"))
cor_table_2020_2022 <- run_partial_analysis_table(dat1, results_year, 2020:2022, controls = c("mental"))

cor_table_all <- bind_rows(
  cor_table_2020,
  cor_table_2021,
  cor_table_2022,
  cor_table_2020_2022) %>% 
  mutate(sig=ifelse(p<0.05, 1, 0)) %>%
  select(IndexType, year, beta, p, sig)

plot_beta_heatmap <- function(data,
                              index_subset,
                              show_y_text = FALSE,
                              legend_title = "β",
                              year_levels = c("2022", "2021", "2020"),
                              fill_limits = c(-0.34, 0.34)) {
  
  df_plot <- data %>%
    filter(year %in% year_levels, IndexType %in% index_subset) %>%
    mutate(
      sig_label = ifelse(sig == 1, "*", ""),
      year = factor(year, levels = year_levels),
      IndexType = factor(IndexType, levels = index_subset))
  
  ggplot(df_plot, aes(x = IndexType, y = year, fill = beta)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sig_label), color = "#ED2023", size = 7, fontface = "bold") +
    scale_fill_gradient2(
      low = "#1E79B4",
      mid = "white",
      high = "#BB2027",
      midpoint = 0,
      limits = fill_limits,
      oob = squish,
      name = legend_title,
      na.value = "grey80") +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = "", y = "", title = "") +
    coord_fixed(ratio = 1) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      axis.text.y = if (show_y_text) element_text() else element_blank(),
      panel.grid = element_blank(),
      legend.position = "none",
      axis.ticks.x = element_line(color = "grey60", linewidth = 0.5),
      axis.ticks.y = element_line(color = "grey60", linewidth = 0.5),
      axis.line.x = element_line(color = "grey60"),
      axis.line.y = element_line(color = "grey60"),
      panel.border = element_rect(color = "grey60", fill = NA, linewidth = 0.5))
}

pp1 <- plot_beta_heatmap(
  data = cor_table_all,
  index_subset = c("C1","C2","C3","C4","C5","C6","C7","C8"),
  show_y_text = TRUE,
  legend_title = "β")

pp2 <- plot_beta_heatmap(
  data = cor_table_all,
  index_subset = c("E1","E2","E3","E4"),
  show_y_text = FALSE,
  legend_title = "β")

pp3 <- plot_beta_heatmap(
  data = cor_table_all,
  index_subset = c("H1","H2","H3","H4","H5","H6","H7","H8"),
  show_y_text = FALSE,
  legend_title = "β")

pp4 <- plot_beta_heatmap(
  data = cor_table_all,
  index_subset = c("V1","V2","V3","V4"),
  show_y_text = FALSE,
  legend_title = "Effect size (β)")

fig1 <- (pp1 | pp2 | pp3 | pp4) +
  plot_layout(widths = c(2, 1, 2, 1))

fig1


