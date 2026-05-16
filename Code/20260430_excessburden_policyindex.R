######################################################################################################
# this file is for analyzing the associations between policy indices and excess mental disorder burden 
# Pei Zhang
# 2026-04-30
######################################################################################################
### packages
packages <- c('foreign','tsModel','lmtest','Epi','splines','vcd','scales','dplyr','ggplot2','ggsci',
              'RColorBrewer','ggpubr','tidyr','Epi','ggtext','countrycode','sf','RColorBrewer',
              'patchwork','gtable','grid','broom','performance','forcats','readr','lubridate',
              'sandwich')
lapply(packages, require, character.only=T) 

country_id <- read.csv("gbd_country_id.csv") 
results_year <- read.csv("F_country_prediction_with_excess_select_noshuffle_seed2026.csv") %>% 
  select(location_name=location, cause_name=cause,
         year,ave_ext=excess, ave_extl=excess_lower,
         ave_exth=excess_upper) %>%
  left_join(country_id, by="location_name")

results_year$location_name[results_year$location_name == "Taiwan"] <- "Taiwan (Province of China)"

dat1 <- read.csv("policy_indices_four.csv")
dat1$Date <- as.Date(as.character(dat1$Date), format="%Y%m%d")
dat1$location_name <- with(dat1, reorder(location_name, ContainmentIndex, FUN = mean))

harmonize_country_names <- function(x) {
  name_map <- c(
    "Taiwan"                      = "Taiwan (Province of China)",
    "South Korea"                 = "Republic of Korea",
    "United States"               = "United States of America",
    "Vietnam"                     = "Viet Nam",
    "Cote d'Ivoire"               = "Côte d'Ivoire",
    "Tanzania"                    = "United Republic of Tanzania",
    "Czech Republic"              = "Czechia",
    "Laos"                        = "Lao People's Democratic Republic",
    "Moldova"                     = "Republic of Moldova",
    "Russia"                      = "Russian Federation",
    "Kyrgyz Republic"             = "Kyrgyzstan",
    "Venezuela"                   = "Venezuela (Bolivarian Republic of)",
    "Slovak Republic"             = "Slovakia",
    "Syria"                       = "Syrian Arab Republic",
    "Cape Verde"                  = "Cabo Verde",
    "Turkey"                      = "Türkiye",
    "Iran"                        = "Iran (Islamic Republic of)",
    "Brunei"                      = "Brunei Darussalam",
    "Bolivia"                     = "Bolivia (Plurinational State of)",
    "Democratic Republic of Congo"= "Democratic Republic of the Congo")
  
  x <- as.character(x)
  ifelse(x %in% names(name_map), unname(name_map[x]), x)
}

prepare_partial_policy_data <- function(dat1,
                                        results_year,
                                        covid_file = "IHME-GBD_2023_DATA-2dadd859-1.csv",
                                        mental_file = "IHME-GBD_2023_DATA-2b3e5b00-1.csv",
                                        year_values = 2022,
                                        policy_vars = c(
                                          "ContainmentIndex",
                                          "EconomicSupportIndex",
                                          "HealthSystemIndex",
                                          "VaccinationIndex"
                                        )) {
  # policy data
  policy_df <- dat1 %>%
    filter(year(Date) %in% year_values) %>%
    group_by(location_name) %>%
    summarise(
      across(all_of(policy_vars), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop") %>%
    mutate(location_name = harmonize_country_names(location_name))
  
  # burden data
  burden_df <- results_year %>%
    filter(year %in% year_values) %>%
    group_by(location_name) %>%
    summarise(
      ave_ext = mean(ave_ext, na.rm = TRUE),
      .groups = "drop") %>%
    mutate(location_name = harmonize_country_names(location_name))
  
  # covid deaths
  covid_df <- read.csv(covid_file) %>%
    filter(year %in% year_values, measure_name == "Deaths") %>%
    group_by(location_name) %>%
    summarise(
      covid = mean(val, na.rm = TRUE),
      .groups = "drop") %>%
    mutate(location_name = harmonize_country_names(location_name))
  
  # baseline mental burden
  mental_df <- read.csv(mental_file) %>%
    filter(year < 2020 & year > 2015, cause_name == "Mental disorders") %>%
    group_by(location_name) %>%
    summarise(
      mental = mean(val, na.rm = TRUE),
      .groups = "drop") %>%
    mutate(location_name = harmonize_country_names(location_name))
  
  # merge
  merged_df <- burden_df %>%
    left_join(covid_df, by = "location_name") %>%
    left_join(mental_df, by = "location_name") %>%
    left_join(policy_df, by = "location_name") %>%
    drop_na()
  
  merged_df <- merged_df %>%
    mutate(
      ContainmentIndex_order     = order(ContainmentIndex),
      EconomicSupportIndex_order = order(EconomicSupportIndex),
      HealthSystemIndex_order    = order(HealthSystemIndex),
      VaccinationIndex_order     = order(VaccinationIndex),
      ext_order                  = order(ave_ext))
  
  # long format
  long_df <- merged_df %>%
    select(all_of(policy_vars), ave_ext, covid, mental) %>%
    pivot_longer(
      cols = all_of(policy_vars),
      names_to = "IndexType",
      values_to = "IndexValue")
  
  long_df$IndexType <- factor(
    long_df$IndexType,
    levels = c(
      "ContainmentIndex",
      "EconomicSupportIndex",
      "HealthSystemIndex",
      "VaccinationIndex"),
    labels = c(
      "Containment Index",
      "Economic Support Index",
      "Health System Index",
      "Vaccination Index"))
  
  list(
    wide_data = merged_df,
    long_data = long_df)
}

plot_partial_index <- function(df,
                               xvar = "IndexValue",
                               yvar = "ave_ext",
                               facet = "IndexType",
                               controls = c("mental"),
                               nrow = 1,
                               scales = "free_x",
                               point_color = "#DDAE7E",
                               line_color = "#2F5275") {
  
  stopifnot(all(c(xvar, yvar, facet, controls) %in% colnames(df)))
  
  # formulas
  z <- paste(controls, collapse = " + ")
  f_y_on_Z <- as.formula(paste(yvar, "~", z))
  f_x_on_Z <- as.formula(paste(xvar, "~", z))
  f_full <- as.formula(paste(yvar, "~", xvar, "+", z))
  
  # standardize helper
  stdz <- function(v) {
    if (!is.numeric(v)) return(v)
    s <- stats::sd(v, na.rm = TRUE)
    if (is.na(s) || s == 0) return(rep(NA_real_, length(v)))
    as.numeric(scale(v))}
  
  # partial residual data
  df_partial <- df %>%
    group_by(.data[[facet]]) %>%
    group_modify(~{
      dat <- .x %>% drop_na(all_of(c(xvar, yvar, controls)))
      if (nrow(dat) < (length(controls) + 3)) {
        return(dat %>% mutate(res_x = NA_real_, res_y = NA_real_))}
      dat$res_y <- resid(lm(f_y_on_Z, data = dat, na.action = na.exclude))
      dat$res_x <- resid(lm(f_x_on_Z, data = dat, na.action = na.exclude))
      dat}) %>%
    ungroup()
  
  # regression summary table
  cor_table <- df %>%
    group_by(.data[[facet]]) %>%
    group_modify(~{
      need <- c(xvar, yvar, controls)
      dat <- .x %>% select(all_of(need)) %>% drop_na()
      
      if (nrow(dat) < (length(controls) + 3)) {
        tibble(
          term = xvar,
          estimate = NA_real_,
          std.error = NA_real_,
          statistic = NA_real_,
          p.value = NA_real_,
          adj_r2 = NA_real_,
          n = nrow(dat))
      } else {
        dat_std <- dat %>% mutate(across(all_of(need), stdz))
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
          n = stats::nobs(m))}
    }) %>%
    ungroup() %>%
    rename(!!facet := 1) %>%
    mutate(
      beta = formatC(estimate, digits = 3, format = "f"),
      p_txt = case_when(
        is.na(p.value)  ~ NA_character_,
        p.value < 0.001 ~ "<0.001",
        TRUE ~ formatC(p.value, digits = 3, format = "f")),
      se = formatC(std.error, digits = 3, format = "f"),
      t = formatC(statistic, digits = 3, format = "f"),
      adj_r2 = formatC(adj_r2, digits = 3, format = "f"),
      sig_star = if_else(!is.na(p.value) & p.value < 0.05, "*", ""),
      label = if_else(
        !is.na(p.value) & p.value < 0.001,
        sprintf('italic(beta)~"="~"%s,"~~italic(p)~"<"~"0.001"', beta),
        sprintf('italic(beta)~"="~"%s,"~~italic(p)~"="~"%s"', beta, p_txt))) %>%
    select(all_of(facet), beta, se, t, p_txt, adj_r2, n, p.value, sig_star, label)
  
  p <- ggplot(df_partial, aes(x = res_x, y = res_y)) +
    geom_point(size = 3, alpha = 0.5, color = point_color) +
    geom_smooth(method = "lm", color = line_color, se = TRUE, linewidth = 1.2) +
    # geom_smooth(method = "loess", color = "#2F5275", se = TRUE, size = 1.2) +
    # geom_smooth(method = "gam", formula = y ~ s(x, k = 4),
    #             color = "#2F5275", se = TRUE, size = 1.2) +
    facet_wrap(
      reformulate(termlabels = facet),
      scales = scales,
      nrow = nrow,
      strip.position = "top") +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.15))) +
    geom_text(
      data = cor_table,
      aes(x = -Inf, y = Inf, label = label),
      inherit.aes = FALSE,
      hjust = -0.1, vjust = 1.6,
      size = 3,
      parse = TRUE,
      color = "black") +
    geom_text(
      data = cor_table %>% filter(!is.na(p.value) & p.value < 0.05),
      aes(x = Inf, y = Inf, label = "*"),
      inherit.aes = FALSE,
      hjust = 4.0, vjust = 1.3,
      size = 7,
      fontface = "bold",
      color = "red") +
    labs(x = "", y = "") +
    theme_classic(base_size = 14) +
    theme(
      strip.background = element_blank(),
      strip.text = element_blank(),
      strip.placement = "outside",
      panel.spacing = unit(1, "lines"),
      panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.5),
      axis.title.y = element_text(size = 10),
      axis.text.y.left = element_text(size = 10),
      axis.ticks.y.left = element_line(color = "grey60", linewidth = 0.5),
      axis.line = element_line(color = "grey70", linewidth = 0.5),
      axis.text.y.right = element_blank(),
      axis.ticks.y.right = element_blank())
  
  list(
    plot = p,
    cor_table = cor_table,
    partial_data = df_partial)
}

make_partial_policy_plot <- function(dat1,
                                     results_year,
                                     year_values = 2022,
                                     controls = c("mental"),
                                     nrow = 1,
                                     covid_file = "IHME-GBD_2023_DATA-2dadd859-1.csv",
                                     mental_file = "IHME-GBD_2023_DATA-2b3e5b00-1.csv") {
  
  dat_list <- prepare_partial_policy_data(
    dat1 = dat1,
    results_year = results_year,
    covid_file = covid_file,
    mental_file = mental_file,
    year_values = year_values)
  
  res <- plot_partial_index(
    df = dat_list$long_data,
    controls = controls,
    nrow = nrow)
  
  list(
    wide_data = dat_list$wide_data,
    long_data = dat_list$long_data,
    plot = res$plot,
    cor_table = res$cor_table,
    partial_data = res$partial_data)
}

# pooled 2020-2022 function
prepare_pooled_policy_data <- function(dat1,
                                       results_year,
                                       covid_file = "IHME-GBD_2023_DATA-2dadd859-1.csv",
                                       mental_file = "IHME-GBD_2023_DATA-2b3e5b00-1.csv",
                                       year_values = 2020:2022,
                                       policy_vars = c(
                                         "ContainmentIndex",
                                         "EconomicSupportIndex",
                                         "HealthSystemIndex",
                                         "VaccinationIndex"
                                       )) {
  # policy data: country-year level
  policy_df <- dat1 %>%
    mutate(year = lubridate::year(Date)) %>%
    filter(year %in% year_values) %>%
    group_by(location_name, year) %>%
    summarise(
      across(all_of(policy_vars), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop") %>%
    mutate(location_name = harmonize_country_names(location_name))
  
  # burden data: country-year level
  burden_df <- results_year %>%
    filter(year %in% year_values) %>%
    group_by(location_name, year) %>%
    summarise(
      ave_ext = mean(ave_ext, na.rm = TRUE),
      .groups = "drop") %>%
    mutate(location_name = harmonize_country_names(location_name))
  
  # covid deaths: country-year level
  covid_df <- read.csv(covid_file) %>%
    filter(year %in% year_values, measure_name == "Deaths") %>%
    group_by(location_name, year) %>%
    summarise(
      covid = mean(val, na.rm = TRUE),
      .groups = "drop") %>%
    mutate(location_name = harmonize_country_names(location_name))
  
  # baseline mental burden: country level
  mental_df <- read.csv(mental_file) %>%
    filter(year < 2020 & year > 2015, cause_name == "Mental disorders") %>%
    group_by(location_name) %>%
    summarise(
      mental = mean(val, na.rm = TRUE),
      .groups = "drop") %>%
    mutate(location_name = harmonize_country_names(location_name))
  
  # merge
  merged_df <- burden_df %>%
    left_join(covid_df, by = c("location_name", "year")) %>%
    left_join(mental_df, by = "location_name") %>%
    left_join(policy_df, by = c("location_name", "year")) %>%
    drop_na() %>%
    mutate(
      year_f = factor(year),
      ContainmentIndex_order     = order(ContainmentIndex),
      EconomicSupportIndex_order = order(EconomicSupportIndex),
      HealthSystemIndex_order    = order(HealthSystemIndex),
      VaccinationIndex_order     = order(VaccinationIndex),
      ext_order                  = order(ave_ext))
  
  # long format
  long_df <- merged_df %>%
    select(location_name, year, year_f, all_of(policy_vars), ave_ext, covid, mental) %>%
    pivot_longer(
      cols = all_of(policy_vars),
      names_to = "IndexType",
      values_to = "IndexValue")
  
  long_df$IndexType <- factor(
    long_df$IndexType,
    levels = c(
      "ContainmentIndex",
      "EconomicSupportIndex",
      "HealthSystemIndex",
      "VaccinationIndex"),
    labels = c(
      "Containment Index",
      "Economic Support Index",
      "Health System Index",
      "Vaccination Index"))
  
  list(
    wide_data = merged_df,
    long_data = long_df)
}

plot_partial_index_clustered <- function(df,
                                         xvar = "IndexValue",
                                         yvar = "ave_ext",
                                         facet = "IndexType",
                                         controls = c("mental", "year_f"),
                                         cluster_var = "location_name",
                                         nrow = 1,
                                         scales = "free_x",
                                         point_color = "#DDAE7E",
                                         line_color = "#2F5275",
                                         show_se = FALSE) {
  
  stopifnot(all(c(xvar, yvar, facet, controls, cluster_var) %in% colnames(df)))
  
  z <- paste(controls, collapse = " + ")
  f_y_on_Z <- as.formula(paste(yvar, "~", z))
  f_x_on_Z <- as.formula(paste(xvar, "~", z))
  f_full <- as.formula(paste(yvar, "~", xvar, "+", z))
  
  stdz <- function(v) {
    if (!is.numeric(v)) return(v)
    s <- stats::sd(v, na.rm = TRUE)
    if (is.na(s) || s == 0) return(rep(NA_real_, length(v)))
    as.numeric(scale(v))}
  
  # partial residuals
  df_partial <- df %>%
    group_by(.data[[facet]]) %>%
    group_modify(~{
      dat <- .x %>% drop_na(all_of(c(xvar, yvar, controls, cluster_var)))
      if (nrow(dat) < (length(controls) + 5)) {
        return(dat %>% mutate(res_x = NA_real_, res_y = NA_real_))}
      dat$res_y <- resid(lm(f_y_on_Z, data = dat, na.action = na.exclude))
      dat$res_x <- resid(lm(f_x_on_Z, data = dat, na.action = na.exclude))
      dat}) %>%
    ungroup()
  
  # regression summary with clustered SE
  cor_table <- df %>%
    group_by(.data[[facet]]) %>%
    group_modify(~{
      need <- c(xvar, yvar, controls, cluster_var)
      dat <- .x %>% select(all_of(need)) %>% drop_na()
      
      if (nrow(dat) < (length(controls) + 5)) {
        return(tibble(
          estimate = NA_real_,
          std.error = NA_real_,
          statistic = NA_real_,
          p.value = NA_real_,
          adj_r2 = NA_real_,
          n = nrow(dat),
          n_cluster = dplyr::n_distinct(dat[[cluster_var]])))}
      
      dat_std <- dat %>%
        mutate(across(where(is.numeric), stdz))
      
      m <- lm(f_full, data = dat_std)
      
      V_cl <- sandwich::vcovCL(m, cluster = dat_std[[cluster_var]], type = "HC1")
      ct <- lmtest::coeftest(m, vcov. = V_cl)
      
      p_col <- grep("^Pr\\(", colnames(ct), value = TRUE)
      stat_col <- setdiff(colnames(ct), c("Estimate", "Std. Error", p_col))[1]
      
      tibble(
        estimate = unname(ct[xvar, "Estimate"]),
        std.error = unname(ct[xvar, "Std. Error"]),
        statistic = unname(ct[xvar, stat_col]),
        p.value = unname(ct[xvar, p_col]),
        adj_r2 = summary(m)$adj.r.squared,
        n = stats::nobs(m),
        n_cluster = dplyr::n_distinct(dat[[cluster_var]]))}) %>%
    ungroup() %>%
    rename(!!facet := 1) %>%
    mutate(
      beta = formatC(estimate, digits = 3, format = "f"),
      p_txt = case_when(
        is.na(p.value)  ~ NA_character_,
        p.value < 0.001 ~ "<0.001",
        TRUE            ~ formatC(p.value, digits = 3, format = "f")),
      se = formatC(std.error, digits = 3, format = "f"),
      stat = formatC(statistic, digits = 3, format = "f"),
      adj_r2 = formatC(adj_r2, digits = 3, format = "f"),
      sig_star = if_else(!is.na(p.value) & p.value < 0.05, "*", ""),
      label = if_else(
        !is.na(p.value) & p.value < 0.001,
        sprintf('italic(beta)~"="~"%s,"~~italic(p)~"<"~"0.001"', beta),
        sprintf('italic(beta)~"="~"%s,"~~italic(p)~"="~"%s"', beta, p_txt))) %>%
    select(all_of(facet), beta, se, stat, p_txt, adj_r2, n, n_cluster, p.value, sig_star, label)
  
  p <- ggplot(df_partial, aes(x = res_x, y = res_y)) +
    geom_point(size = 3, alpha = 0.5, color = point_color) +
    geom_smooth(method = "lm", color = line_color, se = show_se, linewidth = 1.2) +
    facet_wrap(
      reformulate(termlabels = facet),
      scales = scales,
      nrow = nrow,
      strip.position = "top") +
    scale_y_continuous(
      expand = expansion(mult = c(0.05, 0.15))) +
    geom_text(
      data = cor_table,
      aes(x = -Inf, y = Inf, label = label),
      inherit.aes = FALSE,
      hjust = -0.1, vjust = 1.6,
      size = 3,
      parse = TRUE,
      color = "black") +
    geom_text(
      data = cor_table %>% filter(!is.na(p.value) & p.value < 0.05),
      aes(x = Inf, y = Inf, label = "*"),
      inherit.aes = FALSE,
      hjust = 4.0, vjust = 1.3,
      size = 7,
      fontface = "bold",
      color = "red") +
    labs(x = "", y = "") +
    theme_classic(base_size = 14) +
    theme(
      strip.background = element_blank(),
      strip.text = element_blank(),
      strip.placement = "outside",
      panel.spacing = unit(1, "lines"),
      panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.5),
      axis.title.y = element_text(size = 10),
      axis.text.y.left = element_text(size = 10),
      axis.ticks.y.left = element_line(color = "grey60", linewidth = 0.5),
      axis.line = element_line(color = "grey70", linewidth = 0.5),
      axis.text.y.right = element_blank(),
      axis.ticks.y.right = element_blank())
  
  list(
    plot = p,
    cor_table = cor_table,
    partial_data = df_partial)
}

make_pooled_policy_plot <- function(dat1,
                                    results_year,
                                    year_values = 2020:2022,
                                    controls = c("mental", "year_f"),
                                    nrow = 1,
                                    covid_file = "IHME-GBD_2023_DATA-2dadd859-1.csv",
                                    mental_file = "IHME-GBD_2023_DATA-2b3e5b00-1.csv") {
  
  dat_list <- prepare_pooled_policy_data(
    dat1 = dat1,
    results_year = results_year,
    covid_file = covid_file,
    mental_file = mental_file,
    year_values = year_values)
  
  res <- plot_partial_index_clustered(
    df = dat_list$long_data,
    controls = controls,
    cluster_var = "location_name",
    nrow = nrow,
    show_se = T)
  
  list(
    wide_data = dat_list$wide_data,
    long_data = dat_list$long_data,
    plot = res$plot,
    cor_table = res$cor_table,
    partial_data = res$partial_data)
}


res_2020 <- make_partial_policy_plot(
  dat1 = dat1,
  results_year = results_year,
  year_values = 2020,
  controls = c("mental"))
p_2020 <- res_2020$plot
cor_2020 <- res_2020$cor_table

res_2021 <- make_partial_policy_plot(
  dat1 = dat1,
  results_year = results_year,
  year_values = 2021,
  controls = c("mental")
)
p_2021 <- res_2021$plot
cor_2021 <- res_2021$cor_table

res_2022 <- make_partial_policy_plot(
  dat1 = dat1,
  results_year = results_year,
  year_values = 2022,
  controls = c("mental")
)
p_2022 <- res_2022$plot
cor_2022 <- res_2022$cor_table

res_2020_2022 <- make_pooled_policy_plot(
  dat1 = dat1,
  results_year = results_year,
  year_values = 2020:2022,
  controls = c("mental", "year_f")
)
p_2020_2022 <- res_2020_2022$plot
cor_2020_2022 <- res_2020_2022$cor_table

fig1 <- ggarrange(p_2020_2022,p_2020,p_2021,p_2022,nrow = 4,ncol = 1,labels = c("All","2020","2021","2022"),
                  font.label = list(size = 10))

header_row <- ggarrange(
  text_grob("Containment Index", face = "bold", size = 10),
  text_grob("Economic Support Index", face = "bold", size = 10),
  text_grob("Health System Index", face = "bold", size = 10),
  text_grob("Vaccination Index", face = "bold", size = 10),
  ncol = 4)

fig_final <- ggarrange(
  header_row,
  fig1,
  nrow = 2,
  heights = c(0.06, 1))

fig_final












