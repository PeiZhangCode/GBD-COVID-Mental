######################################################################################################
# this file is for calculating excess mental disorder burden using ARIMA/ARIMA-LSTM/ARIMA + ARIMA-LSTM
# ARIMA + ARIMA-LSTM was used in main analysis, while ARIMA and ARIMA-LSTM was for sensitivity analysis
# Pei Zhang
# 2026-04-30
######################################################################################################
### packages
packages <- c('foreign','tsModel','lmtest','Epi','splines','vcd','scales','dplyr','ggplot2','ggsci',
              'RColorBrewer','ggpubr','tidyr','ggtext','lubridate','countrycode','sf','RColorBrewer',
              'patchwork','gtable','grid','broom','performance','forcats','readr','reticulate','forecast', 
              'tensorflow','keras','purrr','progressr')
lapply(packages, require, character.only=T) 

### python env
use_condaenv("r-tf", required = TRUE)
py_config()
py_module_available("tensorflow")
SEED <- 2026
Sys.setenv(TF_DETERMINISTIC_OPS = "1")
Sys.setenv(PYTHONHASHSEED = as.character(SEED))

### set seed
set_all_seeds <- function(seed) {
  seed <- as.integer(seed)
  set.seed(seed)
  
  reticulate::py_run_string(sprintf("
  import os, random
  os.environ['PYTHONHASHSEED'] = '%d'
  random.seed(%d)
  ", seed, seed))
  
  if (reticulate::py_module_available("numpy")) {
    np <- reticulate::import("numpy", delay_load = TRUE)
    np$random$seed(seed)}
  
  try(tensorflow::tf$random$set_seed(seed), silent = TRUE)
  try(keras::k_clear_session(), silent = TRUE)
  
  invisible(seed)
}

set_all_seeds(SEED)

###### 1. parameter selection ######
lstm_residual_forecast_mean <- function(arima_fit, arima_fc, h, lag, units,
                                        lr, epochs, batch_size, seed) {
  seed <- as.integer(seed)
  set_all_seeds(seed)
  
  lag <- as.integer(lag)
  units <- as.integer(units)
  epochs <- as.integer(epochs)
  batch_size <- as.integer(batch_size)
  h <- as.integer(h)
  
  res <- as.numeric(residuals(arima_fit))
  if (length(res) <= lag + 8L) return(as.numeric(arima_fc$mean))
  
  y <- res[(lag + 1L):length(res)]
  x <- sapply(seq_len(lag), function(k) {
    res[(lag + 1L - k):(length(res) - k)]
  })
  x <- as.matrix(x)
  
  x_scaled <- scale(x)
  x_center <- as.numeric(attr(x_scaled, "scaled:center"))
  x_scale <- as.numeric(attr(x_scaled, "scaled:scale"))
  
  if (any(is.na(x_scale)) || any(!is.finite(x_scale)) || any(x_scale == 0)) {
    return(as.numeric(arima_fc$mean))}
  
  y_mean <- mean(y)
  y_sd <- sd(y)
  if (is.na(y_sd) || y_sd == 0) return(as.numeric(arima_fc$mean))
  y_scaled <- (y - y_mean) / y_sd
  
  x_scaled <- array(as.numeric(x_scaled), dim = c(nrow(x_scaled), lag, 1L))
  y_scaled <- as.numeric(y_scaled)
  
  x_tf <- tensorflow::tf$convert_to_tensor(x_scaled, dtype = "float32")
  y_tf <- tensorflow::tf$convert_to_tensor(y_scaled, dtype = "float32")
  
  model <- keras_model_sequential()
  model$add(layer_input(shape = c(lag, 1L)))
  model$add(layer_lstm(units = units))
  model$add(layer_dense(units = 1L))
  
  model$compile(
    optimizer = optimizer_adam(learning_rate = lr),
    loss = "mse")
  
  model$fit(
    x_tf, y_tf,
    epochs = epochs,
    batch_size = batch_size,
    verbose = 0,
    shuffle = FALSE)
  
  last <- rev(tail(res, lag))
  last_x_scaled <- (last - x_center) / x_scale
  last_x_scaled <- as.numeric(last_x_scaled)
  
  pred_y_scaled <- numeric(h)
  
  for (i in seq_len(h)) {
    xin <- array(last_x_scaled, dim = c(1L, lag, 1L))
    xin_tf <- tensorflow::tf$convert_to_tensor(xin, dtype = "float32")
    p <- as.numeric(model$predict(xin_tf, verbose = 0))
    pred_y_scaled[i] <- p
    last_x_scaled <- c(last_x_scaled[-1], p)}
  
  pred_residual <- pred_y_scaled * y_sd + y_mean
  as.numeric(arima_fc$mean) + pred_residual
}

tune_one_setting <- function(data, measure1, location1, sex1, age1, cause1,
                             lag, units, lr, epochs, batch_size, seed) {
  global_data <- subset(
    data,
    measure == measure1 &
      location == location1 &
      sex == sex1 &
      age == age1 &
      cause == cause1)
  
  # train: <= 2015 ; validation: 2016-2019
  forecast_data <- subset(global_data, year <= 2015)
  validation_data <- subset(global_data, year %in% 2016:2019)
  validation_data <- validation_data[order(validation_data$year), ]
  
  if (nrow(forecast_data) == 0 || nrow(validation_data) != 4) return(NULL)
  if (all(is.na(forecast_data$val)) || all(forecast_data$val == 0)) return(NULL)
  
  y_train <- ts(
    forecast_data$val,
    start = min(forecast_data$year),
    frequency = 1)
  
  ar_fit <- auto.arima(y_train, stepwise = FALSE, approximation = FALSE)
  ar_fc <- forecast(ar_fit, h = 4, level = 95)
  
  pred_lstm <- lstm_residual_forecast_mean(
    arima_fit = ar_fit,
    arima_fc = ar_fc,
    h = 4,
    lag = lag,
    units = units,
    lr = lr,
    epochs = epochs,
    batch_size = batch_size,
    seed = seed)
  
  y_true <- as.numeric(validation_data$val)
  
  rmse_lstm <- sqrt(mean((y_true - pred_lstm)^2))
  mae_lstm <- mean(abs(y_true - pred_lstm))
  
  data.frame(
    measure = measure1,
    location = location1,
    sex = sex1,
    age = age1,
    cause = cause1,
    lag = lag,
    units = units,
    lr = lr,
    epochs = epochs,
    batch_size = batch_size,
    rmse_lstm = rmse_lstm,
    mae_lstm = mae_lstm,
    stringsAsFactors = FALSE)
}

tune_global_hyperparams <- function(data,
                                    lag_grid = c(1, 2, 3),
                                    lr_grid = c(0.01, 0.05),
                                    epochs_grid = c(60, 100, 200, 300),
                                    units = 32,
                                    batch_size = 8,
                                    seed) {
  data <- data[order(data$year), ]
  
  measure1 <- data$measure[1]
  location1 <- data$location[1]
  sex1 <- data$sex[1]
  age1 <- data$age[1]
  cause1 <- data$cause[1]
  
  grid <- expand.grid(
    lag = lag_grid,
    lr = lr_grid,
    epochs = epochs_grid,
    stringsAsFactors = FALSE)
  
  res_list <- vector("list", nrow(grid))
  
  for (i in seq_len(nrow(grid))) {
    cat(
      "Running", i, "/", nrow(grid),
      " lag =", grid$lag[i],
      " lr =", grid$lr[i],
      " epochs =", grid$epochs[i], "\n")
    
    res_list[[i]] <- tune_one_setting(
      data = data,
      measure1 = measure1,
      location1 = location1,
      sex1 = sex1,
      age1 = age1,
      cause1 = cause1,
      lag = grid$lag[i],
      units = units,
      lr = grid$lr[i],
      epochs = grid$epochs[i],
      batch_size = batch_size,
      seed = seed)
    gc()}
  
  res_df <- bind_rows(res_list)
  
  best_row <- res_df %>%
    arrange(rmse_lstm, mae_lstm) %>%
    slice(1)
  
  list(all_results = res_df, best_params = best_row)
}

dat <- read.csv("IHME-GBD_2023_DATA-a33adbeb-1.csv") %>%
  select(location_id, location_name, measure_name, sex_name, age_name, cause_name, year, val) %>%
  filter(cause_name == "Mental disorders")

names(dat) <- c("location_id", "location", "measure", "sex", "age", "cause", "year", "val")

dat1 <- dat %>%
  mutate(
    year = as.integer(year),
    val = as.numeric(val),
    location = as.character(location)) %>%
  arrange(year)

tune_out <- tune_global_hyperparams(
  data = dat1,
  lag_grid = c(1, 2, 3),
  lr_grid = c(0.01, 0.05),
  epochs_grid = c(60, 100, 200, 300),
  units = 32,
  batch_size = 8,
  seed = SEED)

tune_results <- tune_out$all_results
print(tune_results)

best_params <- tune_out$best_params
print(best_params)

# write.csv(best_params, "best_params.csv", row.names = FALSE)

###### 2. model selection & predict 2020-2023 ######
model_valid <- function(data, measure1, location1, sex1, age1, cause1,
                        lag, units, lr, epochs, batch_size, seed) {
  global_data <- subset(data,
    measure == measure1 & location == location1 & sex == sex1 & age == age1 & cause == cause1)
  
  # train: <=2015 ; validate: 2016-2019
  forecast_data <- subset(global_data, year <= 2015)
  validation_data <- subset(global_data, year %in% 2016:2019)
  validation_data <- validation_data[order(validation_data$year), ]
  
  if (nrow(forecast_data) == 0 || nrow(validation_data) != 4) return(NULL)
  if (all(is.na(forecast_data$val)) || all(forecast_data$val == 0)) return(NULL)
  
  y_train <- ts(forecast_data$val, start = min(forecast_data$year), frequency = 1)
  
  ar_fit <- auto.arima(y_train, stepwise = FALSE, approximation = FALSE)
  ar_fc <- forecast(ar_fit, h = 4, level = 95)
  pred_arima <- as.numeric(ar_fc$mean)
  
  pred_lstm <- lstm_residual_forecast_mean(
    ar_fit, ar_fc, h = 4,
    lag = lag, units = units, lr = lr,
    epochs = epochs, batch_size = batch_size,
    seed = seed)
  
  y_true <- as.numeric(validation_data$val)
  
  rmse_arima <- sqrt(mean((y_true - pred_arima)^2))
  mae_arima <- mean(abs(y_true - pred_arima))
  
  rmse_lstm <- sqrt(mean((y_true - pred_lstm)^2))
  mae_lstm <- mean(abs(y_true - pred_lstm))
  
  model_best <- ifelse(rmse_lstm < rmse_arima & mae_lstm < mae_arima, "ARIMA-LSTM", "ARIMA")
  
  data.frame(
    measure = measure1,
    location = location1,
    sex = sex1,
    age = age1,
    cause = cause1,
    rmse_arima = rmse_arima,
    mae_arima = mae_arima,
    rmse_lstm = rmse_lstm,
    mae_lstm = mae_lstm,
    model_best = model_best,
    stringsAsFactors = FALSE)
}

arima_model <- function(data) {
  data1 <- subset(data, year <= 2019)
  data1 <- data1[order(data1$year), ]
  dt1 <- ts(data1$val, start = min(data1$year), frequency = 1)
  
  fit <- auto.arima(dt1, stepwise = FALSE, approximation = FALSE)
  fc <- forecast(fit, h = 4, level = 95)
  
  data.frame(
    model = "ARIMA",
    year = 2020:2023,
    mean = as.numeric(fc$mean),
    pre_lower = as.numeric(fc$lower[, 1]),
    pre_upper = as.numeric(fc$upper[, 1]),
    stringsAsFactors = FALSE)
}

arima_lstm_model <- function(data, lag, units, lr, epochs, batch_size, seed) {
  data1 <- subset(data, year <= 2019)
  data1 <- data1[order(data1$year), ]
  dt1 <- ts(data1$val, start = min(data1$year), frequency = 1)
  
  fit <- auto.arima(dt1, stepwise = FALSE, approximation = FALSE)
  fc <- forecast(fit, h = 4, level = 95)
  
  mean_hybrid <- lstm_residual_forecast_mean(
    fit, fc, h = 4,
    lag = lag, units = units, lr = lr,
    epochs = epochs, batch_size = batch_size,
    seed = seed)
  
  data.frame(
    model = "ARIMA-LSTM",
    year = 2020:2023,
    mean = as.numeric(mean_hybrid),
    pre_lower = as.numeric(fc$lower[, 1]),
    pre_upper = as.numeric(fc$upper[, 1]),
    stringsAsFactors = FALSE)
}

process_model_group <- function(data,
                                model_strategy = c("select", "arima", "arima_lstm"),
                                lag, units, lr, epochs, batch_size, seed) {
  model_strategy <- match.arg(model_strategy)
  
  measure_val <- data$measure[1]
  sex_val <- data$sex[1]
  cause_val <- data$cause[1]
  age_val <- data$age[1]
  location_val <- data$location[1]
  
  if (model_strategy == "select") {
    model_result <- model_valid(
      data = data,
      measure1 = measure_val,
      location1 = location_val,
      sex1 = sex_val,
      age1 = age_val,
      cause1 = cause_val,
      lag = lag,
      units = units,
      lr = lr,
      epochs = epochs,
      batch_size = batch_size,
      seed = seed)
    
    if (is.null(model_result)) {
      return(list(model_validation = NULL, prediction = NULL))}
    
    chosen_model <- model_result$model_best[1]
    
  } else if (model_strategy == "arima") {
    
    model_result <- data.frame(
      measure = measure_val,
      location = location_val,
      sex = sex_val,
      age = age_val,
      cause = cause_val,
      rmse_arima = NA_real_,
      mae_arima = NA_real_,
      rmse_lstm = NA_real_,
      mae_lstm = NA_real_,
      model_best = "ARIMA",
      stringsAsFactors = FALSE)
    
    chosen_model <- "ARIMA"
    
  } else if (model_strategy == "arima_lstm") {
    
    model_result <- data.frame(
      measure = measure_val,
      location = location_val,
      sex = sex_val,
      age = age_val,
      cause = cause_val,
      rmse_arima = NA_real_,
      mae_arima = NA_real_,
      rmse_lstm = NA_real_,
      mae_lstm = NA_real_,
      model_best = "ARIMA-LSTM",
      stringsAsFactors = FALSE)
    
    chosen_model <- "ARIMA-LSTM"
    }
  
  if (chosen_model == "ARIMA") {
    pred_result <- arima_model(data)
  } else {
    pred_result <- arima_lstm_model(
      data,
      lag = lag,
      units = units,
      lr = lr,
      epochs = epochs,
      batch_size = batch_size,
      seed = seed)
  }
  
  pred_result$measure <- measure_val
  pred_result$location <- location_val
  pred_result$sex <- sex_val
  pred_result$cause <- cause_val
  pred_result$age <- age_val
  
  list(model_validation = model_result, prediction = pred_result)
}

dat <- read.csv("IHME-GBD_2023_DATA-2b3e5b00-1.csv")
dat1 <- dat %>%
  select(location_id, location_name, measure_name, sex_name, age_name, cause_name, year, val)
names(dat1) <- c("location_id", "location_name", "measure", "sex", "age", "cause", "year", "val")
dat1$location <- dat1$location_name

dat1 <- dat1 %>%
  mutate(
    year = as.integer(year),
    val = as.numeric(val),
    location = as.character(location))

out_dir <- "output/result_clean"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

model_strategy <- "select" # "select" / "arima" / "arima_lstm"

groups <- dat1 %>%
  group_by(measure, sex, cause, location, age) %>%
  group_split()

cat("Total groups:", length(groups), "\n")

chunk_size <- 50
chunks <- split(groups, ceiling(seq_along(groups) / chunk_size))

all_val_list <- list()
all_pred_list <- list()

handlers(global = TRUE)
with_progress({
  p <- progressor(along = seq_along(chunks))
  
  for (batch_id in seq_along(chunks)) {
    message(
      "Running batch ", batch_id, " / ", length(chunks),
      " (groups: ", length(chunks[[batch_id]]), ")")
    
    res_list <- lapply(chunks[[batch_id]], function(df_one_group) {
      out <- process_model_group(
        df_one_group,
        model_strategy = model_strategy,
        lag = 2,
        units = 32,
        lr = 0.01,
        epochs = 100,
        batch_size = 8,
        seed = SEED)
      
      if (is.data.frame(out) && all(c("model_validation", "prediction") %in% names(out))) {
        val <- out$model_validation[[1]]
        pred <- out$prediction[[1]]
      } else if (is.list(out) && all(c("model_validation", "prediction") %in% names(out))) {
        val <- out$model_validation
        pred <- out$prediction
      } else {
        val <- NULL
        pred <- NULL}
      list(val = val, pred = pred)
    })
    
    batch_val <- bind_rows(lapply(res_list, `[[`, "val"))
    batch_pred <- bind_rows(lapply(res_list, `[[`, "pred"))
    
    all_val_list[[batch_id]] <- batch_val
    all_pred_list[[batch_id]] <- batch_pred
    
    rm(res_list, batch_val, batch_pred)
    gc()
    p()
  }
})

best_model_df <- bind_rows(all_val_list)
prediction_df <- bind_rows(all_pred_list)

# write.csv(best_model_df, file = file.path(out_dir, "best_model.csv"), row.names = FALSE)

###### 3. calculate excess for countries #######
pred <- prediction_df

obs <- dat1 %>%
  mutate(
    year = as.integer(year),
    val = as.numeric(val)
  ) %>%
  filter(year %in% 2020:2023) %>%
  select(measure, sex, cause, age, location, year, observed = val)

pred2 <- pred %>%
  mutate(
    year = as.integer(year),
    mean = as.numeric(mean)
  ) %>%
  select(
    measure, sex, cause, age, location, year,
    model, model_best = model, predicted = mean,
    pre_lower, pre_upper)

pred_excess <- pred2 %>%
  left_join(obs, by = c("measure", "sex", "cause", "age", "location", "year")) %>%
  mutate(
    excess = observed - predicted,
    excess_lower = observed - as.numeric(pre_upper),
    excess_upper = observed - as.numeric(pre_lower))

cat("Rows with missing observed:", sum(is.na(pred_excess$observed)), "\n")

# write.csv(
#   pred_excess,
#   file = file.path(out_dir, paste0("F_country_prediction_with_excess_", model_strategy, "_noshuffle_seed", SEED, ".csv")),
#   row.names = FALSE
# )

### excess mental health burden due to covid-19 during 2020-2023 (consider population)
mental_post <- pred_excess %>%
  select(location, year, predicted,
         extra_dalys=excess, 
         extra_low=excess_lower, 
         extra_high=excess_upper)

pop <- read.csv("IHME-GBD_2023_DATA-3f959b33-1.csv") %>% 
  select(location = location_name, year, popn=val) %>%
  mutate(popn = popn/100000)

mental_post2 <- mental_post %>%
  left_join(pop, by = c("location", "year"))

summary_dalys <- mental_post2 %>%
  filter(year %in% 2020:2023) %>%
  group_by(location) %>%
  summarise(
    avg_extra_rate = round(sum(extra_dalys * popn, na.rm = TRUE) / sum(popn, na.rm = TRUE), 1),
    avg_extral_rate = round(sum(extra_low   * popn, na.rm = TRUE) / sum(popn, na.rm = TRUE), 1),
    avg_extrah_rate = round(sum(extra_high  * popn, na.rm = TRUE) / sum(popn, na.rm = TRUE), 1),
    .groups = "drop")

summary_dalys <- summary_dalys %>%
  mutate(
    dalys_summary = ifelse(
      is.na(avg_extra_rate) | is.na(avg_extral_rate) | is.na(avg_extrah_rate),
      NA,
      sprintf("%.1f\n(%.1f, %.1f)", avg_extra_rate, avg_extrah_rate, avg_extral_rate))
  ) %>%
  arrange(desc(avg_extra_rate)) %>%
  select(location, dalys_summary)

# write.csv(summary_dalys, file = "tableS1_country_excess_burden_F_seed.csv")

###### 4. calculate excess for global #######
dat <- read.csv("IHME-GBD_2023_DATA-a33adbeb-1.csv") %>%
  select(location_id, location_name, measure_name, sex_name, age_name, cause_name, year, val) %>%
  filter(cause_name == "Mental disorders")

names(dat) <- c("location_id", "location", "measure", "sex", "age", "cause", "year", "val")

dat1 <- dat %>%
  mutate(
    year = as.integer(year),
    val = as.numeric(val),
    location = as.character(location)
  ) %>%
  arrange(year)

model_strategy <- "select"  # "arima" / "arima_lstm" / "select"

out <- process_model_group(
  dat1,
  model_strategy = model_strategy,
  lag = 2,
  units = 32,
  lr = 0.01,
  epochs = 100,
  batch_size = 8,
  seed = SEED)

best_model_df <- out$model_validation
prediction_df <- out$prediction

obs <- dat1 %>%
  filter(year %in% 2020:2023) %>%
  select(measure, sex, cause, age, location, year, observed = val)

pred2 <- prediction_df %>%
  mutate(
    year = as.integer(year),
    predicted = as.numeric(mean)
  ) %>%
  select(measure, sex, cause, age, location, year, model, predicted, pre_lower, pre_upper)

pred_excess <- pred2 %>%
  left_join(obs, by = c("measure", "sex", "cause", "age", "location", "year")) %>%
  mutate(
    excess = observed - predicted,
    excess_lower = observed - as.numeric(pre_upper),
    excess_upper = observed - as.numeric(pre_lower))

pred_excess

### global excess mental health burden due to covid-19 during 2020-2023 (consider population change)
mental_post <- pred_excess %>%
  select(year, model, predicted,
         extra_dalys=excess, 
         extra_low=excess_lower, 
         extra_high=excess_upper)

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


