# Clear Workspace
#remove(list=ls())

# Load Libraries ----------------------------------------------------------
library(aTSA); library(astsa); library(car); library(chron); library(clock)
library(cowplot); library(datasets); library(DEoptim); library(dplyr)
library(fable); library(fBasics); library(fImport); library(forecast)
library(ForecastComb); library(fpp); library(fpp2); library(fpp3)
library(gam); library(GenSA); library(ggplot2); library(GGally)
library(glmnet); library(gridExtra); library(highfrequency); library(kableExtra)
library(knitr); library(lavaan); library(lme4); library(lubridate)
library(MASS); library(mice); library(Metrics); library(mvtnorm)
library(plotly); library(PortfolioAnalytics); library(quantmod)
library(randomForest); library(readxl); library(reshape); library(reshape2)
library(robust); library(ROI); library(ROI.plugin.glpk); library(ROI.plugin.quadprog)
library(sandwich); library(sn); library(stats); library(timeSeries)
library(tibble); library(tidyr); library(tidyverse); library(tseries)
library(tsibble); library(xtable); library(xts); library(zoo)

# Load External Scripts
# 
#  ---------------------------------------------------
source("/.../sse_v3.lib")
source("/.../for.perf.R")

# Load and Prepare Data ---------------------------------------------------
df <- read_xlsx("/.../dati_pun.xlsx", sheet = "Sheet1")
df$Date <- as.Date(df$Date, format = "%Y%m%d")
df$Hour <- as.numeric(df$Hour)

# Split the data frame by Hour into lists of data frames
df_list <- split(df, df$Hour)
df1 <- subset(df, Date >= as.Date("2021-01-01") & Date <= as.Date("2023-01-31"))
df1_list <- split(df1, df1$Hour)
df2 <- subset(df, Date >= as.Date("2023-02-01") & Date <= as.Date("2023-04-30"))
df2_list <- split(df2, df2$Hour)

# Define colors for visualizations
my_colors <- c("#4C72B0", "#DD8452")

# Functions ----------------------------------------------------------------
# Rolling forecast function
rolling_forecast <- function(hour, p, d, q, P, D, Q) {
  forecast_values <- c()
  for (i in 0:88) {
    model <- Arima(df_list[[hour]]$PUN[(1+i):(length(df1_list[[hour]]$PUN) + i)], 
                   order = c(p, d, q), seasonal = list(order = c(P, D, Q), period = 7), include.drift = FALSE)
    forecast_values[i+1] <- forecast(model)$mean[1]
  }
  return(forecast_values)
}

# SMAPE calculation function
smape <- function(a, f) {
  return (1 / length(a) * sum(2 * abs(f - a) / (abs(a) + abs(f)) * 100))
}

# RMSE calculation function
rmse <- function(actual, forecast) {
  sqrt(mean((actual - forecast)^2))
}

# Combined error metrics function
error_metrics <- function(actual, forecast) {
  rmse_val <- rmse(actual, forecast)
  mase_val <- mase(actual, forecast)
  smape_val <- smape(actual, forecast)
  
  metrics_df <- data.frame(
    RMSE = rmse_val,
    MASE = mase_val,
    SMAPE = smape_val
  )
  
  return(metrics_df)
}

# Lagged features creation function
create_lagged_features <- function(data, num_lags) {
  for (lag in 1:num_lags) {
    lag_name <- paste0("PUN_lag_", lag)
    data <- data %>% mutate(!!lag_name := lag_func(PUN, lag))
  }
  return(data)
}

# Stationarity Tests -------------------------------------------------------
# Data frame for Dickey-Fuller test results
results <- data.frame(Series = character(), p_value = numeric(), p_value_diff = numeric(), stringsAsFactors = FALSE)
for (i in 1:length(df1_list)) {
  series <- log(df1_list[[i]]$PUN)
  result <- adf.test(series)
  series_diff <- diff(series)
  result_diff <- adf.test(series_diff)
  results <- rbind(results, data.frame(Series = paste0("Time Series ", i), p_value = result$p.value, p_value_diff = result_diff$p.value, stringsAsFactors = FALSE))
}
latex_table <- xtable(results, caption = "Dickey-Fuller Test Results", label = "table:df-test")
print(latex_table, include.rownames = FALSE)

# Data frame for combined Phillips-Perron test results
results_combined <- data.frame(Series = character(), p_value_original = numeric(), p_value_differenced = numeric(), stringsAsFactors = FALSE)
for (i in 1:length(df1_list)) {
  series <- log((df1_list[[i]]$PUN))
  result_original <- pp.test(series)
  series_diff <- diff(series)
  result_differenced <- pp.test(series_diff)
  results_combined <- rbind(results_combined, data.frame(Series = paste0("Time Series ", i), p_value_original = result_original$p.value, p_value_differenced = result_differenced$p.value, stringsAsFactors = FALSE))
}
latex_table_combined <- xtable(results_combined, caption = "Phillips-Perron Test Results", label = "table:pp-test-combined")
print(latex_table_combined, include.rownames = FALSE)


# Plots and Visualization --------------------------------------------------

# Create plots for each hour with data from df1_list
plot_list <- list()
for (i in seq_along(df1_list)) {
  p <- ggplot(df1_list[[i]], aes(x = Date, y = PUN)) +
    geom_line(aes(color = factor(Hour))) +
    geom_ribbon(aes(ymin = 0, ymax = PUN, fill = factor(Hour)), alpha = 0.2) +
    geom_smooth(aes(group = i), color = "black", se = FALSE) +
    labs(x = "Date", y = "PUN", title = paste0("Hourly PUN data for hour ", i),
         subtitle = "January 1, 2021 - January 31, 2023") +
    theme(plot.title = element_text(size=10, hjust = 0.5),
          plot.subtitle = element_text(size=6, hjust = 0.5),
          legend.position = "none",
          axis.text.x = element_blank()) +
    scale_x_date(date_labels = "%m/%d/%y", limits = c(as.Date("2021-01-01"), as.Date("2023-01-31"))) +
    scale_color_manual(values = my_colors) +
    ylim(0, 1000)
  
  plot_list[[i]] <- p
}

# Display and arrange multiple plots using grid.arrange
for (j in seq(1, length(plot_list), by = 6)) {
  grid.arrange(grobs = plot_list[j:min(j+5, length(plot_list))], ncol = 3)
}

# Boxplot of PUN values ----------------------------------------------------

# Convert PUN to numeric and create box plots
df1_list <- lapply(df1_list, function(x) { x$PUN <- as.numeric(x$PUN); x })

plot_list <- list()
for (i in seq_along(df1_list)) {
  p <- ggplot(df1_list[[i]], aes(x = "", y = PUN)) +
    geom_boxplot(fill = my_colors[i %% 2 + 1], color = "black", alpha = 0.8) +
    geom_hline(yintercept = quantile(df1_list[[i]]$PUN, probs = c(0.25, 0.75)) + c(-1.5, 1.5) * IQR(df1_list[[i]]$PUN), linetype = "dashed", color = "red") +
    labs(title = paste0("Hour ", i)) +
    theme_classic() +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          plot.margin = unit(c(0.5, 0.5, 1, 1), "lines")) +
    coord_cartesian(ylim = c(0, 1000))
  
  plot_list[[i]] <- p
}

# Combine and display box plots
combined_plot <- do.call(gridExtra::grid.arrange, c(plot_list, ncol = 12))
print(combined_plot)

# ARIMA Modeling and Error Metrics -----------------------------------------

library(forecast)

# Create forecasts using ARIMA models
window_size <- length(df_list[[1]]$PUN) - 89
ForecastArima <- data.frame(matrix(NA, nrow = 89, ncol = 24))
colnames(ForecastArima) <- paste0("ForecastArima", 1:24)

for (j in 1:24) {
  model_order <- if (j >= 2 & j <= 19) c(7, 0, 0) else c(7, 1, 0)
  data <- log(df_list[[j]]$PUN)
  
  for (i in 1:89) {
    model <- Arima(data[i:(window_size+i-1)], order = model_order, include.drift = FALSE)
    ForecastArima[i, j] <- exp(forecast(model)$mean[1])
  }
  print(j)
}

# Error Metrics for ARIMA
error_metrics_ARIMA <- data.frame(Series = integer(), Model = character(), RMSE = numeric(), MASE = numeric(), SMAPE = numeric())

model_arima <- c("ARIMA(7,1,0)", rep("ARIMA(7,0,0)", 18), "ARIMA(7,1,0)")

for (hour in 1:24) {
  error_metrics_1 <- error_metrics(df2_list[[hour]]$PUN, ForecastArima[,hour])
  error_metrics_ARIMA <- rbind(error_metrics_ARIMA, cbind(Series = hour, Model = model_arima[hour], error_metrics_1))
}

# Display and format error metrics using kable
kable(error_metrics_ARIMA, caption = "Error Metrics for ARIMA", format = "latex", booktabs = TRUE) %>%
  kable_styling(font_size = 10)

window_size <- length(df_list[[20]]$PUN)-89


# Create an empty dataframe for the comparison results
comparison <- data.frame(
  Series = integer(), 
  Lowest_RMSE_Model = character(), 
  Lowest_MASE_Model = character(), 
  Lowest_SMAPE_Model = character(), 
  stringsAsFactors = FALSE
)

# Define a function that will compare models for a given metric
compare_models <- function(series, metric) {
  models <- list(error_metrics_hw, error_metrics_Holt, error_metrics_ETS)
  model_names <- c("hw", "Holt", "ETS")
  model_values <- sapply(models, function(model) model[model$Series == series, metric])
  names(model_values) <- model_names
  return(names(which.min(model_values)))
}

#found based on aic criteria
ets_models <- list("MNN", "MNN", "MMN", "MMN", "MNN", "MNN", "MAM", "MAM", "MNM", "MNM", "MNM", "MAM", "MMM", "MMM", "MAM", "MAM", "MNM", "MAA", "MAM", "MNA", "MNN", "MNN", "MMM", "MNM")
column_names <- paste0("ForecastETS", 1:24)
ForecastETS<- data.frame(matrix(NA, nrow = 0, ncol = length(column_names)))
colnames(ForecastETS) <- column_names

for (hour in 1:24) {
  timeseries <- df_list[[hour]]$PUN
  forecasts_EXPS <- c()
  for (i in 1:89) {
    train_data <- ts(timeseries[(length(timeseries) - window_size - 88 + i - 1):(length(timeseries) - 89 + i - 1)], frequency = 7)
    model <- forecast::ets(train_data, model = ets_models[hour])
    ForecastETS[i,hour] <- forecast::forecast(model, h = 1)$mean
  }
  print(hour)
  print(ets_models[hour])}


window_size <- length(df_list[[20]]$PUN)-89


#Tabelle LateX AppendiX
library(kableExtra)


######################################################################
####################     SPLINE ######################################
######################################################################

window_size <- length(df_list[[1]]$PUN) - 89
column_names <- paste0("ForecastSPLINE", 1:24)
ForecastSPLINE <- data.frame(matrix(NA, nrow = 89, ncol = length(column_names)))
colnames(ForecastSPLINE) <- column_names

# For each time series...
for (j in 1:24) {
  # Plot the series
  plot(df_list[[j]]$PUN)
  # Fit the spline
  spline_fit <- smooth.spline(time(df_list[[j]]$Date), df_list[[j]]$PUN, spar = 1)
  # Generate a sequence of points within the range of the original data
  time_points <- seq(min(time(df_list[[j]]$Date)), max(time(df_list[[j]]$Date)))
  # Evaluate the spline at the exact same time points as the original series
  spline_points <- predict(spline_fit, time(df_list[[j]]$Date))
  # Subtract the spline values from the original series
  test1spline <- df_list[[j]]$PUN - spline_points$y
  # Initialize an empty vector to hold the forecasts
  for1spline <- vector()
  # For each window in the series...
  for (i in 1:(length(test1spline) - window_size)) {
    # Fit an ARIMA(7,0,0) model
    model <- Arima(test1spline[i:(i + window_size - 1)], order = c(7,0,0), include.drift = FALSE)
    
    # Forecast the next value
    for1spline <- c(for1spline, forecast(model)$mean[1])
  }
  
  # Add the spline values back to the forecasts
  prova <- spline_points$y[(window_size + 1):length(spline_points$y)] + for1spline
  
  # Store the forecasts in the dataframe
  ForecastSPLINE[,j] <- prova
  print(j)
}


###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
###############################################################################
# Determine the number of forecasts for each time series
num_forecasts <- 89
# Initialize a matrix to store the forecasts for each time series
forecasts_matrix <- matrix(nrow = num_forecasts, ncol = 24)

# For each time series
for (j in 1:24) {
  # Extract the time series data
  data <- data.frame(PUN = df_list[[j]]$PUN, Date = df_list[[j]]$Date)
  # Create the lagged features
  data <- create_lagged_features(data, 7)
  # Remove rows with missing values (these will be the first few rows due to the creation of lagged features)
  data <- na.omit(data)
  # Initial training data from January 1st, 2021 to January 31st, 2023
  start_date <- as.Date("2021-01-01")
  end_date <- as.Date("2023-01-31")
  # Initialize an empty vector to hold the forecasts for this time series
  forecasts <- vector("numeric", length = nrow(data[data$Date >= as.Date("2023-02-01"), ]))
  # For each day in the test set, retrain the model on all data up to that day, and forecast the next day
  for (i in 1:length(forecasts)) {
    # All data up to the current day
    train <- data[data$Date >= start_date & data$Date <= end_date, ]
    # Retrain the model
    model <- randomForest(PUN ~ ., data = train, ntree = 1000)
    # Make a one-step ahead forecast
    if (i <= length(forecasts)) {
      forecasts[i] <- predict(model, newdata = data[data$Date == end_date + 1, ])
    }
    # Shift the training window
    start_date <- start_date + 1
    end_date <- end_date + 1
  }
  # Store the forecasts for this time series in the matrix
  forecasts_matrix[, j] <- forecasts
  print(j)
}
forecastRandomForest=forecasts_matrix


ForecastCombination <- data.frame(matrix(NA, nrow = nrow(ForecastArima), ncol = 24))


result_table <- matrix(NA, nrow = 24, ncol = 6)
colnames(result_table) <- c("Column Index", "ARIMA", "ETS","SPLINE","Random Forest", "Combination")

for (column_index in 1:24) {
  result_arima <- for.perf(df_list[[column_index]]$PUN, ForecastArima[, column_index])$U.theil
  result_ets <- for.perf(df_list[[column_index]]$PUN, ForecastETS[, column_index])$U.theil
  result_spline <- for.perf(df_list[[column_index]]$PUN, ForecastSPLINE[, column_index])$U.theil
  result_RF =for.perf(df_list[[column_index]]$PUN,forecastRandomForest[,column_index])$U.theil
  result_combination <- for.perf(df_list[[column_index]]$PUN, ForecastCombination[, column_index])$U.theil
  
  result_table[column_index, 1] <- column_index
  result_table[column_index, 2] <- round(result_arima, 2)
  result_table[column_index, 3] <- round(result_ets, 2)
  result_table[column_index,4] = round(result_spline, 2)
  result_table[column_index, 5] <- round(result_RF, 2)
  result_table[column_index, 6] <- round(result_combination, 2)
}

kable(result_table, format = "latex", booktabs = TRUE, caption = "RMSPE Results")
colMeans(result_table)


result_table <- matrix(NA, nrow = 24, ncol = 6)
colnames(result_table) <- c("Column Index", "ARIMA", "ETS","SPLINE","Random Forest", "Combination")

for (column_index in 1:24) {
  result_arima <- for.perf(df_list[[column_index]]$PUN, ForecastArima[, column_index])$mase
  result_ets <- for.perf(df_list[[column_index]]$PUN, ForecastETS[, column_index])$mase
  result_spline <- for.perf(df_list[[column_index]]$PUN, ForecastSPLINE[, column_index])$mase
  result_RF =for.perf(df_list[[column_index]]$PUN,forecastRandomForest[,column_index])$mase
  result_combination <- for.perf(df_list[[column_index]]$PUN, ForecastCombination[, column_index])$mase
  
  result_table[column_index, 1] <- column_index
  result_table[column_index, 2] <- round(result_arima, 2)
  result_table[column_index, 3] <- round(result_ets, 2)
  result_table[column_index,4] = round(result_spline, 2)
  result_table[column_index, 5] <- round(result_RF, 2)
  result_table[column_index, 6] <- round(result_combination, 2)
}

kable(result_table, format = "latex", booktabs = TRUE, caption = "RMSPE Results")
colMeans(result_table)


library(ggplot2)
library(purrr)
library(tidyr)
# Assuming `hours` is a vector containing the hour values (1 to 24)
hours <- 1:24

# Create a list of ggplot objects
plots <- map(hours, ~{
  
  hour <- .
  
  data_to_plot <- data.frame(
    Date = df2_list[[hour]]$Date,
    Actual = df2_list[[hour]]$PUN,
    ForecastArima = ForecastArima[, hour],
    ForecastETS = ForecastETS[, hour],
    ForecastSPLINE = ForecastSPLINE[,hour],
    ForecastCombination = ForecastCombination[,hour]
  )
  
  # Reshape data to long format
  data_to_plot_long <- data_to_plot %>%
    pivot_longer(cols = -Date, names_to = "Method", values_to = "PUN")
  
  ggplot(data_to_plot_long, aes(x = Date, y = PUN, color = Method)) +
    geom_line(aes(linetype = Method)) +
    labs(title = paste("Hour", hour)) +
    theme_minimal() +
    scale_color_manual(values = c("black", "blue", "red", "green", "purple")) +
    scale_linetype_manual(values = c("solid", "dashed", "dotted", "dotdash", "longdash"))
})



# Display the plot for the first hour
print(plots[[1]])

# Display the plot for the second hour
print(plots[[24]])


##########################################################################
##########################################################################
##########################################################################
##########################################################################
##########################################################################
# Create lagged features
lag_func <- function(x, n = 1) c(rep(NA, n), x)[1:length(x)]

create_lagged_features <- function(data, num_lags) {
  for (lag in 1:num_lags) {
    lag_name <- paste0("PUN_lag_", lag)
    data <- mutate(data, !!lag_name := lag_func(data$PUN, lag))
  }
  return(data)
}

data <- data.frame(PUN = df_list[[21]]$PUN,Date=df_list[[21]]$Date)
data <- create_lagged_features(data, 7)

# Remove rows with missing values (these will be the first few rows due to the creation of lagged features)
data <- na.omit(data)

# Initial training data from January 1st, 2021 to January 31st, 2023
start_date <- as.Date("2021-01-01")
end_date <- as.Date("2023-01-31")

# Initialize an empty vector to hold the forecasts
forecasts <- vector("numeric", length = nrow(data[data$Date >= as.Date("2023-02-01"), ]))

# For each day in the test set, retrain the model on all data up to that day, and forecast the next day
for (i in 1:length(forecasts)) {
  # All data up to the current day
  train <- data[data$Date >= start_date & data$Date <= end_date, ]
  # Retrain the model
  model <- randomForest(PUN ~ ., data = train, ntree = 1000, mtry=7)
  # Make a one-step ahead forecast
  if (i <= length(forecasts)) {
    forecasts[i] <- predict(model, newdata = data[data$Date == end_date + 1, ])
  }
  # Shift the training window
  start_date <- start_date + 1
  end_date <- end_date + 1
}
forecast_RF <- forecasts
for.perf(df_list[[21]]$PUN, forecast_RF)  # Please provide the definition or source of the for.perf function

##############################################################################
##############################################################################
##############################################################################
##############################################################################





# ARIMA 
error_measures_ARIMA <- data.frame()
model_arima <- character(24)

# First 7 elements: "ARIMA(7,1,0)"
model_arima[1] <- "ARIMA(7,1,0)"

# Elements 8 to 18: "ARIMA(7,0,0)"
model_arima[2:19] <- "ARIMA(7,0,0)"

# Elements 19 to 24: "ARIMA(7,1,0)"
model_arima[20:24] <- "ARIMA(7,1,0)"
# For each model
error_measures_ARIMA <- data.frame()
for (j in 1:24) {
  # Get the actual and predicted values
  actual <- df_list[[j]]$PUN
  predicted <- ForecastArima[,j]
  
  # Calculate error metrics
  MSE <- for.perf(actual, predicted)$MSE
  RMSPE <- for.perf(actual, predicted)$RMSPE
  MAPE <- for.perf(actual, predicted)$MAPE
  MASE <- for.perf(actual, predicted)$mase
  
  # Store the error metrics in the dataframe
  error_measures_ARIMA <- rbind(error_measures_ARIMA, data.frame(
    Series = j,
    Model = model_arima[j],
    MSE = MSE,
    RMSPE = RMSPE,
    MAPE = MAPE,
    MASE = MASE
  ))
}

# Print the error metrics in a LaTeX table format
library(knitr)
library(kableExtra)
kable(error_measures_ARIMA, caption = "Error Metrics for ARIMA", format = "latex", booktabs = TRUE) %>%
  kable_styling(font_size = 10)
########################################################################
########################################################################
########################################################################
########################################################################
########################################################################

error_measures_RandomForest <- data.frame()
for (j in 1:24) {
  # Get the actual and predicted values
  actual <- df_list[[j]]$PUN
  predicted <- ForecastRandomForest[,j]
  
  # Calculate error metrics
  MSE <- for.perf(actual, predicted)$MSE
  RMSPE <- for.perf(actual, predicted)$RMSPE
  MAPE <- for.perf(actual, predicted)$MAPE
  MASE <- for.perf(actual, predicted)$mase
  
  # Store the error metrics in the dataframe
  error_measures_RandomForest <- rbind(error_measures_RandomForest, data.frame(
    Series = j,
    Model = "Random Forest",
    MSE = MSE,
    RMSPE = RMSPE,
    MAPE = MAPE,
    MASE = MASE
  ))
}

# Print the error metrics in a LaTeX table format
library(knitr)
library(kableExtra)
kable(error_measures_RandomForest, caption = "Error Metrics for Random Forest", format = "latex", booktabs = TRUE) %>%
  kable_styling(font_size = 10)


########################################################################
########################################################################
########################################################################
########################################################################
########################################################################

error_measures_ETS <- data.frame()
for (j in 1:24) {
  # Get the actual and predicted values
  actual <- df_list[[j]]$PUN
  predicted <- ForecastETS[,j]
  
  # Calculate error metrics
  MSE <- for.perf(actual, predicted)$MSE
  RMSPE <- for.perf(actual, predicted)$RMSPE
  MAPE <- for.perf(actual, predicted)$MAPE
  MASE <- for.perf(actual, predicted)$mase
  
  # Store the error metrics in the dataframe
  error_measures_ETS <- rbind(error_measures_ETS, data.frame(
    Series = j,
    Model = paste0(ets_models[j],j),
    MSE = MSE,
    RMSPE = RMSPE,
    MAPE = MAPE,
    MASE = MASE
  ))
}

# Print the error metrics in a LaTeX table format
library(knitr)
library(kableExtra)
kable(error_measures_ETS, caption = "Error Metrics for ETS", format = "latex", booktabs = TRUE) %>%
  kable_styling(font_size = 10)


########################################################################
########################################################################
########################################################################
########################################################################
########################################################################

error_measures_SPLINE <- data.frame()
for (j in 1:24) {
  # Get the actual and predicted values
  actual <- df_list[[j]]$PUN
  predicted <- ForecastSPLINE[,j]
  
  # Calculate error metrics
  MSE <- for.perf(actual, predicted)$MSE
  RMSPE <- for.perf(actual, predicted)$RMSPE
  MAPE <- for.perf(actual, predicted)$MAPE
  MASE <- for.perf(actual, predicted)$mase
  
  # Store the error metrics in the dataframe
  error_measures_SPLINE <- rbind(error_measures_SPLINE, data.frame(
    Series = j,
    Model = "SPLINE - ARIMA (7,0,0)",
    MSE = MSE,
    RMSPE = RMSPE,
    MAPE = MAPE,
    MASE = MASE
  ))
}

# Print the error metrics in a LaTeX table format
library(knitr)
library(kableExtra)
kable(error_measures_SPLINE, caption = "Error Metrics for SPLINE", format = "latex", booktabs = TRUE) %>%
  kable_styling(font_size = 10)



#save(ForecastArima,ForecastRandomForest,ForecastETS,ForecastSPLINE,file="Forecasts_thesis")


ForecastCombination_CSR <- matrix(nrow = nrow(df2_list[[1]]), ncol = 24)
ForecastCombination_LAD <- matrix(nrow = nrow(df2_list[[1]]), ncol = 24)
ForecastCombination_NG <- matrix(nrow = nrow(df2_list[[1]]), ncol = 24)
ForecastCombination_OLS <- matrix(nrow = nrow(df2_list[[1]]), ncol = 24)
ForecastCombination_SA <- matrix(nrow = nrow(df2_list[[1]]), ncol = 24)

for (hour in 1:24) {
  forecasts <- cbind(ForecastArima[, hour], ForecastETS[, hour], ForecastRandomForest[, hour], ForecastSPLINE[, hour])
  x <- foreccomb(df2_list[[hour]]$PUN, forecasts)
  
  # Combination using comb_CSR()
  a_CSR <- comb_CSR(x)
  ForecastCombination_CSR[, hour] <- a_CSR$Fitted[, 1]
  
  # Combination using comb_LAD()
  a_LAD <- comb_LAD(x)
  ForecastCombination_LAD[, hour] <- a_LAD$Fitted
  
  # Combination using comb_NG()
  a_NG <- comb_NG(x)
  ForecastCombination_NG[, hour] <- a_NG$Fitted
  
  # Combination using comb_OLS()
  a_OLS <- comb_OLS(x)
  ForecastCombination_OLS[, hour] <- a_OLS$Fitted
  
  # Combination using comb_SA()
  a_SA <- comb_SA(x)
  ForecastCombination_SA[, hour] <- a_SA$Fitted
}

error_measures_CSR <- data.frame()

for (j in 1:24) {
  actual <- df_list[[j]]$PUN
  predicted <- ForecastCombination_CSR[, j]
  
  MSE <- for.perf(actual, predicted)$MSE
  RMSPE <- for.perf(actual, predicted)$RMSPE
  MAPE <- for.perf(actual, predicted)$MAPE
  MASE <- for.perf(actual, predicted)$mase
  
  error_measures_CSR <- rbind(error_measures_CSR, data.frame(
    Series = j,
    Model = "CSR",
    MSE = MSE,
    RMSPE = RMSPE,
    MAPE = MAPE,
    MASE = MASE
  ))
}

#####
#####
#####
#####

error_measures_LAD <- data.frame()

for (j in 1:24) {
  actual <- df_list[[j]]$PUN
  predicted <- ForecastCombination_LAD[, j]
  
  MSE <- for.perf(actual, predicted)$MSE
  RMSPE <- for.perf(actual, predicted)$RMSPE
  MAPE <- for.perf(actual, predicted)$MAPE
  MASE <- for.perf(actual, predicted)$mase
  
  error_measures_LAD <- rbind(error_measures_LAD, data.frame(
    Series = j,
    Model = "LAD",
    MSE = MSE,
    RMSPE = RMSPE,
    MAPE = MAPE,
    MASE = MASE
  ))
}

#######
#######
#######
#######

error_measures_NG <- data.frame()

for (j in 1:24) {
  actual <- df_list[[j]]$PUN
  predicted <- ForecastCombination_NG[, j]
  
  MSE <- for.perf(actual, predicted)$MSE
  RMSPE <- for.perf(actual, predicted)$RMSPE
  MAPE <- for.perf(actual, predicted)$MAPE
  MASE <- for.perf(actual, predicted)$mase
  
  error_measures_NG <- rbind(error_measures_NG, data.frame(
    Series = j,
    Model = "NG",
    MSE = MSE,
    RMSPE = RMSPE,
    MAPE = MAPE,
    MASE = MASE
  ))
}

# Print the error measures for NG forecasts
error_measures_NG

########################################################################
########################################################################
########################################################################
########################################################################
########################################################################

error_measures_OLS <- data.frame()

for (j in 1:24) {
  actual <- df_list[[j]]$PUN
  predicted <- ForecastCombination_OLS[, j]
  
  MSE <- for.perf(actual, predicted)$MSE
  RMSPE <- for.perf(actual, predicted)$RMSPE
  MAPE <- for.perf(actual, predicted)$MAPE
  MASE <- for.perf(actual, predicted)$mase
  
  error_measures_OLS <- rbind(error_measures_OLS, data.frame(
    Series = j,
    Model = "OLS",
    MSE = MSE,
    RMSPE = RMSPE,
    MAPE = MAPE,
    MASE = MASE
  ))
}

# Print the error measures for OLS forecasts
error_measures_OLS

########################################################################
########################################################################
########################################################################
########################################################################
########################################################################

error_measures_SA <- data.frame()

for (j in 1:24) {
  actual <- df_list[[j]]$PUN
  predicted <- ForecastCombination_SA[, j]
  
  MSE <- for.perf(actual, predicted)$MSE
  RMSPE <- for.perf(actual, predicted)$RMSPE
  MAPE <- for.perf(actual, predicted)$MAPE
  MASE <- for.perf(actual, predicted)$mase
  
  error_measures_SA <- rbind(error_measures_SA, data.frame(
    Series = j,
    Model = "SA",
    MSE = MSE,
    RMSPE = RMSPE,
    MAPE = MAPE,
    MASE = MASE
  ))
}

# Print the error measures for SA forecasts
error_measures_SA

# Calculate the column means
means_ARIMA <- colMeans(error_measures_ARIMA[3:6])
means_ETS <- colMeans(error_measures_ETS[3:6])
means_SPLINE <- colMeans(error_measures_SPLINE[3:6])
means_RandomForest <- colMeans(error_measures_RandomForest[3:6])
means_CSR <- colMeans(error_measures_CSR[3:6])
means_LAD <- colMeans(error_measures_LAD[3:6])
means_NG <- colMeans(error_measures_NG[3:6])
means_OLS <- colMeans(error_measures_OLS[3:6])
means_SA <- colMeans(error_measures_SA[3:6])

# Create a data frame with the means
means_df <- rbind(means_ARIMA, means_ETS, means_SPLINE, means_RandomForest, means_CSR, means_LAD, means_NG, means_OLS, means_SA)
rownames(means_df) <- c("ARIMA", "ETS", "SPLINE", "RandomForest", "CSR", "LAD", "NG", "OLS", "SA")

# Convert the data frame to a LaTeX table
library(knitr)
latex_table <- kable(means_df, "latex", booktabs = TRUE)

# Print the LaTeX table
print(latex_table)



#################################
#################################
#################################
# Objective function to minimize
calc_rmse <- function(weights, actual, forecasts) {
  weights <- weights / sum(weights)  # Normalize weights
  weighted_forecast <- forecasts %*% weights
  rmse <- sqrt(mean((actual - weighted_forecast)^2))
  return(rmse)
}

# Function to find the best combination
find_best_combination <- function(actual, forecasts) {
  num_forecasts <- ncol(forecasts)
  
  # Initialize weights
  initial_weights <- rep(1/num_forecasts, num_forecasts)
  
  # Use optim function to minimize rmse
  result <- optim(initial_weights, calc_rmse, method = "L-BFGS-B", lower = rep(0, num_forecasts), upper = rep(1, num_forecasts), actual = actual, forecasts = forecasts)
  
  # Normalize the weights to sum to 1
  optimal_weights <- result$par / sum(result$par)
  
  # Return the optimal weights
  return(optimal_weights)
}

for (hour in 1 :24){
  forecasts <- cbind(ForecastArima[, hour], ForecastETS[, hour], ForecastRandomForest[, hour], ForecastSPLINE[, hour])
  actual <- df2_list[[hour]]$PUN
  best_weights <- find_best_combination(actual, forecasts)
  print(best_weights)}

library(ggplot2)


library(ggplot2)

# Create an empty list to store the plots
plot_list <- list()

for (hour in 1:24) {
  
  data1 <- data.frame(
    Date = df2_list[[hour]]$Date, 
    Actual = df2_list[[hour]]$PUN, 
    ARIMA = ForecastArima[, hour],
    ETS = ForecastETS[, hour],
    RandomForest = ForecastRandomForest[, hour],
    SPLINE = ForecastSPLINE[, hour],
    CSR = ForecastCombination_CSR[, hour],
    LAD = ForecastCombination_LAD[, hour],
    NG = ForecastCombination_NG[, hour],
    OLS = ForecastCombination_OLS[, hour],
    SA = ForecastCombination_SA[, hour]
  )
  
  # Transform data from wide to long format
  data_long <- reshape2::melt(data1, id.vars = "Date")
  
  plot_list[[hour]] <- ggplot(data_long, aes(x = Date, y = value, color = variable, linetype = variable)) +
    geom_line() +
    scale_color_manual(values = c("Actual" = "black",
                                  "ARIMA" = "blue",
                                  "ETS" = "red",
                                  "RandomForest" = "green",
                                  "SPLINE" = "purple",
                                  "CSR" = "orange",
                                  "LAD" = "brown",
                                  "NG" = "pink",
                                  "OLS" = "gray",
                                  "SA" = "yellow")) +
    scale_linetype_manual(values = c("Actual" = "solid",
                                     "ARIMA" = "solid",
                                     "ETS" = "solid",
                                     "RandomForest" = "solid",
                                     "SPLINE" = "solid",
                                     "CSR" = "dashed",
                                     "LAD" = "dashed",
                                     "NG" = "dashed",
                                     "OLS" = "dashed",
                                     "SA" = "dashed")) +
    labs(title = paste("Hour", hour),
         y = "Value",
         color = "Forecast",
         linetype = "Forecast") +
    theme_minimal() +
    theme(legend.position = "bottom")
}

plot_list[[18]]


# Create an empty list to store the plots
plot_list_combinations <- list()

for (hour in 1:24) {
  
  data1 <- data.frame(
    Date = df2_list[[hour]]$Date, 
    Actual = df2_list[[hour]]$PUN, 
    CSR = ForecastCombination_CSR[, hour],
    LAD = ForecastCombination_LAD[, hour],
    NG = ForecastCombination_NG[, hour],
    OLS = ForecastCombination_OLS[, hour],
    SA = ForecastCombination_SA[, hour]
  )
  
  # Transform data from wide to long format
  data_long <- reshape2::melt(data1, id.vars = "Date")
  
  plot_list_combinations[[hour]] <- ggplot(data_long, aes(x = Date, y = value, color = variable, linetype = variable)) +
    geom_line() +
    scale_color_manual(values = c("Actual" = "black",
                                  "CSR" = "blue",
                                  "LAD" = "red",
                                  "NG" = "green",
                                  "OLS" = "purple",
                                  "SA" = "orange")) +
    scale_linetype_manual(values = c("Actual" = "solid",
                                     "CSR" = "dashed",
                                     "LAD" = "dashed",
                                     "NG" = "dashed",
                                     "OLS" = "dashed",
                                     "SA" = "dashed")) +
    labs(title = paste("Hour", hour),
         y = "Value",
         color = "Forecast",
         linetype = "Forecast") +
    theme_minimal() +
    theme(legend.position = "bottom")
}

plot_list_combinations[[1]]


library(ggplot2)
library(dplyr)

# Preparing data for ggplot
data_list <- list()
for(hour in 1:24){
  df <- data.frame(
    Time = 1:length(df2_list[[hour]]$PUN),
    Actual = df2_list[[hour]]$PUN,
    Forecast_LAD = ForecastCombination_LAD[, hour],
    Forecast_ETS = ForecastETS[, hour]
  )
  
  df <- df %>%
    gather(key = "Type", value = "Value", -Time)
  
  data_list[[hour]] <- df
}

# Plotting
plots_list <- list()
for(hour in 1:24){
  plots_list[[hour]] <- ggplot(data_list[[hour]], aes(x = Time, y = Value, color = Type)) +
    geom_line() +
    scale_color_manual(values = c("Actual" = "black", "Forecast_LAD" = "blue", "Forecast_ETS" = "red")) +
    labs(title = paste("Hour", hour), y = "Value", x = "Time") +
    theme_minimal()
}


library(ggplot2)
library(dplyr)
library(RColorBrewer)

# Define custom color palette
my_palette <- brewer.pal(3, "Set1")

# Preparing data for ggplot
data_list <- list()
for(hour in 1:24){
  df <- data.frame(
    Time = 1:length(df2_list[[hour]]$PUN),
    Actual = df2_list[[hour]]$PUN,
    Forecast_LAD = ForecastCombination_LAD[, hour],
    Forecast_ETS = ForecastETS[, hour]
  )
  
  df <- df %>%
    gather(key = "Type", value = "Value", -Time)
  
  data_list[[hour]] <- df
}

# Plotting
plots_list <- list()
for(hour in 1:24){
  plots_list[[hour]] <- ggplot(data_list[[hour]], aes(x = Time, y = Value, color = Type)) +
    geom_line(size = 1) +
    scale_color_manual(values = my_palette) +
    labs(title = paste("Hour", hour), y = "Value", x = "Time") +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "#FCFCFC"),
      panel.grid.major = element_line(color = "black", size = 0.2),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(size = 1),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
}


# Storing the plots in a list
plots_list[[1]]
grid.arrange(grobs = plots_list[1:6], ncol = 3)
grid.arrange(grobs = plots_list[7:12], ncol = 3)
grid.arrange(grobs = plots_list[13:18], ncol = 3)
grid.arrange(grobs = plots_list[19:24], ncol = 3)



library(forecast)
library(dplyr)
library(forecast)
library(dplyr)

# List of forecast combinations

forecast_combinations <- list(
  #CSR = ForecastCombination_CSR
  LAD = ForecastCombination_LAD
  #NG = ForecastCombination_NG,
  #OLS = ForecastCombination_OLS
)

# List of individual forecasts
forecast_individuals <- list(
  #Arima = ForecastArima
  #ETS = ForecastETS
  #RandomForest = ForecastRandomForest
  SPLINE = ForecastSPLINE
)

# Create an empty data frame to store results
results <- data.frame()

# For each hour
for (hour in 1:24) {
  
  # For each combination method
  for (combination_name in names(forecast_combinations)) {
    
    # For each individual method
    for (individual_name in names(forecast_individuals)) {
      
      # Ensure that all the series have the same length
      n <- min(length(df2_list[[hour]]$PUN),
               length(forecast_combinations[[combination_name]][, hour]),
               length(forecast_individuals[[individual_name]][, hour]))
      
      # Create time series objects
      x1 <- ts(forecast_combinations[[combination_name]][, hour][1:n])
      x2 <- ts(forecast_individuals[[individual_name]][, hour][1:n])
      
      # Perform DM test
      dm_test <- dm.test((x1-df2_list[[hour]]$PUN), (x2-df2_list[[hour]]$PUN), h = 1, power = 1)
      
      # Create a data frame row with the results
      result_row <- data.frame(
        Hour = hour,
        Combination = combination_name,
        Individual = individual_name,
        DM_Stat = dm_test$statistic,
        p_value = dm_test$p.value
      )
      
      # Append to the results
      results <- bind_rows(results, result_row)
    }
  }
}

# Print the results
print(results)

# To create a LaTeX table from the `results` dataframe, you can use the `xtable` package:
library(xtable)
print(xtable(results), type = "latex")


#####################################################################
#####################################################################
#####################################################################
# Dates to be displayed on x-axis
dates <- c("2023-02-01", "2023-03-01", "2023-04-01", "2023-04-30")
dates <- c("02-01", "03-01", "04-01", "04-30")
# Corresponding time points in your data
breaks <- round(seq(from = 1, to = 89, length.out = length(dates)))

# Modify the plot code as follows
plots_list <- list()
for(hour in 1:24){
  plots_list[[hour]] <- ggplot(data_list[[hour]], aes(x = Time, y = Value, color = Type)) +
    geom_line(size = 0.8) +
    scale_color_manual(values = my_palette) +
    labs(title = paste("Hour", hour), y = "Value", x = "Date") +  # Changed x label to "Date"
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "#FCFCFC"),
      panel.grid.major = element_line(color = "black", size = 0.2),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(size = 1),
      legend.position = "bottom",
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    ) +
    scale_x_continuous(breaks = breaks, labels = dates)  # Manually specified breaks and labels for x-axis
}

grid.arrange(grobs = plots_list[1:6], ncol = 3)
grid.arrange(grobs = plots_list[7:12], ncol = 3)
grid.arrange(grobs = plots_list[13:18], ncol = 3)
grid.arrange(grobs = plots_list[19:24], ncol = 3)


######################################################################
######################################################################
######################################################################
######################################################################
calc_mse <- function(weights, actual, forecasts) {
  weights <- weights / sum(weights)  # Normalize weights
  weighted_forecast <- forecasts %*% weights
  mse=for.perf(actual,weighted_forecast)$MSE
  return(mse)
}

# Function to find the best combination
find_best_combination <- function(actual, forecasts) {
  num_forecasts <- ncol(forecasts)
  
  # Initialize weights
  initial_weights <- rep(1/num_forecasts, num_forecasts)
  
  # Use optim function to minimize rmse
  result <- optim(initial_weights, calc_mse, method = "L-BFGS-B", lower = rep(0, num_forecasts), upper = rep(1, num_forecasts), actual = actual, forecasts = forecasts)
  
  # Normalize the weights to sum to 1
  optimal_weights <- result$par / sum(result$par)
  
  # Return the optimal weights
  return(optimal_weights)
}

# Create an empty data frame to store the results
weighted_forecastsMSE <- data.frame(NA)
dim(weighted_forecastsMSE)
for (hour in 1:24) {
  forecasts <- cbind(ForecastArima[, hour], ForecastETS[, hour], ForecastRandomForest[, hour], ForecastSPLINE[, hour])
  actual <- df_list[[hour]]$PUN
  best_weights <- find_best_combination(actual, forecasts)
  
  # Calculate the weighted forecasts
  weighted_forecast <- forecasts * best_weights
  
  # Add the weighted forecasts to the dataframe
  weighted_forecastsMSE <- cbind(weighted_forecastsMSE, rowSums(weighted_forecast))
}

dim(weighted_forecastsMSE)
weighted_forecastsMSE=weighted_forecastsMSE[,2:25]
# Set the column names
colnames(weighted_forecastsMSE) <- paste0("Hour_", 1:24)

# Print the weighted forecasts
print(weighted_forecastsMSE)

error_measures_weighted <- data.frame()

for (j in 1:24) {
  actual <- df_list[[j]]$PUN
  predicted <- weighted_forecastsMSE[, j]
  
  MSE <- for.perf(actual, predicted)$MSE
  RMSPE <- for.perf(actual, predicted)$RMSPE
  MAPE <- for.perf(actual, predicted)$MAPE
  MASE <- for.perf(actual, predicted)$mase
  
  error_measures_weighted <- rbind(error_measures_weighted, data.frame(
    Series = j,
    Model = "Weighted",
    MSE = MSE,
    RMSPE = RMSPE,
    MAPE = MAPE,
    MASE = MASE
  ))
}

# Print the error measures for weighted forecasts
error_measures_weighted
colMeans(error_measures_weighted[3:6])


#######################################################################
#######################################################################
#######################################################################
#######################################################################
#######################################################################


calc_mase <- function(weights, actual, forecasts) {
  weights <- weights / sum(weights)  # Normalize weights
  weighted_forecast <- forecasts %*% weights
  mase=for.perf(actual,weighted_forecast)$mase
  return(mase)
}

# Function to find the best combination
find_best_combination <- function(actual, forecasts) {
  num_forecasts <- ncol(forecasts)
  
  # Initialize weights
  initial_weights <- rep(1/num_forecasts, num_forecasts)
  
  # Use optim function to minimize rmase
  result <- optim(initial_weights, calc_mase, method = "L-BFGS-B", lower = rep(0, num_forecasts), upper = rep(1, num_forecasts), actual = actual, forecasts = forecasts)
  
  # Normalize the weights to sum to 1
  optimal_weights <- result$par / sum(result$par)
  
  # Return the optimal weights
  return(optimal_weights)
}

# Create an empty data frame to store the results
weighted_forecastsmase <- data.frame(NA)
dim(weighted_forecastsmase)
for (hour in 1:24) {
  forecasts <- cbind(ForecastArima[, hour], ForecastETS[, hour], ForecastRandomForest[, hour], ForecastSPLINE[, hour])
  actual <- df_list[[hour]]$PUN
  best_weights <- find_best_combination(actual, forecasts)
  
  # Calculate the weighted forecasts
  weighted_forecast <- forecasts * best_weights
  print(best_weights)
  # Add the weighted forecasts to the dataframe
  weighted_forecastsmase <- cbind(weighted_forecastsmase, rowSums(weighted_forecast))
}

dim(weighted_forecastsmase)
weighted_forecastsmase=weighted_forecastsmase[,2:25]
# Set the column names
colnames(weighted_forecastsmase) <- paste0("Hour_", 1:24)


error_measures_weighted <- data.frame()

for (j in 1:24) {
  actual <- df_list[[j]]$PUN
  predicted <- weighted_forecastsmase[, j]
  
  MSE <- for.perf(actual, predicted)$MSE
  RMSPE <- for.perf(actual, predicted)$RMSPE
  MAPE <- for.perf(actual, predicted)$MAPE
  MASE <- for.perf(actual, predicted)$mase
  
  error_measures_weighted <- rbind(error_measures_weighted, data.frame(
    Series = j,
    Model = "Weighted",
    MSE = MSE,
    RMSPE = RMSPE,
    MAPE = MAPE,
    MASE = MASE
  ))
}

