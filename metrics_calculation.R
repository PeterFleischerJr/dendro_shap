# Metrics extraction script for dendro data

# Load required libraries
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(ranger)
library(Metrics)

# Path to the CSV data file
path_to_csv <- "wd_2015_ok.csv"
if (!file.exists(path_to_csv)) {
  stop("File not found: ", path_to_csv)
}

# Read the data
raw_data <- read_csv(path_to_csv)

# Data preprocessing
processed <- raw_data %>%
  mutate(DateTime = as.POSIXct(paste(Date, Time), format = "%m/%d/%y %H:%M"),
         Date = as.Date(DateTime)) %>%
  pivot_longer(
    cols = starts_with("JD"),
    names_to = "Individual",
    values_to = "Value"
  ) %>%
  drop_na(Date, Individual, TV, swp, VPD, Value) %>%
  group_by(Date, Individual) %>%
  summarise(
    across(c(TV, swp, VPD, Value), mean, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Year = year(Date))

# Helper function to compute metrics
calc_metrics <- function(actual, predicted, n_params) {
  n <- length(actual)
  residuals <- actual - predicted
  rss <- sum(residuals^2)
  tss <- sum((actual - mean(actual))^2)
  r2 <- 1 - rss / tss
  rmse_val <- sqrt(mean(residuals^2))
  aic_val <- n * log(rss / n) + 2 * n_params
  list(R2 = r2, RMSE = rmse_val, AIC = aic_val, Baseline = mean(actual))
}

results <- data.frame()
individuals <- unique(processed$Individual)
years <- unique(processed$Year)

# Loop over individuals and years
for (ind in individuals) {
  for (yr in years) {
    subset_data <- processed %>% filter(Individual == ind, Year == yr)
    if (nrow(subset_data) < 10) {
      next
    }
    rf_model <- ranger(
      Value ~ TV + swp + VPD,
      data = subset_data,
      num.trees = 100,
      importance = "impurity"
    )
    predictions <- predict(rf_model, subset_data)$predictions
    m <- calc_metrics(
      subset_data$Value,
      predictions,
      length(rf_model$variable.importance)
    )
    results <- rbind(
      results,
      data.frame(
        Individual = ind,
        Year = yr,
        Baseline = m$Baseline,
        R2 = m$R2,
        RMSE = m$RMSE,
        AIC = m$AIC
      )
    )
  }
}

# Save metrics to CSV
write_csv(results, "model_metrics_results.csv")

cat("Metrics extraction complete. Results saved to model_metrics_results.csv\n")
