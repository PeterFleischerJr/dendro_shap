############################################################
# RF Metrics Calculation + Plotting Script                  #
############################################################

# 1) Load required libraries
library(readr)       # read_csv, write_csv
library(dplyr)       # data manipulation
library(tidyr)       # pivot_longer
library(lubridate)   # date parsing, year()
library(ranger)      # random forest
library(fastshap)    # SHAP values
library(ggplot2)     # plotting
library(RColorBrewer) # palettes
library(Metrics)     # rmse()

# 2) Read data
# Allow supplying a custom CSV path via command line. Fallback to the
# bundled data file if no argument is provided.
args <- commandArgs(trailingOnly = TRUE)
path_to_csv <- if (length(args) > 0) args[1] else "wd_2015_ok.csv"
if (!file.exists(path_to_csv)) stop("File not found: ", path_to_csv)
raw_data <- read_csv(path_to_csv)

# 3) Preprocess: pivot and daily aggregation
daily_data <- raw_data %>%
  mutate(
    DateTime = as.POSIXct(paste(Date, Time), format = "%m/%d/%y %H:%M"),
    Date     = as.Date(DateTime),
    Year     = year(DateTime)
  ) %>%
  pivot_longer(
    cols      = starts_with("JD"),
    names_to  = "Individual",
    values_to = "Value"
  ) %>%
  drop_na(Date, Individual, TV, swp, VPD, Value) %>%
  group_by(Date, Year, Individual) %>%
  summarise(
    TV    = mean(TV, na.rm = TRUE),
    swp   = mean(swp, na.rm = TRUE),
    VPD   = mean(VPD, na.rm = TRUE),
    Value = mean(Value, na.rm = TRUE),
    .groups = "drop"
  )

# 4) Metrics calculation function
calculate_metrics <- function(actual, predicted, num_params) {
  n      <- length(actual)
  resid  <- actual - predicted
  RSS    <- sum(resid^2)
  TSS    <- sum((actual - mean(actual))^2)
  R2     <- 1 - (RSS/TSS)
  RMSE   <- sqrt(mean(resid^2))
  AIC    <- n * log(RSS / n) + 2 * num_params
  baseline <- mean(actual)
  tibble(R2 = R2, RMSE = RMSE, AIC = AIC, Baseline = baseline)
}

# 5) Initialize container for metrics
metrics_list <- list()

# 6) LOOP: Individuals across ALL years -----------------------
for (ind in unique(daily_data$Individual)) {
  df <- filter(daily_data, Individual == ind)
  if (nrow(df) < 10) next
  
  rf   <- ranger(Value ~ TV + swp + VPD, data = df, importance = "impurity")
  preds <- predict(rf, df)$predictions
  m <- calculate_metrics(df$Value, preds, length(rf$variable.importance)) %>%
    mutate(Type = "Individual_AllYears", Group = ind)
  metrics_list[[paste0("IndAll_", ind)]] <- m
}

# 7) LOOP: Years across ALL individuals ------------------------
for (yr in unique(daily_data$Year)) {
  df <- filter(daily_data, Year == yr)
  if (nrow(df) < 10) next
  
  rf   <- ranger(Value ~ TV + swp + VPD, data = df, importance = "impurity")
  preds <- predict(rf, df)$predictions
  m <- calculate_metrics(df$Value, preds, length(rf$variable.importance)) %>%
    mutate(Type = "Year_AllIndividuals", Group = as.character(yr))
  metrics_list[[paste0("YearAll_", yr)]] <- m
}

# 8) LOOP: Each (Individual, Year) combo ----------------------
# Also open PDF for SHAP plots
pdf("All_Individual_Year_Plots.pdf", width = 10, height = 6)
for (ind in unique(daily_data$Individual)) {
  for (yr in unique(daily_data$Year)) {
    df <- daily_data %>% filter(Individual == ind, Year == yr) %>% arrange(Date)
    if (nrow(df) < 10) next
    
    # Random Forest & metrics
    rf    <- ranger(Value ~ TV + swp + VPD, data = df)
    preds <- predict(rf, df)$predictions
    m <- calculate_metrics(df$Value, preds, length(rf$variable.importance)) %>%
      mutate(Type = "Individual_Year", Individual = ind, Year = yr)
    metrics_list[[paste0("IndYr_", ind, "_", yr)]] <- m
    
    # SHAP values & plot
    shap_vals <- fastshap::explain(
      object       = rf,
      X            = df %>% select(TV, swp, VPD),
      pred_wrapper = function(m, newdata) predict(m, data = newdata)$predictions,
      nsim         = 10,
      adjust       = TRUE
    )
    shap_df <- bind_cols(Date = df$Date, as.data.frame(shap_vals))
    shap_long <- pivot_longer(shap_df, cols = c("TV","swp","VPD"),
                              names_to = "Predictor", values_to = "SHAP")
    
    baseline <- m$Baseline
    p <- ggplot() +
      geom_col(data = shap_long, aes(x = Date, y = SHAP, fill = Predictor),
               position = "stack", alpha = 0.8) +
      geom_line(data = df, aes(x = Date, y = Value), color = "black", size = 1) +
      geom_hline(yintercept = baseline, linetype = "dashed", color = "blue") +
      scale_fill_brewer(palette = "Set1") +
      labs(title    = paste(ind, "—", yr),
           subtitle = paste0("R²=", round(m$R2, 3),
                             " | RMSE=", round(m$RMSE, 3)),
           x        = "Date", y = "SHAP & Value") +
      theme_minimal() + theme(legend.position = "top",
                              axis.text.x = element_text(angle = 45, hjust = 1))
    print(p)
  }
}
dev.off()

# 9) Combine metrics & save CSV ------------------------------
metrics_results <- bind_rows(metrics_list)
write_csv(metrics_results, "model_metrics_results.csv")
print(metrics_results)

cat("Done.\n", 
    "• CSV -> model_metrics_results.csv\n", 
    "• PDF -> All_Individual_Year_Plots.pdf\n")
