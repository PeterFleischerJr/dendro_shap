

# Load necessary libraries
library(ggplot2)
library(RColorBrewer)
library(patchwork) # For combining plots
library(readr)
library(dplyr)
library(ranger)
library(fastshap)
library(lubridate)
library(doParallel)
library(tidyr)

# Setup parallel processing
cores <- parallel::detectCores() - 1
cl <- makeCluster(cores)
registerDoParallel(cl)

# Load data
data <- read_csv("/home/peter/Downloads/wd_2015_ok.csv") #wd_2015_custom_avg_updated

# Define predictors and targets
predictors <- c("TV", "swp", "VPD")
targets <- c("JD1", "JD2", "JD3", "JD4", "JD5")  # adjusted for all targets JD1-JD5

# Parse datetime and calculate daily averages
daily_data <- data %>%
  mutate(DateTime = as.POSIXct(paste(Date, Time), format='%m/%d/%y %H:%M'),
         Date = as.Date(DateTime)) %>%
  drop_na(DateTime, all_of(predictors), all_of(targets)) %>%
  group_by(Date) %>%
  summarise(across(all_of(c(predictors, targets)), mean, na.rm = TRUE)) %>%
  ungroup()






# Start PDF output
pdf("/home/peter/Downloads/All_SHAP_Combined_Plots_3env_NoScaling.pdf", width = 14, height = 16)

for (target in targets) {
  
  cat("\nGenerating visualization for target:", target, "\n")
  
  daily_target_data <- daily_data %>% drop_na(all_of(target)) %>% mutate(Year = year(Date))
  years <- unique(daily_target_data$Year)
  
  for (current_year in years) {
    
    cat("\nPlotting Year:", current_year, "for target:", target, "\n")
    
    yearly_data <- daily_target_data %>% filter(Year == current_year)
    
    rf_model <- ranger(
      formula    = as.formula(paste(target, "~ .")),
      data       = yearly_data %>% select(all_of(target), all_of(predictors)),
      num.trees  = 100,
      importance = "impurity"
    )
    
    predictions <- predict(rf_model, yearly_data)$predictions
    R2_value <- cor(yearly_data[[target]], predictions)^2
    
    shap_values <- fastshap::explain(
      object       = rf_model,
      X            = yearly_data %>% select(all_of(predictors)),
      pred_wrapper = function(m, newdata) predict(m, data = newdata)$predictions,
      nsim         = 10,
      parallel     = FALSE,
      adjust       = TRUE
    )
    
    shap_df <- yearly_data %>% select(Date) %>% bind_cols(shap_values)
    
    daily_sums <- shap_df %>%
      mutate(Day = Date) %>%   
      group_by(Day) %>%
      summarise(across(all_of(predictors), sum, na.rm = TRUE)) %>%
      arrange(Day)
    
    daily_long_df <- daily_sums %>%
      pivot_longer(cols = all_of(predictors), names_to = "Variable", values_to = "SHAP")
    
    actual_target_df <- yearly_data %>%
      select(Day = Date, ActualTarget = all_of(target))
    
    baseline <- mean(actual_target_df$ActualTarget, na.rm = TRUE)
    
    # Upper plot: SHAP, Actual, and Baseline (No standard scaling)
    p_shap <- ggplot() +
      geom_bar(data = daily_long_df, aes(x = Day, y = SHAP, fill = Variable),
               stat = "identity", position = "stack", alpha = 0.8) +
      geom_line(data = actual_target_df, aes(x = Day, y = ActualTarget),
                color = "black", linewidth = 1, group = 1) +
      geom_hline(yintercept = baseline, color = "blue", linetype = "dashed", linewidth = 1) +
      scale_fill_brewer(palette = "Set1") +
      labs(title = paste("Daily SHAP Values and", target, "for Year", current_year),
           subtitle = paste("Baseline (mean", target, "):", round(baseline, 3), "| R²:", round(R2_value, 3)),
           x = "",
           y = paste("SHAP Sum and", target),
           fill = "Predictors") +
      theme_minimal() +
      theme(axis.text.x = element_blank(),
            legend.position = "top")
    
    # Environmental variables plotted separately without scaling
    p_swp <- ggplot(yearly_data, aes(x = Date, y = swp)) +
      geom_line(color = "darkgreen") +
      labs(title = "Soil Water Potential (swp)", x = "", y = "swp") +
      theme_minimal() +
      theme(axis.text.x = element_blank())
    
    p_tv <- ggplot(yearly_data, aes(x = Date, y = TV)) +
      geom_line(color = "red") +
      labs(title = "Air Temperature (TV)", x = "", y = "°C") +
      theme_minimal() +
      theme(axis.text.x = element_blank())
    
    p_vpd <- ggplot(yearly_data, aes(x = Date, y = VPD)) +
      geom_line(color = "purple") +
      labs(title = "Vapor Pressure Deficit (VPD)", x = "Day", y = "VPD") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    # Combine all plots vertically without shared y-axis
    combined_plot <- p_shap / p_swp / p_tv / p_vpd +
      plot_layout(heights = c(3, 1, 1, 1))
    
    print(combined_plot)
  }
}

# Close PDF output
dev.off()









############################################################
# 1) Load libraries
############################################################
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(ranger)
library(fastshap)
library(ggplot2)
library(RColorBrewer)

############################################################
# 2) Read your CSV
############################################################
# Adjust path to your actual file location
path_to_csv <- "/home/peter/Downloads/wd_2015_ok.csv"
if (!file.exists(path_to_csv)) {
  stop("File not found: ", path_to_csv)
}

raw_data <- read_csv(path_to_csv)

############################################################
# 3) Parse Date + Time => POSIXct => daily Date
############################################################
# Your sample: "4/1/15" for Date, "00:20" for Time => %m/%d/%y %H:%M
raw_data <- raw_data %>%
  mutate(
    DateTime = as.POSIXct(
      paste(Date, Time),
      format = "%m/%d/%y %H:%M"
    ),
    Date = as.Date(DateTime)
  )

############################################################
# 4) Pivot JD1..JD5 from wide to long
#    so we get columns: Date, Time, TV, swp, VPD, Individual, Value
############################################################
# If you also want JD6, just add it to 'cols = c(JD1,JD2,JD3,JD4,JD5,JD6)'
long_data <- raw_data %>%
  pivot_longer(
    cols = c(JD1, JD2, JD3, JD4, JD5),  # ignoring JD6
    names_to = "Individual",            # becomes JD1..JD5
    values_to = "Value"
  )

############################################################
# 5) Aggregate to daily means by (Date, Individual)
############################################################
daily_data <- long_data %>%
  drop_na(Date, Individual, TV, swp, VPD, Value) %>%
  group_by(Date, Individual) %>%
  summarise(
    TV    = mean(TV, na.rm=TRUE),
    swp   = mean(swp, na.rm=TRUE),
    VPD   = mean(VPD, na.rm=TRUE),
    Value = mean(Value, na.rm=TRUE),  # your measurement
    .groups = "drop"
  )

cat("Daily data shape => rows:", nrow(daily_data), ", cols:", ncol(daily_data), "\n")
head(daily_data)

############################################################
# 6) Model 1: Individuals across all years
############################################################
# We'll add a 'Year' column from the daily 'Date'
daily_data <- daily_data %>%
  mutate(Year = year(Date))

cat("\n=== MODEL 1: INDIVIDUALS ACROSS ALL YEARS ===\n")

# We do 'Value ~ TV + swp + VPD' for each 'Individual'
# The script loops over JD1..JD5 as 'Individual'
indivs <- unique(daily_data$Individual)

for (ind in indivs) {
  cat("\nFitting across all years for Individual:", ind, "\n")
  
  subset_data <- daily_data %>%
    filter(Individual == ind)
  
  if (nrow(subset_data) == 0) {
    cat("No rows for individual:", ind, "\n")
    next
  }
  
  # Random Forest => Value ~ TV + swp + VPD
  rf_model <- ranger(
    formula    = Value ~ TV + swp + VPD,
    data       = subset_data,
    num.trees  = 100,
    importance = "impurity"
  )
  
  preds    <- predict(rf_model, subset_data)$predictions
  R2_value <- cor(subset_data$Value, preds)^2
  
  cat("R² for individual", ind, "=", round(R2_value, 3), "\n")
  
  # SHAP
  shap_vals <- fastshap::explain(
    object = rf_model,
    X = subset_data %>% select(TV, swp, VPD),
    pred_wrapper = function(m, newdata) predict(m, data=newdata)$predictions,
    nsim         = 10,
    parallel     = FALSE,
    adjust       = TRUE
  )
  
  # Prepare for plotting
  shap_df <- subset_data %>%
    select(Date) %>%
    bind_cols(shap_vals)
  
  shap_long <- shap_df %>%
    pivot_longer(
      cols = c("TV","swp","VPD"),
      names_to = "Predictor",
      values_to = "SHAP"
    )
  
  # Plot
  ggplot(shap_long, aes(x = Date, y = SHAP, fill = Predictor)) +
    geom_bar(stat="identity", position="stack") +
    labs(
      title    = paste("Individual", ind, "| SHAP over time"),
      subtitle = paste("R² =", round(R2_value, 3)),
      x        = "Date",
      y        = "SHAP Value",
      fill     = "Predictor"
    ) +
    scale_fill_brewer(palette="Set2") +
    theme_minimal() +
    theme(
      axis.text.x     = element_text(angle=45, hjust=1),
      legend.position = "top"
    )
}

############################################################
# 7) Model 2: Each Year across all Individuals
############################################################
cat("\n=== MODEL 2: YEARS ACROSS ALL INDIVIDUALS ===\n")

yrs <- unique(daily_data$Year)
for (yr in yrs) {
  cat("\nFitting for Year:", yr, "\n")
  
  subset_data <- daily_data %>%
    filter(Year == yr)
  
  if (nrow(subset_data) == 0) {
    cat("No rows for year", yr, "\n")
    next
  }
  
  # Random Forest => Value ~ TV + swp + VPD
  rf_model <- ranger(
    formula    = Value ~ TV + swp + VPD,
    data       = subset_data,
    num.trees  = 100,
    importance = "impurity"
  )
  
  preds    <- predict(rf_model, subset_data)$predictions
  R2_value <- cor(subset_data$Value, preds)^2
  cat("R² for year", yr, "=", round(R2_value, 3), "\n")
  
  # SHAP
  shap_vals <- fastshap::explain(
    object = rf_model,
    X = subset_data %>% select(TV, swp, VPD),
    pred_wrapper = function(m, newdata) predict(m, data=newdata)$predictions,
    nsim         = 10,
    parallel     = FALSE,
    adjust       = TRUE
  )
  
  shap_df <- subset_data %>% select(Date) %>% bind_cols(shap_vals)
  
  shap_long <- shap_df %>%
    pivot_longer(
      cols = c("TV","swp","VPD"),
      names_to = "Predictor",
      values_to = "SHAP"
    )
  
  # Plot
  ggplot(shap_long, aes(x=Date, y=SHAP, fill=Predictor)) +
    geom_bar(stat="identity", position="stack") +
    labs(
      title    = paste("Year", yr, "| SHAP over time"),
      subtitle = paste("R² =", round(R2_value, 3)),
      x        = "Date",
      y        = "SHAP Value",
      fill     = "Predictor"
    ) +
    scale_fill_brewer(palette="Set1") +
    theme_minimal() +
    theme(
      axis.text.x     = element_text(angle=45, hjust=1),
      legend.position = "top"
    )
}

cat("\n=== ALL DONE! ===\n")








############################################################
# 1) Load libraries
############################################################
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(ranger)
library(fastshap)
library(ggplot2)
library(RColorBrewer)

############################################################
# 2) Read your CSV (adjust path as needed)
############################################################
path_to_csv <- "/home/peter/Downloads/wd_2015_ok.csv"
if (!file.exists(path_to_csv)) {
  stop("File not found: ", path_to_csv)
}
raw_data <- read_csv(path_to_csv)

############################################################
# 3) Parse Date + Time => POSIXct => daily Date
############################################################
raw_data <- raw_data %>%
  mutate(
    DateTime = as.POSIXct(paste(Date, Time), format = "%m/%d/%y %H:%M"),
    Date = as.Date(DateTime)
  )

############################################################
# 4) Pivot JD1..JD5 from wide to long
############################################################
long_data <- raw_data %>%
  pivot_longer(
    cols = c(JD1, JD2, JD3, JD4, JD5),
    names_to = "Individual",
    values_to = "Value"
  )

############################################################
# 5) Aggregate to daily means by (Date, Individual)
############################################################
daily_data <- long_data %>%
  drop_na(Date, Individual, TV, swp, VPD, Value) %>%
  group_by(Date, Individual) %>%
  summarise(
    TV    = mean(TV, na.rm=TRUE),
    swp   = mean(swp, na.rm=TRUE),
    VPD   = mean(VPD, na.rm=TRUE),
    Value = mean(Value, na.rm=TRUE),
    .groups = "drop"
  )

cat("Daily data shape => rows:", nrow(daily_data), ", cols:", ncol(daily_data), "\n")
head(daily_data)

############################################################
# Define your predictor columns and formula
############################################################
predictors <- c("TV", "swp", "VPD")  # what you use for modeling

############################################################
# MODEL 1: Individuals across all years
############################################################
cat("\n=== MODEL 1: INDIVIDUALS ACROSS ALL YEARS ===\n")

# Output plots to PDF
pdf("Model1_Individuals.pdf", width = 12, height = 7)

# Identify unique individuals (JD1..JD5 in the pivot)
indivs <- unique(daily_data$Individual)

for (ind in indivs) {
  cat("\nFitting across all years for Individual:", ind, "\n")
  
  subset_data <- daily_data %>%
    filter(Individual == ind)
  
  if (nrow(subset_data) == 0) {
    cat("No rows for individual:", ind, "\n")
    next
  }
  
  # Fit Random Forest: Value ~ TV + swp + VPD
  rf_model <- ranger(
    formula    = Value ~ TV + swp + VPD,
    data       = subset_data,
    num.trees  = 100,
    importance = "impurity"
  )
  
  preds    <- predict(rf_model, subset_data)$predictions
  R2_value <- cor(subset_data$Value, preds)^2
  
  cat("R² for individual", ind, "=", round(R2_value, 3), "\n")
  
  # SHAP
  shap_vals <- fastshap::explain(
    object = rf_model,
    X = subset_data %>% select(TV, swp, VPD),
    pred_wrapper = function(m, newdata) predict(m, data=newdata)$predictions,
    nsim         = 10,
    parallel     = FALSE,
    adjust       = TRUE
  )
  
  # Prepare for plotting
  shap_df <- subset_data %>%
    select(Date) %>%
    bind_cols(shap_vals)
  
  shap_long <- shap_df %>%
    pivot_longer(
      cols = c("TV","swp","VPD"),
      names_to = "Predictor",
      values_to = "SHAP"
    )
  
  # Plot
  p <- ggplot(shap_long, aes(x = Date, y = SHAP, fill = Predictor)) +
    geom_bar(stat="identity", position="stack") +
    labs(
      title    = paste("Individual", ind, "| SHAP over time"),
      subtitle = paste("R² =", round(R2_value, 3)),
      x        = "Date",
      y        = "SHAP Value",
      fill     = "Predictor"
    ) +
    scale_fill_brewer(palette="Set2") +
    theme_minimal() +
    theme(
      axis.text.x     = element_text(angle=45, hjust=1),
      legend.position = "top"
    )
  
  print(p)
}

dev.off()  # close Model1 PDF

############################################################
# MODEL 2: Each Year across all Individuals
############################################################
cat("\n=== MODEL 2: YEARS ACROSS ALL INDIVIDUALS ===\n")

# Add a 'Year' column
daily_data <- daily_data %>%
  mutate(Year = year(Date))

pdf("Model2_Years.pdf", width = 12, height = 7)

yrs <- unique(daily_data$Year)
for (yr in yrs) {
  cat("\nFitting for Year:", yr, "\n")
  
  subset_data <- daily_data %>%
    filter(Year == yr)
  
  if (nrow(subset_data) == 0) {
    cat("No rows for year", yr, "\n")
    next
  }
  
  rf_model <- ranger(
    formula    = Value ~ TV + swp + VPD,
    data       = subset_data,
    num.trees  = 100,
    importance = "impurity"
  )
  
  preds    <- predict(rf_model, subset_data)$predictions
  R2_value <- cor(subset_data$Value, preds)^2
  cat("R² for year", yr, "=", round(R2_value, 3), "\n")
  
  # SHAP
  shap_vals <- fastshap::explain(
    object = rf_model,
    X = subset_data %>% select(TV, swp, VPD),
    pred_wrapper = function(m, newdata) predict(m, data=newdata)$predictions,
    nsim         = 10,
    parallel     = FALSE,
    adjust       = TRUE
  )
  
  shap_df <- subset_data %>% select(Date) %>% bind_cols(shap_vals)
  
  shap_long <- shap_df %>%
    pivot_longer(
      cols = c("TV","swp","VPD"),
      names_to = "Predictor",
      values_to = "SHAP"
    )
  
  p <- ggplot(shap_long, aes(x=Date, y=SHAP, fill=Predictor)) +
    geom_bar(stat="identity", position="stack") +
    labs(
      title    = paste("Year", yr, "| SHAP over time"),
      subtitle = paste("R² =", round(R2_value, 3)),
      x        = "Date",
      y        = "SHAP Value",
      fill     = "Predictor"
    ) +
    scale_fill_brewer(palette="Set1") +
    theme_minimal() +
    theme(
      axis.text.x     = element_text(angle=45, hjust=1),
      legend.position = "top"
    )
  
  print(p)
}

dev.off()  # close Model2 PDF

cat("\n=== ALL DONE! Plots are saved to Model1_Individuals.pdf and Model2_Years.pdf ===\n")


# Load required libraries
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(ranger)
library(fastshap)
library(ggplot2)
library(RColorBrewer)
library(Metrics)

# Load data
path_to_csv <- "/home/peter/Downloads/wd_2015_ok.csv"
data <- read_csv(path_to_csv)

# Data preprocessing
data <- data %>%
  mutate(DateTime = as.POSIXct(paste(Date, Time), format='%m/%d/%y %H:%M'),
         Date = as.Date(DateTime)) %>%
  pivot_longer(cols = starts_with("JD"), names_to = "Individual", values_to = "Value") %>%
  drop_na(Date, Individual, TV, swp, VPD, Value) %>%
  group_by(Date, Individual) %>%
  summarise(across(c(TV, swp, VPD, Value), mean, na.rm = TRUE), .groups = 'drop') %>%
  mutate(Year = year(Date))

# Create empty dataframe to store results
results_df <- data.frame()

# Loop through Individuals and Years to extract metrics
individuals <- unique(data$Individual)
years <- unique(data$Year)

for (ind in individuals) {
  for (yr in years) {
    subset_data <- data %>% filter(Individual == ind, Year == yr)
    
    if (nrow(subset_data) < 10) next
    
    baseline <- mean(subset_data$Value, na.rm = TRUE)
    
    rf_model <- ranger(Value ~ TV + swp + VPD,
                       data = subset_data,
                       num.trees = 100,
                       importance = "impurity")
    
    predictions <- predict(rf_model, subset_data)$predictions
    
    # Metrics
    R2 <- cor(subset_data$Value, predictions)^2
    RMSE <- rmse(subset_data$Value, predictions)
    
    # AIC approximation using residual sum of squares
    n <- nrow(subset_data)
    RSS <- sum((subset_data$Value - predictions)^2)
    AIC_value <- n * log(RSS/n) + 2 * length(rf_model$variable.importance)
    
    # Append results
    results_df <- rbind(results_df, data.frame(
      Individual = ind,
      Year = yr,
      Baseline = baseline,
      R2 = R2,
      RMSE = RMSE,
      AIC = AIC_value
    ))
  }
}

# Save results to CSV
write_csv(results_df, "model_metrics11_results.csv")

# Print results
print(results_df)














# Load required libraries
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(ranger)
library(fastshap)
library(ggplot2)
library(RColorBrewer)
library(Metrics)

# Load data
path_to_csv <- "/home/peter/Downloads/wd_2015_ok.csv"
data <- read_csv(path_to_csv)

# Data preprocessing
data <- data %>%
  mutate(DateTime = as.POSIXct(paste(Date, Time), format='%m/%d/%y %H:%M'),
         Date = as.Date(DateTime)) %>%
  pivot_longer(cols = starts_with("JD"), names_to = "Individual", values_to = "Value") %>%
  drop_na(Date, Individual, TV, swp, VPD, Value) %>%
  group_by(Date, Individual) %>%
  summarise(across(c(TV, swp, VPD, Value), mean, na.rm = TRUE), .groups = 'drop') %>%
  mutate(Year = year(Date))

# Create empty dataframe to store results
results_df <- data.frame()

# Loop through Individuals and Years to extract metrics
individuals <- unique(data$Individual)
years <- unique(data$Year)

for (ind in individuals) {
  for (yr in years) {
    subset_data <- data %>% filter(Individual == ind, Year == yr)
    
    if (nrow(subset_data) < 10) next
    
    baseline <- mean(subset_data$Value, na.rm = TRUE)
    
    rf_model <- ranger(Value ~ TV + swp + VPD,
                       data = subset_data,
                       num.trees = 100,
                       importance = "impurity")
    
    predictions <- predict(rf_model, subset_data)$predictions
    
    # Metrics
    R2 <- cor(subset_data$Value, predictions)^2
    RMSE <- rmse(subset_data$Value, predictions)
    
    # AIC approximation using residual sum of squares
    n <- nrow(subset_data)
    RSS <- sum((subset_data$Value - predictions)^2)
    AIC_value <- n * log(RSS/n) + 2 * length(rf_model$variable.importance)
    
    # Append results
    results_df <- rbind(results_df, data.frame(
      Individual = ind,
      Year = yr,
      Baseline = baseline,
      R2 = R2,
      RMSE = RMSE,
      AIC = AIC_value
    ))
  }
}

# Save results to CSV
write_csv(results_df, "model_metrics_results22.csv")

# Print results
print(results_df)























































# Load required libraries
library(readr)
library(dplyr)
library(tidyr)
library(lubridate)
library(ranger)
library(Metrics)
library(broom)

# Load data
data <- read_csv("/home/peter/Downloads/wd_2015_ok.csv")

# Prepare data
data <- data %>%
  mutate(DateTime = as.POSIXct(paste(Date, Time), format = '%m/%d/%y %H:%M'),
         Date = as.Date(DateTime),
         Year = year(Date)) %>%
  drop_na(Date, JD1, JD2, JD3, JD4, JD5, TV, swp, VPD)

# Pivot JD1-JD5 to long format
data_long <- data %>%
  pivot_longer(cols = starts_with("JD"), names_to = "Individual", values_to = "Value") %>%
  drop_na(Value)

# Daily aggregation
daily_data <- data_long %>%
  group_by(Date, Year, Individual) %>%
  summarise(
    TV = mean(TV, na.rm = TRUE),
    swp = mean(swp, na.rm = TRUE),
    VPD = mean(VPD, na.rm = TRUE),
    Value = mean(Value, na.rm = TRUE),
    .groups = "drop"
  )

# Define metric extraction function
calculate_metrics <- function(actual, predicted, num_params) {
  n <- length(actual)
  residuals <- actual - predicted
  RSS <- sum(residuals^2)
  TSS <- sum((actual - mean(actual))^2)
  R2 <- 1 - RSS/TSS
  RMSE <- sqrt(mean(residuals^2))
  AIC <- n * log(RSS/n) + 2 * num_params
  list(R2 = R2, RMSE = RMSE, AIC = AIC, Baseline = mean(actual))
}

# Results storage
metrics_results <- data.frame()

### ======================================
### MODEL 1: Individuals across ALL YEARS
### ======================================
for(ind in unique(daily_data$Individual)){
  
  cat("\nModeling Individual:", ind, "across all years\n")
  
  subset_data <- daily_data %>% filter(Individual == ind)
  
  # Fit Random Forest
  rf_model <- ranger(Value ~ TV + swp + VPD, data = subset_data, importance = "impurity")
  
  preds <- predict(rf_model, subset_data)$predictions
  
  metrics <- calculate_metrics(subset_data$Value, preds, num_params = length(rf_model$variable.importance))
  
  metrics_results <- rbind(metrics_results, data.frame(
    Type = "Individual_AllYears",
    Group = ind,
    R2 = metrics$R2,
    RMSE = metrics$RMSE,
    AIC = metrics$AIC,
    Baseline = metrics$Baseline
  ))
}

### ======================================
### MODEL 2: Each YEAR across ALL INDIVIDUALS
### ======================================
for(yr in unique(daily_data$Year)){
  
  cat("\nModeling Year:", yr, "across all individuals\n")
  
  subset_data <- daily_data %>% filter(Year == yr)
  
  # Fit Random Forest
  rf_model <- ranger(Value ~ TV + swp + VPD, data = subset_data, importance = "impurity")
  
  preds <- predict(rf_model, subset_data)$predictions
  
  metrics <- calculate_metrics(subset_data$Value, preds, num_params = length(rf_model$variable.importance))
  
  metrics_results <- rbind(metrics_results, data.frame(
    Type = "Year_AllIndividuals",
    Group = yr,
    R2 = metrics$R2,
    RMSE = metrics$RMSE,
    AIC = metrics$AIC,
    Baseline = metrics$Baseline
  ))
}

# Save results to CSV
write_csv(metrics_results, "model_metrics_results.csv")

cat("\nMetrics extraction complete. Results saved to model_metrics_results.csv\n")

