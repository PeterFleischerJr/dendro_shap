# dendro_shap

This repository contains code for processing dendrometer data and extracting model metrics using Random Forests. The main script [`metrics_calculation.R`](metrics_calculation.R) reads a CSV file named `wd_2015_ok.csv`, computes daily means, fits models for each individual and year, and writes a `model_metrics_results.csv` with R², RMSE, AIC, and baseline values.

To run the script:

```R
Rscript metrics_calculation.R
```

Ensure that the input CSV is located in the working directory.
