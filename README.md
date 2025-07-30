# dendroShap

This repository contains an R package for exploring SHAP values using dendrogram-based methods.

## Overview

The package provides helper functions to aggregate SHAP contributions and visualize them as dendrograms. This can help interpret complex tree-based models.

## Installation

```R
devtools::install_local("dendroShap")
```

## Basic Usage

```R
library(dendroShap)

# Example of loading SHAP values and building a dendrogram
shap_matrix <- matrix(rnorm(100), nrow = 10)
result <- build_dendrogram(shap_matrix)
plot(result)
```

