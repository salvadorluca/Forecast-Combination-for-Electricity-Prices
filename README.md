# Forecast Combinations for Electricity Prices: The Case of the Italian Market

This repository contains the R scripts and datasets used in my Bachelor's thesis, which explores forecast combinations for electricity prices in the Italian market. The thesis evaluates the effectiveness of various forecasting models and their combinations to enhance the accuracy of electricity price predictions.

## Abstract
This Bachelor's thesis delves into the utility of forecast combinations to enhance the accuracy of electricity price predictions in the Italian market. A variety of forecasting models are employed and assessed for their individual and combined effectiveness. The study features time series forecasting models such as ARIMA, Exponential Smoothing, Random Forest, and Spline models. These models are individually evaluated for their predictive performance and then combined using various weighting methods. The objective is to leverage the strengths of multiple forecasts through model averaging, a technique that helps offset the weaknesses inherent in individual models by blending their forecasts. This approach aims to produce more reliable and accurate predictions for electricity prices, facilitating better decision-making in energy market operations.

---

## Repository Structure

- **/Additional files/**:
  - Contains all auxiliary files necessary for running the main script.
  
- **Thesis.R**:
  - This R script includes all the data analysis, model implementation, and forecasting combinations used in the thesis. It is the primary script for replicating the study's results.

- **Thesis.pdf**:
  - A PDF version of the complete Bachelor’s thesis, which discusses the methodologies, findings, and implications of the research.

- **dati_pun.xlsx**:
  - The Excel spreadsheet containing the dataset of electricity prices used in the thesis.

## Setup and Running

### Prerequisites

Ensure you have R installed on your machine. Additionally, the following R libraries are required to run the scripts:

- `aTSA`, `astsa`, `car`, `chron`, `clock`, `cowplot`, `datasets`, `DEoptim`, `dplyr`, `fable`, `fBasics`, `fImport`, `forecast`, `ForecastComb`, `fpp`, `fpp2`, `fpp3`, `gam`, `GenSA`, `ggplot2`, `GGally`, `glmnet`, `gridExtra`, `highfrequency`, `kableExtra`, `knitr`, `lavaan`, `lme4`, `lubridate`, `MASS`, `mice`, `Metrics`, `mvtnorm`, `plotly`, `PortfolioAnalytics`, `quantmod`, `randomForest`, `readxl`, `reshape`, `reshape2`, `robust`, `ROI`, `ROI.plugin.glpk`, `ROI.plugin.quadprog`, `sandwich`, `sn`, `stats`, `timeSeries`, `tibble`, `tidyr`, `tidyverse`, `tseries`, `tsibble`, `xtable`, `xts`, `zoo`

### Running the Scripts

1. Clone this repository to your local machine.
2. Open the scripts in RStudio or your preferred R environment.
3. Set the working directory to the folder containing the scripts.
4. Run the scripts individually to perform the analysis as described in the thesis.

### Next Steps
- **Code Review and Optimization**: Review and refactor the code for efficiency and readability.
- **Libraries**: Check whether all libraries are useful
