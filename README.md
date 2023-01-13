# Implicit Learning of a Kinematically Complex Multi-Articular Motor Skill Analysis

This repository contains the analysis code for Solomon et al. (2023), a study investigating if a complex movement can be learned implicitly.

This study uses [TraceLab](https://github.com/LBRF/TraceLab).


## Experiment Code

The experiment code used to collect the data is found in the [Implicit_CME](https://github.com/LBRF-Projects/Implicit_CME_Solomon2023_Experiment) repository.

## Analysis Requirements

All dependencies for these scripts can be installed by running the following line:

```r
install.packages(c("tidyverse", "TSEntropies", "dtw","brms","tidybayes","emmeans","parametes","modelr"))
```

The scripts were run in R 4.1 for publication. They may work with older versions of R but are not guaranteed to function correctly.


## Analysis Code Usage

The raw data for the project can be found on the [Open Science Framework](https://osf.io/v45pq/).

1. Place the implicit knowledge test results from OSF (`*.csv`) in a created directory `Tracelab_Analysis/_Data/`.
2. Extract the (`task.zip`) from OSF in the created directory `Tracelab_Analysis/_Data/`.
3. Extract the (`figure.zip`) from OSF in the created directory `Tracelab_Analysis/_Data/`.
4. Open a new R session and set the working directory to the root `TraceLab_Analysis/` folder (or whatever you've renamed it to) using `setwd()` or the RStudio menu.
5. Run one of the following commands in the R terminal:

```r
source('./_Scripts/0_import.R') # imports task and figure data
source('./_Scripts/1_preprocessing.R') # imports and preprocesses data
source('./_Scripts/2_processing.R') # generates descriptives and runs statistical models
```

Running the processing script will run the preprocessing and the import scripts, so in most cases you just want to run the third line.

The code is an exact copy used for the publication. Future commits will be made to optimize the code's performance.

