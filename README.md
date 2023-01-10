# Implicit Learning of a Kinematically Complex Multi-Articular Motor Skill

This repository contains the experiment code for Solomon et al. (2023), a study investigating if a complex movement can be learned implicitly.

This study uses [TraceLab](https://github.com/LBRF/TraceLab).

## Requirements

All dependencies for these scripts can be installed by running the following line:

```r
install.packages(c("tidyverse", "TSEntropies", "dtw","brms","tidybayes","emmeans","parametes","modelr"))
```

The scripts were run in R 4.1 for publication. They may work with older versions of R but are not guaranteed to function correctly.


## Usage

The raw data for the mproject can be found on the [Open Science Framework](https://osf.io/v45pq/).

1. Place all task data from OSF (`*.txt`)  in `_Data/task/`.
2. Place all figure data from OSF (`*.zip`)  in `_Data/figure/`.
2. Place the implicit knowledge test results from OSF (`*.csv`) in `_Data/`.
3. Open a new R session and set the working directory to the root `TraceLabAnalysis/` folder (or whatever you've renamed it to) using `setwd()` or the RStudio menu.
4. Run one of the following commands in the R terminal:

```r
source('./_Scripts/0_import.R') # imports task and figure data
source('./_Scripts/1_preprocessing.R') # imports and preprocesses data
source('./_Scripts/2_processing.R') # generates descriptives and runs statistical models
```

Running the processing script will run the preprocessing script and the import script, so in most cases you just want to run the thira line.

