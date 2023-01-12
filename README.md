# Implicit Learning of a Kinematically Complex Multi-Articular Motor Skill

This repository contains the experiment anad analysis code for Solomon et al. (2023), a study investigating if a complex movement can be learned implicitly.

This study uses [TraceLab](https://github.com/LBRF/TraceLab).


## Experiment Code

The experiment code used to collect the data is found in the `Tracelab_Program` directory. The code is not functional as Tracelab has undergone several updates since the data was collected for the project. This directory is strictly for documentation purposes. For a functional version of the program, please visit [TraceLab](https://github.com/LBRF/TraceLab).


## Analysis Requirements

All dependencies for these scripts can be installed by running the following line:

```r
install.packages(c("tidyverse", "TSEntropies", "dtw","brms","tidybayes","emmeans","parametes","modelr"))
```

The scripts were run in R 4.1 for publication. They may work with older versions of R but are not guaranteed to function correctly.


## Analysis Code Usage

The raw data for the mproject can be found on the [Open Science Framework](https://osf.io/v45pq/).

1. Place the implicit knowledge test results from OSF (`*.csv`) in a created directory `Tracelab_Analysis/_Data/`.
2. Place all task data (`*.txt`) from OSF (`task.zip`) in a created directory `Tracelab_Analysis/_Data/task/`.
3. Place all the participant directories from OSF (`figure.zip`) in a created directory `Tracelab_Analysis/_Data/figure/`.
5. Open a new R session and set the working directory to the root `TraceLab_Analysis/` folder (or whatever you've renamed it to) using `setwd()` or the RStudio menu.
6. Run one of the following commands in the R terminal:

```r
source('./_Scripts/0_import.R') # imports task and figure data
source('./_Scripts/1_preprocessing.R') # imports and preprocesses data
source('./_Scripts/2_processing.R') # generates descriptives and runs statistical models
```

Running the processing script will run the preprocessing script and the import script, so in most cases you just want to run the third line.

The code is an exact copy used for the publication. Future commits will be made to optimize the code's performance.

