# scarf-QTL

Code for the manuscript:

**A retrospective association test for static and dynamic single-cell eQTLs**

## Overview

This repository contains code for scarf-QTL, a statistical framework for genome-wide mapping of static and dynamic single-cell eQTLs. The method combines a functional mixed model with a retrospective association test to detect genetic effects that vary along pseudotime.

The repository includes scripts for real-data analysis and simulation studies. Large raw data files and intermediate results are not included in this repository.

## Repository structure

The repository is organized into real-data and simulation pipelines.

- `code/scarf_QTL_function.R`: core functions for scarf-QTL and the pseudobulk baseline. The scarf-QTL implementation includes model fitting, association testing, and effect estimation.
- `code/realdata/`: scripts for reproducing the real-data analyses in the manuscript.
- `code/simulation/`: scripts for reproducing the simulation studies in the manuscript.
- `data/`: input data directory.
- `results/`: output directory for intermediate and summary results.

### Real-data pipeline
The real-data analysis scripts should be run in the following order:

- `1_prepare_data.R`
- `2_fit_null_hypothesis.R`
- `3_retrospective_test.R`
- `4_estimation.R`
- `5_cluster_and_summary.R`
- `6_plot_realdata.R`

### Simulation pipeline
The simulation scripts should be run in the following order:

- `1_generate_data.R`
- `2_test.R`
- `3_permutation.R`
- `4_runtime.R`
- `5_summary.R`
- `6_plot_simulation.R`

Each script contains detailed comments on the required inputs, outputs, and usage.

P.S. A lightweight toy example will be added in a future update to illustrate the workflow on a small dataset.
