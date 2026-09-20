# scarf-QTL

Code for the manuscript:

**scarf-QTL: Incorporating cell states for mapping eQTLs in single-cell studies**

## Overview

This repository contains code for scarf-QTL, a statistical framework for genome-wide mapping of single-cell eQTLs while incorporating cell states. The method combines a functional mixed model with retrospective association tests to identify both static and cell-state-dependent genetic effects on gene expression.

The repository includes scripts for real-data analysis, simulation studies, and a minimal reproducible example. Large raw data files and intermediate results are not included in this repository.

## Repository structure

The repository is organized into real-data and simulation pipelines.

- `code/scarf_QTL_function.R`: core functions for scarf-QTL and the pseudobulk baseline. The scarf-QTL implementation includes model fitting, association testing, and effect estimation.
- `code/realdata/`: scripts for reproducing the real-data analyses in the manuscript.
- `code/simulation/`: scripts for reproducing the simulation studies in the manuscript.
- `code/run_example.R`: a minimal reproducible example illustrating the basic scarf-QTL workflow on a synthetic dataset.
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

## Example

A minimal reproducible example is provided in `code/run_example.R`.

The example uses a fully synthetic single-cell eQTL dataset with 1,000 individuals, five genes, and 3,000 synthetic cis-SNPs per gene. A cell-state-dependent genetic effect is introduced for one gene-SNP pair to demonstrate the retrospective association tests. The synthetic example data can be generated directly using `code/run_example.R` and are also provided as `data/toy_example_data.RData` for convenience.

From the root directory of the repository, run:

```r
source("code/run_example.R")
