# scarf-QTL

Code for the manuscript:

**A retrospective association test for static and dynamic single-cell eQTLs**

## Overview

This repository contains code for scarf-QTL, a statistical framework for genome-wide mapping of static and dynamic single-cell eQTLs. The method combines a functional mixed model with a retrospective association test to detect genetic effects that vary along pseudotime.

The repository includes scripts for real-data analysis, simulation studies, and a lightweight toy example. Large raw data files and intermediate results are not included in this repository.

## Repository structure

The repository is organized into real-data and simulation pipelines.

- `code/scarf_QTL_function.R`: core functions for scarf-QTL and the pseudobulk baseline. The scarf-QTL implementation includes model fitting, association testing, and effect estimation.
- `code/realdata/`: scripts for reproducing the real-data analyses in the manuscript.
- `code/simulation/`: scripts for reproducing the simulation studies in the manuscript.
- `example/`: a lightweight toy example illustrating the basic scarf-QTL workflow on a small dataset.
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

## Toy example

A lightweight toy example is provided in `example/` to illustrate the basic usage of scarf-QTL on a small dataset.

- `example/produce_toy_example_data.R`: generates the toy dataset and saves it as `example/toy_example_data.RData`
- `example/run_toy_example.R`: loads the toy data and runs scarf-QTL on the example dataset

This toy example is intended to demonstrate the basic usage of scarf-QTL, rather than reproduce the full manuscript analyses.
