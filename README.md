# scarf-QTL

Code for the manuscript:

**A retrospective association test for static and dynamic single-cell eQTLs**

## Overview

This repository contains code for scarf-QTL, a statistical framework for genome-wide mapping of static and dynamic single-cell eQTLs. The method combines a functional mixed model with a retrospective association test to detect genetic effects that vary along pseudotime.

The repository includes scripts for real-data analysis and simulation studies. Large raw data files and intermediate results are not included in this repository.

## Repository structure

```text
code/
  realdata/       real-data analysis scripts
  simulation/     simulation scripts

data/
  input data placeholders or small example files only

results/
  realdata/       summary outputs and figures for real-data analysis
  simulation/     summary outputs and figures for simulations
