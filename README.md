# BP-CNNI.2025

This repository contains the analysis code used in our study investigating how brain network dynamics, derived from resting-state fMRI, change before and after interventions in patients with depression and healthy controls.  
We employed a Hidden Markov Model (HMM) approach to extract brain dynamics features and conducted group-level comparisons.

## Repository structure

- `data/` — [User-provided] Preprocessed timeseries data (rows = time points, columns = ROIs)
- `scripts/` — Analysis scripts (HMM, feature extraction, calc metastates)
- `results/` — outputs of HMM analysis
- `features/` — HMM features
- `README.md` — This document
- `LICENSE` — License information

## Requirements

- **MATLAB**  
  Required toolbox: [HMM-MAR toolbox](https://github.com/OHBA-analysis/HMM-MAR) 

- **R**  
  Required R packages (specified in each script)

## Usage

1. Place your preprocessed fMRI timeseries data into the `data/` folder.  
   The expected format is subject-wise timeseries `.csv` files.

2. Run the analysis scripts in order:  
   `scripts/01_run_HMM_MAR.m`  
   `scripts/02_extract_HMM_features.m`  
   `scripts/03_calc_metastates.R`  
