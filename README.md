# Astrocyte-dependent spatial coding analysis

MATLAB code used in *“Selective roles of astrocytes in the formation and stabilization of new hippocampal place fields”* for one‑photon CA1 calcium imaging analysis.

## Code overview

`bayesian_decoder.m` : Bayesian decoding of the animal’s position from CA1 population calcium activity (single session).

`GetSpaInfo.m` : Computes Skaggs spatial information (bits/event) from occupancy maps and smoothed rate maps for individual neurons.

`pv_correlation_mouseMedian.m` : Computes mouse‑level population vector (PV) correlations between two sessions using CellReg‑matched neurons.

## Demo data overview

The repository includes a sample file `information_data.mat`, which is used when you run `bayesian_decoder` from the main folder. This demo file contains:

- `SmoothMat_Total` – spatial rate maps for each CA1 neuron  
- `spkMat_Total` – spike/calcium event timestamps for each neuron  
- `posTable_array` – time‑stamped animal positions  
- `calmea` – mean event rates for all neurons  

These rate maps can also be used to test `GetSpaInfo.m` by providing an occupancy map (`OCCMat`) and a smoothed rate map (`SkaggsrateMat`) derived from `SmoothMat_Total`.

## System requirements

- Software  
  - MATLAB R2022b  
  - Statistics and Machine Learning Toolbox  
- Operating systems  
  - Developed and tested on 64‑bit Microsoft Windows 11 Education (Version 10.0, build 26100)  
- Hardware  
  - Minimum 8 GB RAM; ≥16 GB RAM recommended  
  - No non‑standard hardware required (miniscope only needed for data acquisition)

## Installation and demo

1. Download this repository as a ZIP from GitHub and unzip it to a local folder.  
2. In MATLAB, set the unzipped folder as your current folder, or add it to the path:
3. Run the demo:bayesian_decoder *Run this command from the repository’s main folder to execute the position‑decoding demo on the bundled sample dataset; on a standard desktop, it typically finishes within about 1 hour.
