# Astrocyte-dependent spatial coding analysis

This repository contains MATLAB code used in the study *“Selective roles of astrocytes in the formation and stabilization of new hippocampal place fields”*. The scripts analyze one‑photon CA1 calcium imaging data to quantify how astrocytic Gq activation affects the precision and stability of hippocampal spatial maps.

## Code overview

`bayesian_decoder.m` : Performs Bayesian decoding of the animal’s position from ensemble calcium activity within a single session. It takes precomputed spatial rate maps and calcium event trains, applies a naïve Bayes decoder, and returns decoded positions, decoding error, and related summary statistics.

`bayesian_decoder_demo.m` : Lightweight demo script that calls `bayesian_decoder.m` with reduced `num_iterations` and `num_shuffles` so that users and reviewers can quickly reproduce decoding results on example data. Only the data path (and optional parameters) need to be edited to run this demo on a new dataset.

`pv_correlation_mouseMedian.m` : Standalone script that computes mouse‑level population vector (PV) correlation between two sessions. The script (1) loads CellReg results to identify neurons tracked across sessions, (2) selects neurons that are active in both sessions (e.g. mean event rate ≥ 0.01 Hz), (3) builds population rate matrices by stacking each neuron’s spatial rate map into a vector, and (4) for each spatial bin, computes the Pearson correlation between the population activity vectors of session 1 and session 2. Bins outside the arena or with no activity are excluded, and the distribution of bin‑wise PV correlations is summarized at the mouse level using the median and interquartile range to provide a robust, non‑parametric measure of spatial map stability.

## Demo

For a quick demonstration:

- Run `bayesian_decoder_demo.m` to see how Bayesian decoding is applied to a single session with computationally light settings (small numbers of iterations and shuffles).
- Run `pv_correlation_mouseMedian.m`, select the experiment folder when prompted, and set the two session indices in the script; the code will output the mouse‑level median PV correlation and save all bin‑wise correlations to a `.mat` file.

The full analysis code is available at:

GitHub: https://github.com/rkcldk/YSL-lab-Yoo-et-al
