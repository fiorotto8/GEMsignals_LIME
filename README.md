# GEM Wave Analysis Toolkit

## Overview

This toolkit provides a set of tools for analyzing Gas Electron Multiplier (GEM) detector waveforms, including the calculation of charge integrals, noise analysis, and visualization of waveform data. It is designed to work with ROOT and uproot for efficient data handling and analysis.

## Dependencies

- Python 3.x
- ROOT (with PyROOT enabled)
- uproot
- numpy
- tqdm

## Usage

- Reading and Analyzing Waveforms: Use the GEMwave class to read individual GEM detector waveforms from data files, perform signal and noise integration, and generate plots for analysis.
- Histogram and Graph Generation: The toolkit includes functions to create histograms and graphs of the analyzed data, facilitating the visual inspection of signal characteristics and noise levels.
- Batch Processing: An example script demonstrates how to process multiple waveforms, generate summary histograms, and write the results to a ROOT file for further analysis.
