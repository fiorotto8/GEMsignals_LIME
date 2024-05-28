# Import necessary libraries
import os
import uproot
import numpy as np
import ROOT
import re
from tqdm import tqdm
import GEMClass as GEM

f = "reco/reco_run50898_3D.root"
rootFile=GEM.ROOTFile(f)
#! you can use this function to instantiate all the GEM objects
GEMwaves=GEM.process_gem_data(rootFile)
# Now you can access the GEMwave objects as follows:
# GEMwaves["GEM1"], GEMwaves["GEM2"], GEMwaves["GEM3"]
#! Or you may cycle over them with this:
"""
for gem_id, waveforms in root_file.gem_data.items():
    for i, y in enumerate(waveforms):
        gem_wave = GEMwave(y, f"{gem_id}_wave_{i}")
"""

# Create a new ROOT file for output
main = ROOT.TFile("analGEM_test.root", "RECREATE")
# Iterate over the processed_gem_waves and call GetTGraph every 100 waveforms
for gem_id, gem_wave_list in GEMwaves.items():
    for i, gem_wave in enumerate(gem_wave_list):
        if i % 250 == 0:
            tgraph = gem_wave.GetTGraph(write=True)
            # Print or use the TGraph as needed
            print(f"Generated TGraph for {gem_wave.name}: {tgraph}")