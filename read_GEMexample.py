import os
import uproot
import numpy as np
import ROOT
import re
from tqdm import tqdm
import GEMClass as GEM

f = "reco/reco_run50898_3D.root"
run = GEM.get_numbers_from_filename(f)
try:
    events = uproot.open(f+":Events")
    events_GEM = uproot.open(f+":GEM_Events")
except Exception as e:
    print(f"Failed to open (maybe empty) {f}: {e}")
    exit()

main = ROOT.TFile("analGEM.root", "RECREATE")
for gem_name in ["GEM1", "GEM2", "GEM3"]:
    main.mkdir(gem_name)

allY = events_GEM["pmt_fullWaveform_Y"].array(library="np")
GEMid = events_GEM["pmt_wf_channel"].array(library="np")
integrals = {"GEM1": [], "GEM2": [], "GEM3": [],"GEM1_noise":[],"GEM1_signal":[]}

threshold = 100E-12  # Example threshold, adjust as needed

accumulated_waveforms = {"GEM1": None, "GEM2": None, "GEM3": None}
waveform_counts = {"GEM1": 0, "GEM2": 0, "GEM3": 0}

print("Processing waves...")
for i, y in tqdm(enumerate(allY)):
    gem_key = f"GEM{GEMid[i]-4}"  # Adjust based on GEMid mapping
    tempName = gem_key
    wave = GEM.GEMwave(y, f"{tempName} - run{run} - wave{i}")
    main.cd(tempName)

    # Obtain integrals
    chargeIntegral = wave.GetChargeIntegral()
    noiseIntegral = wave.NoiseIntegral()

    # Check if both integrals are below the threshold before appending
    if np.abs(chargeIntegral) <= threshold and np.abs(noiseIntegral) <= threshold:
        integrals[tempName].extend([chargeIntegral, noiseIntegral])
        if tempName=="GEM1":
            integrals["GEM1_noise"].append(noiseIntegral)
            integrals["GEM1_signal"].append(chargeIntegral)

    # Accumulate waveforms
    if accumulated_waveforms[gem_key] is None:
        accumulated_waveforms[gem_key] = y
    else:
        accumulated_waveforms[gem_key] += y
    waveform_counts[gem_key] += 1

    if i < 1000:
        wave.GetTGraph(write=True)

for gem_name, vals in integrals.items():
    main.cd()
    GEM.hist(vals, f"{gem_name} Integrals (C)", channels=1000)

average_gem_waves = {}

for gem_name, waveform_sum in accumulated_waveforms.items():
    if waveform_counts[gem_name] > 0:
        average_waveform = waveform_sum / waveform_counts[gem_name]
        average_gem_waves[gem_name] = GEM.GEMwave(average_waveform, f"{gem_name} Average")

for gem_name, gem_wave in average_gem_waves.items():
    graph = gem_wave.GetTGraph(write=True)

main.Close()