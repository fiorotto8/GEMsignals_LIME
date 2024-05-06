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
    main.mkdir(f"SIGNAL/{gem_name}")
    main.mkdir(f"FFT/{gem_name}")
    main.mkdir(f"INTEGRAL/{gem_name}")

allY = events_GEM["pmt_fullWaveform_Y"].array(library="np")
GEMid = events_GEM["pmt_wf_channel"].array(library="np")
integrals = {"GEM1": [], "GEM2": [], "GEM3": [],"GEM1_noise":[],"GEM1_signal":[]}

threshold = 100E-12  # Example threshold, adjust as needed

accumulated_waveforms = {"GEM1": None, "GEM2": None, "GEM3": None}
waveform_counts = {"GEM1": 0, "GEM2": 0, "GEM3": 0}

fft_power_sums = {"GEM1": None, "GEM2": None, "GEM3": None}
fft_counts = {"GEM1": 0, "GEM2": 0, "GEM3": 0}

int_sums = {gem_key: None for gem_key in ["GEM1", "GEM2", "GEM3"]}
int_counts = {gem_key: 0 for gem_key in ["GEM1", "GEM2", "GEM3"]}

print("Processing waves...")
for i, y in tqdm(enumerate(allY)):
    gem_key = f"GEM{GEMid[i]-4}"  # Adjust based on GEMid mapping
    tempName = gem_key
    wave = GEM.GEMwave(y, f"{tempName} - run{run} - wave{i}")
    main.cd(f"SIGNAL/{tempName}")

    # Obtain integrals
    #wave.ApplyBandPassFilter()
    #wave.ApplyLowPassFilter()
    chargeIntegral = wave.GetChargeIntegral()
    noiseIntegral = wave.NoiseIntegral()
    if i < 100: wave.GetTGraph(write=True)

#### SIGNALS
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

### FFTs
    # Compute FFT and accumulate
    freq, power = wave.GetFFT()
    if fft_power_sums[gem_key] is None:
        fft_power_sums[gem_key] = power
    else:
        fft_power_sums[gem_key] += power
    fft_counts[gem_key] += 1

    main.cd(f"FFT/{tempName}")
    if i< 100: GEM.graph(freq,power,"Freq (Hz)", "Power",f"FFT_wave_{i}")

#### INTEGRATE
    x, y_int = wave.x, wave.ApplyIntegratorWithDecay()

    if int_sums[gem_key] is None:
        int_sums[gem_key] = y_int.copy()  # Make a copy to avoid reference issues
    else:
        int_sums[gem_key] += y_int
    int_counts[gem_key] += 1

    main.cd(f"INTEGRAL/{tempName}")
    if i < 100: GEM.graph(x, y_int, "time (s)", "Charge (C)", f"int_wave_{i}")

############################################################################################

### SIGNALS
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

### FFTs
average_fft = {}
for gem_name, power_sum in fft_power_sums.items():
    if fft_counts[gem_name] > 0:
        average_power = power_sum / fft_counts[gem_name]
        average_fft[gem_name] = (freq, average_power)  # Store frequency and average power
for gem_name, (freq, average_power) in average_fft.items():
    # Using your custom graph function to create and write the plot
    # This assumes you have a color, marker style, etc., that you prefer
    graph = GEM.graph(
        freq, average_power,"Frequency (Hz)", "Average Power", name=f"{gem_name} Average FFT")

#### INTEGRATE
average_int = {}
for gem_name, y_sum in int_sums.items():
    if int_counts[gem_name] > 0:
        average_y = y_sum / int_counts[gem_name]
        average_int[gem_name] = (x, average_y)  # Store time and average integrated output

for gem_name, (x, average_y) in average_int.items():
    # Using your custom graph function to create and write the plot
    GEM.graph(
        x, average_y, "time (s)", "Average Integrated Charge (C)", name=f"{gem_name} Average Integration")

main.Close()