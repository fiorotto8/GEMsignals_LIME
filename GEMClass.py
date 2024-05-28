import os
import uproot
import numpy as np
import ROOT
import re
from tqdm import tqdm
from scipy.signal import butter, filtfilt, savgol_filter
from scipy.ndimage import gaussian_filter

def get_files_in_folder(folder_path):
    return [os.path.join(folder_path, f) for f in os.listdir(folder_path) if os.path.isfile(os.path.join(folder_path, f))]
def get_sc_integral(file, cuts=None):
    try:
        events = uproot.open(file+":Events")
        if cuts is None:
            cutIntegral = events.arrays(["sc_integral"], library="np")
        else:
            cutIntegral = events.arrays(["sc_integral"], cuts, library="np")
        return cutIntegral['sc_integral']
    except Exception as e:
        print(f"Failed to open (maybe empty) {file}: {e}")
        return []
def grapherr(x, y, ex, ey, x_string, y_string, name=None, color=4, markerstyle=22, markersize=2, write=True):
    plot = ROOT.TGraphErrors(len(x), np.array(x, dtype="d"), np.array(y, dtype="d"), np.array(ex, dtype="d"), np.array(ey, dtype="d"))
    plot_title = name if name else y_string + " vs " + x_string
    plot.SetNameTitle(plot_title, plot_title)
    plot.GetXaxis().SetTitle(x_string)
    plot.GetYaxis().SetTitle(y_string)
    plot.SetMarkerColor(color)
    plot.SetMarkerStyle(markerstyle)
    plot.SetMarkerSize(markersize)
    if write:
        plot.Write()
    return plot
def create_fill_TH2(name, x_name, y_name, z_name, x_vals, y_vals, weights, x_bins=15, y_bins=15, write=True):
    hist = ROOT.TH2F(name, name, x_bins, 0.99*np.min(x_vals), 1.01*np.max(x_vals), y_bins, 0.99*np.min(y_vals), 1.01*np.max(y_vals))
    for x, y, weight in zip(x_vals, y_vals, weights):
        hist.Fill(x, y, weight)
    hist.GetXaxis().SetTitle(x_name)
    hist.GetYaxis().SetTitle(y_name)
    hist.GetZaxis().SetTitle(z_name)
    if write:
        hist.Write()
    return hist
def get_numbers_from_filename(filename):
    return re.search(r'\d+', filename).group(0)
def fill_h(histo_name, array):
    for x in range (len(array)):
        histo_name.Fill((np.array(array[x] ,dtype="d")))
def hist(list, x_name, channels=100, linecolor=4, linewidth=4,write=True):
    array=np.array(list ,dtype="d")
    hist=ROOT.TH1D(x_name,x_name,channels,0.99*np.min(array),1.01*np.max(array))
    fill_h(hist,array)
    hist.SetLineColor(linecolor)
    hist.SetLineWidth(linewidth)
    hist.GetXaxis().SetTitle(x_name)
    hist.GetYaxis().SetTitle("Entries")
    if write==True: hist.Write()
    #hist.SetStats(False)
    hist.GetYaxis().SetMaxDigits(3);
    hist.GetXaxis().SetMaxDigits(3);
    return hist
def graph(x,y,x_string, y_string,name=None, color=4, markerstyle=22, markersize=2,write=True):
        plot = ROOT.TGraph(len(x),  np.array(x  ,dtype="d")  ,   np.array(y  ,dtype="d") )
        if name is None: plot.SetNameTitle(y_string+" vs "+x_string,y_string+" vs "+x_string)
        else: plot.SetNameTitle(name, name)
        plot.GetXaxis().SetTitle(x_string)
        plot.GetYaxis().SetTitle(y_string)
        plot.SetMarkerColor(color)#blue
        plot.SetMarkerStyle(markerstyle)
        plot.SetMarkerSize(markersize)
        if write==True: plot.Write()
        return plot
def nparr(list):
    return np.array(list, dtype="d")

# You can use this function to directly create all the GEMwave objects from the root files if needed!
def process_gem_data(root_file):
    gem_waves_obj = { "GEM1": [], "GEM2": [], "GEM3": [] }
    for gem_id, waveforms in root_file.gem_data.items():
        for i, y in enumerate(waveforms):
            gem_wave = GEMwave(y, f"{gem_id}_wave_{i}")
            gem_waves_obj[gem_id].append(gem_wave)
            #print(f"Processed {gem_wave.name} with {gem_wave.samples} samples")
    return gem_waves_obj

class ROOTFile:
    def __init__(self, path):
        """
        Constructor to initialize the ROOTFile object.

        Parameters:
        path (str): Path to the ROOT file.

        It has one important attribute which is a dictionary:
        - self.gem_data: A dictionary with keys "GEM1", "GEM2", "GEM3" containing the sampled voltages arrays
                        for each respective GEM ID.
        """
        # Store the path to the ROOT file
        self.path = path
        # Extract the run number from the file path using regular expressions
        self.run_num = re.search(r'\d+', path).group(0)
        try:
            # Open the ROOT file and access the GEM events tree
            events_GEM = uproot.open(path + ":GEM_Events")
        except Exception as e:
            # Print an error message and exit if the file cannot be opened
            print(f"Failed to open (maybe empty) {path}: {e}")
            exit()
        # Extract the sampled voltage data (waveforms) from the GEM events tree
        allY = events_GEM["pmt_fullWaveform_Y"].array(library="np")
        # Extract the GEM ID data from the GEM events tree
        GEMid_temp = events_GEM["pmt_wf_channel"].array(library="np")
        # Adjust GEM IDs by subtracting 4
        adjusted_GEMID = GEMid_temp - 4
        # Define the valid GEM IDs
        valid_gem_ids = [1, 2, 3]
        # Check if all elements in adjusted_GEMID are either 1, 2, or 3
        invalid_gem_ids = ~np.isin(adjusted_GEMID, valid_gem_ids)
        # If there are any invalid GEM IDs, print an error message and exit
        if np.any(invalid_gem_ids):
            print("Something wrong with the GEM number and the relative ADC channel")
            exit()
        # Create a dictionary to store the sampled voltage arrays for each GEM ID
        self.gem_data = {"GEM1": [], "GEM2": [], "GEM3": []}
        # Populate the dictionary with the corresponding allY arrays
        for i, gem_id in enumerate(adjusted_GEMID):
            if gem_id == 1:
                self.gem_data["GEM1"].append(allY[i])
            elif gem_id == 2:
                self.gem_data["GEM2"].append(allY[i])
            elif gem_id == 3:
                self.gem_data["GEM3"].append(allY[i])
        # Convert lists to numpy arrays for efficient processing
        for key in self.gem_data:
            self.gem_data[key] = np.array(self.gem_data[key])

class GEMwave:
    def __init__(self, y, name, invert=False, dynRange=4096, ampGain=10):
        """
        Constructor to initialize the GEMwave object.

        Parameters:
        y (array): The sampled voltage array.
        name (str): Name of the waveform.
        invert (bool): Whether to invert the waveform. Default is False.
        dynRange (int): Dynamic range for normalization. Default is 4096.
        ampGain (int): Amplification gain for normalization. Default is 10.
        """
        self.name = name  # Name of the waveform.
        # Invert the waveform if requested. Normalize by amplification gain and dynamic range.
        if not invert:
            self.y = (np.array(y) / ampGain) / dynRange
        else:
            self.y = -1 * (np.array(y) / ampGain) / dynRange
        self.scopeImpedence = 50  # Ohm, impedance of the oscilloscope used.

        """ I'm really not sure about that....
        # Adjust voltage for 50 ohm input impedance of the oscilloscope
        self.y = self.y * (50 + self.scopeImpedence) / 50
        """

        self.samples = len(self.y)  # Total number of samples in the waveform.
        self.dt = self.dtdefinition()  # Determine the time interval between samples.
        # Create an x-axis time vector from 0 to the total time span of the samples.
        self.x = np.linspace(0, (self.samples - 1) * self.dt, self.samples)

#### A BUNCH OF FILTERS
    def ApplyBandPassFilter(self, lowcut=100E6, highcut=300E6, order=10):
        # Convert the cutoff frequencies to normalized form
        nyquist = 0.5 / self.dt  # Nyquist frequency
        low = lowcut / nyquist
        high = highcut / nyquist
        # Create bandpass Butterworth filter
        b, a = butter(order, [low, high], btype='band')
        # Apply the filter
        self.y = filtfilt(b, a, self.y)

    def ApplyLowPassFilter(self, cutoff=2E8, order=5):
        # Calculate Nyquist frequency
        nyquist = 0.5 / self.dt
        # Check if the cutoff frequency is too high
        if cutoff >= nyquist:
            raise ValueError("Cutoff frequency must be less than half the sampling rate (Nyquist frequency).")
        # Calculate normalized cutoff frequency
        normal_cutoff = cutoff / nyquist
        # Generate low-pass Butterworth filter coefficients
        from scipy.signal import butter, filtfilt
        b, a = butter(order, normal_cutoff, btype='low', analog=False)
        # Apply the filter to self.y using filtfilt to avoid phase shift
        self.y = filtfilt(b, a, self.y)

    def ApplyGaussianFilter(self, sigma=5):
        # Apply Gaussian filter with specified sigma
        self.y = gaussian_filter(self.y, sigma=sigma)

    def ApplyMovingAverage(self, window_size=20):
        # Apply moving average filter with specified window size
        if window_size % 2 == 0:
            window_size += 1  # Ensure window size is odd
        window = np.ones(window_size) / window_size
        self.y = np.convolve(self.y, window, 'same')

    def ApplySavGolFilter(self, window_length=51, polyorder=5):
        # Apply the Savitzky-Golay filter to self.y and update it.
        # Ensure window_length is odd and less than the number of samples.
        if window_length % 2 == 0:
            window_length += 1  # Make sure the window length is odd
        if window_length > self.samples:
            print("Window length must be less than number of samples.")
            return
        self.y = savgol_filter(self.y, window_length, polyorder)

#### OTHER METHODS
    def GetChargeIntegral(self, range=[0.25E-6, 0.5E-6]):
        # Calculate the integral of the charge within a specified time range.
        sel_y = self.y[(self.x > range[0]) & (self.x < range[1])]
        return (np.sum(sel_y) * self.dt) / self.scopeImpedence

    def GetTGraph(self, write=False):
        # Generate a graph of the waveform.
        return graph(self.x, self.y, "time(s)", "voltage(V)", self.name, write=write)

    def dtdefinition(self):
        # Define the sampling interval based on the number of samples.
        if self.samples == 1024:
            dt = 1.33E-9  # seconds, specific for 1024 samples.
        elif self.samples == 4000:
            dt = 4E-9  # seconds, specific for 4000 samples.
        else:
            print("Unknown sampling --> dt set to 0!")
            dt = 0
        return dt

    def NoiseIntegral(self, range=[0, 0.2E-6]):
        # Calculate the integral of the noise within a specified time range.
        sel_y = self.y[(self.x > range[0]) & (self.x < range[1])]
        return (np.sum(sel_y) * self.dt) / self.scopeImpedence

    def GetFFT(self):
        # Compute the FFT of self.y
        fft_result = np.fft.fft(self.y)
        # Compute the corresponding frequencies
        freq = np.fft.fftfreq(self.samples, d=self.dt)
        # Only take the half of the spectrum corresponding to positive frequencies
        positive_freq_indices = freq > 0
        freq = freq[positive_freq_indices]
        fft_result = fft_result[positive_freq_indices]
        # Calculate the power spectrum (magnitude squared)
        power = np.abs(fft_result) ** 2
        # Return the frequency and power
        return freq, power

    def ApplyIntegratorWithDecay(self, feedback_capacitance=1E-6, feedback_resistor=100):
        """
        Simulates a charge amplifier/integrator with a decay determined by the feedback components.
        This version outputs the integrated charge directly.
        No gain implemented to get the full charge
        Parameters:
        - feedback_capacitance (float): Feedback capacitance in farads.
        - feedback_resistor (float): Feedback resistor in ohms.
        """
        tau = feedback_resistor * feedback_capacitance  # Time constant
        decay_factor = np.exp(-self.dt / tau)  # Calculate decay factor for each sample
        accumulated_charge = np.zeros_like(self.y)  # Initialize accumulated charge array
        # Perform integration with decay to accumulate charge
        for i in range(1, len(self.y)):
            accumulated_charge[i] = accumulated_charge[i - 1] * decay_factor + self.y[i] * self.dt

        return accumulated_charge