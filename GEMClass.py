import os
import uproot
import numpy as np
import ROOT
import re
from tqdm import tqdm

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

class GEMwave:
    def __init__(self, y,name,invert=False,dynRange=4096,ampGain=10):
        self.name=name
        if invert==False: self.y=(nparr(y)/ampGain)/dynRange
        else:self.y=-1*(nparr(y)/ampGain)/dynRange
        self.scopeImpedence = 50#ohm
        self.samples=len(self.y)
        self.dt=self.dtdefinition()
        self.x=np.linspace(0,(self.samples-1)*self.dt,self.samples)

    def GetChargeIntegral(self,range=[0.25E-6,0.5E-6]):
        sel_y = self.y[(self.x > range[0]) & (self.x < range[1])]
        return (np.sum(sel_y)*self.dt)/self.scopeImpedence

    def GetTGraph(self,write=False):
        return graph(self.x,self.y,"time(s)","voltage(V)",self.name,write=write)

    def dtdefinition(self):
        if self.samples==1024: dt=1.33E-9#seconds
        elif self.samples==4000: dt=4E-9#seconds
        else:
            print("Unkown sampling --> dt set to 0!")
            dt=0
        return dt

    def NoiseIntegral(self,range=[0,0.2E-6]):
        sel_y = self.y[(self.x > range[0]) & (self.x < range[1])]
        return (np.sum(sel_y)*self.dt)/self.scopeImpedence
