import ROOT
from ROOT import TFile
from array import array
from ROOT import TChain
from ROOT import gStyle
from ROOT import TStyle

import numpy as np
import matplotlib.pyplot as plt
#import mplhep as hep
import os

# Enable multi-threading
ROOT.ROOT.EnableImplicitMT()
 
# Include necessary header
muonSV_histograms_header_path = "nano_plot_data_cms_paper_python.h"
 
ROOT.gInterpreter.Declare('#include "{}"'.format(muonSV_histograms_header_path))

#verbosity = ROOT.Experimental.RLogScopedVerbosity(ROOT.Detail.RDF.RDFLogChannel(), ROOT.Experimental.ELogLevel.kInfo)

data = TChain();

'''
#for k in range(1, 4740):     #0.49/fb of era B
for k in range(1, 395):
    str_input_file = "/vols/cms/mc3909/bparkProductionAll_V1p0/ParkingBPH1_Run2018B-05May2019-v2_MINIAOD_v1p0_generationSync"+"/output_"+ str(k)+".root"
    if os.path.isfile(str_input_file) == True and k != 9 and k != 20 and k != 21 and k != 39 and k != 44 and k != 54 and k != 63 and k != 109 and k != 116 and k != 145 and k != 174 and k != 185 and k != 186 and k != 195 and k != 200 and k != 201 and k != 215 and k != 220 and k != 225 and k != 227 and k != 228 and k != 229 and k != 257 and k != 262 and k != 304 and k != 308 and k != 325 and k != 359 and k != 377 and k != 388 and k != 391 :
        data.Add(str_input_file+"/Events"); # Add the root file to the TChain chain
        print("adding ", str_input_file) 

       
#for k in range(1, 4569):     #0.47/fb of era C
for k in range(1, 318):      
    str_input_file_C = "/vols/cms/mc3909/bparkProductionAll_V1p0/ParkingBPH1_Run2018C-05May2019-v1_MINIAOD_v1p0_generationSync"+"/output_"+ str(k)+".root"
    if os.path.isfile(str_input_file_C) == True and k != 145 and k != 157 and k != 189 and k != 203 and k != 235 and k != 275:
        data.Add(str_input_file_C+"/Events")  # Add the root file to the TChain chain
        print("adding ", str_input_file_C)
        
for k in range(1, 1074):
    str_input_file_D = "/vols/cms/mc3909/bparkProductionAll_V1p0/tmp"+"/output_"+ str(k)+".root"
    if os.path.isfile(str_input_file_D) == True:
        data.Add(str_input_file_D+"/Events")  # Add the root file to the TChain chain
        print("adding ", str_input_file_D)
'''

for k in range(1, 1074):
    str_input_file_D = "/vols/cms/mc3909/bparkProductionAll_V1p3/tmp/ParkingBPH1_Run2018D-05May2019promptD-v1_MINIAOD_v1p3_generationSync"+"/output_"+ str(k)+".root"
    if os.path.isfile(str_input_file_D) == True:
        data.Add(str_input_file_D+"/Events")  # Add the root file to the TChain chain
        print("adding ", str_input_file_D)


#data.Add("/vols/cms/mc3909/bparkProductionAll_V1p3/tmp/ParkingBPH1_Run2018B-05May2019-v2_MINIAOD_v1p3_generationSync/*.root/Events")

df_data = ROOT.RDataFrame(data)

triggers = "L1_SingleMu6er1p5, L1_SingleMu7er1p5, L1_SingleMu8er1p5, L1_SingleMu9er1p5, L1_SingleMu10er1p5, L1_SingleMu12er1p5, L1_SingleMu14er1p5, L1_SingleMu16er1p5, L1_SingleMu18er1p5, L1_SingleMu22, HLT_Mu7_IP4_part0, HLT_Mu7_IP4_part1, HLT_Mu7_IP4_part2, HLT_Mu7_IP4_part3, HLT_Mu7_IP4_part4, HLT_Mu8_IP3_part0, HLT_Mu8_IP3_part1, HLT_Mu8_IP3_part2, HLT_Mu8_IP3_part3, HLT_Mu8_IP3_part4, HLT_Mu8_IP5_part0, HLT_Mu8_IP5_part1, HLT_Mu8_IP5_part2, HLT_Mu8_IP5_part3, HLT_Mu8_IP5_part4, HLT_Mu8_IP6_part0, HLT_Mu8_IP6_part1, HLT_Mu8_IP6_part2, HLT_Mu8_IP6_part3, HLT_Mu8_IP6_part4,  HLT_Mu9_IP4_part0, HLT_Mu9_IP4_part1, HLT_Mu9_IP4_part2, HLT_Mu9_IP4_part3, HLT_Mu9_IP4_part4, HLT_Mu9_IP5_part0, HLT_Mu9_IP5_part1, HLT_Mu9_IP5_part2, HLT_Mu9_IP5_part3, HLT_Mu9_IP5_part4, HLT_Mu9_IP6_part0, HLT_Mu9_IP6_part1, HLT_Mu9_IP6_part2, HLT_Mu9_IP6_part3, HLT_Mu9_IP6_part4, HLT_Mu12_IP6_part0, HLT_Mu12_IP6_part1, HLT_Mu12_IP6_part2, HLT_Mu12_IP6_part3, HLT_Mu12_IP6_part4"

df_with_triggers = df_data.Define("L1_7p0_HLT_7p0_IP4p0", "trigger_splitting(7.0, 7.0, 4.0, "+triggers+")").Define("L1_7p0_HLT_8p0_IP3p0", "trigger_splitting(7.0, 8.0, 3.0, "+triggers+")").Define("L1_8p0_HLT_7p0_IP4p0", "trigger_splitting(8.0, 7.0, 4.0, "+triggers+")").Define("L1_8p0_HLT_8p0_IP5p0", "trigger_splitting(8.0, 8.0, 5.0, "+triggers+")").Define("L1_8p0_HLT_9p0_IP5p0", "trigger_splitting(8.0, 9.0, 5.0, "+triggers+")").Define("L1_8p0_HLT_9p0_IP6p0", "trigger_splitting(8.0, 9.0, 6.0, "+triggers+")").Define("L1_9p0_HLT_8p0_IP5p0", "trigger_splitting(9.0, 8.0, 5.0, "+triggers+")").Define("L1_9p0_HLT_9p0_IP5p0", "trigger_splitting(9.0, 9.0, 5.0, "+triggers+")").Define("L1_9p0_HLT_9p0_IP6p0", "trigger_splitting(9.0, 9.0, 6.0, "+triggers+")").Define("L1_10p0_HLT_9p0_IP5p0", "trigger_splitting(10.0, 9.0, 5.0, "+triggers+")").Define("L1_10p0_HLT_9p0_IP6p0", "trigger_splitting(10.0, 9.0, 6.0, "+triggers+")").Define("L1_12p0_HLT_12p0_IP6p0", "trigger_splitting(12.0, 12.0, 6.0, "+triggers+")")

df_L1_7p0_HLT_7p0_IP4p0 = df_with_triggers.Filter("L1_7p0_HLT_7p0_IP4p0[0] == 1")
df_L1_7p0_HLT_8p0_IP3p0 = df_with_triggers.Filter("L1_7p0_HLT_8p0_IP3p0[0] == 1")
df_L1_8p0_HLT_7p0_IP4p0 = df_with_triggers.Filter("L1_8p0_HLT_7p0_IP4p0[0] == 1")
df_L1_8p0_HLT_8p0_IP5p0 = df_with_triggers.Filter("L1_8p0_HLT_8p0_IP5p0[0] == 1")
df_L1_8p0_HLT_9p0_IP5p0 = df_with_triggers.Filter("L1_8p0_HLT_9p0_IP5p0[0] == 1")
df_L1_8p0_HLT_9p0_IP6p0 = df_with_triggers.Filter("L1_8p0_HLT_9p0_IP6p0[0] == 1")
df_L1_9p0_HLT_8p0_IP5p0 = df_with_triggers.Filter("L1_9p0_HLT_8p0_IP5p0[0] == 1")
df_L1_9p0_HLT_9p0_IP5p0 = df_with_triggers.Filter("L1_9p0_HLT_9p0_IP5p0[0] == 1")
df_L1_9p0_HLT_9p0_IP6p0 = df_with_triggers.Filter("L1_9p0_HLT_9p0_IP6p0[0] == 1")
df_L1_10p0_HLT_9p0_IP5p0 = df_with_triggers.Filter("L1_10p0_HLT_9p0_IP5p0[0] == 1")
df_L1_10p0_HLT_9p0_IP6p0 = df_with_triggers.Filter("L1_10p0_HLT_9p0_IP6p0[0] == 1")
df_L1_12p0_HLT_12p0_IP6p0 = df_with_triggers.Filter("L1_12p0_HLT_12p0_IP6p0[0] == 1")


'''
print("Entries in L1_7p0_HLT_7p0_IP4p0 = ", df_L1_7p0_HLT_7p0_IP4p0.Count().GetValue())
print("Entries in L1_7p0_HLT_8p0_IP3p0 = ", df_L1_7p0_HLT_8p0_IP3p0.Count().GetValue())
print("Entries in L1_8p0_HLT_7p0_IP4p0 = ", df_L1_8p0_HLT_7p0_IP4p0.Count().GetValue())
'''


#df_trigger = df_data.Define("L1_7p0_HLT_7p0_IP4p0", "trigger_splitting(7.0, 7.0, 4.0, L1_SingleMu6er1p5, L1_SingleMu7er1p5, L1_SingleMu8er1p5, L1_SingleMu9er1p5, L1_SingleMu10er1p5, L1_SingleMu12er1p5, L1_SingleMu14er1p5, L1_SingleMu16er1p5, L1_SingleMu18er1p5, L1_SingleMu22, HLT_Mu7_IP4_part0, HLT_Mu7_IP4_part1, HLT_Mu7_IP4_part2, HLT_Mu7_IP4_part3, HLT_Mu7_IP4_part4, HLT_Mu8_IP3_part0, HLT_Mu8_IP3_part1, HLT_Mu8_IP3_part2, HLT_Mu8_IP3_part3, HLT_Mu8_IP3_part4, HLT_Mu8_IP5_part0, HLT_Mu8_IP5_part1, HLT_Mu8_IP5_part2, HLT_Mu8_IP5_part3, HLT_Mu8_IP5_part4, HLT_Mu8_IP6_part0, HLT_Mu8_IP6_part1, HLT_Mu8_IP6_part2, HLT_Mu8_IP6_part3, HLT_Mu8_IP6_part4,  HLT_Mu9_IP4_part0, HLT_Mu9_IP4_part1, HLT_Mu9_IP4_part2, HLT_Mu9_IP4_part3, HLT_Mu9_IP4_part4, HLT_Mu9_IP5_part0, HLT_Mu9_IP5_part1, HLT_Mu9_IP5_part2, HLT_Mu9_IP5_part3, HLT_Mu9_IP5_part4, HLT_Mu9_IP6_part0, HLT_Mu9_IP6_part1, HLT_Mu9_IP6_part2, HLT_Mu9_IP6_part3, HLT_Mu9_IP6_part4, HLT_Mu12_IP6_part0, HLT_Mu12_IP6_part1, HLT_Mu12_IP6_part2, HLT_Mu12_IP6_part3, HLT_Mu12_IP6_part4)")


#df_muonBPark_SV = df_data.Filter("nMuonBPark >= 2").Filter("nSV >= 1")
df_muonBPark_SV = df_data.Filter("nMuonBPark >= 2").Filter("nmuonSV >= 1")
#df_muonBPark_SV = df_data.Filter("nMuonBPark == 2").Filter("nmuonSV == 1")
#df_muonBPark_SV = df_data.Filter("nMuonBPark == 2").Filter("MuonBPark_tightId[0] == 1 && MuonBPark_tightId[1] == 1").Filter("nSV == 1")


#df_tightid = df_muonBPark_SV.Define("muonSV_mass_tightid", "SV_mass_tightid(MuonBPark_isTriggering, MuonBPark_tightId, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, SV_x, SV_y, SV_z, SV_mass, SV_chi2)")
df_tightid = df_muonBPark_SV.Define("muonSV_mass_tightid", "SV_mass_tightid(MuonBPark_matched_dr, MuonBPark_pt, MuonBPark_isTriggering, MuonBPark_dxy, MuonBPark_softId, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, muonSV_x, muonSV_y, muonSV_z, muonSV_mass, muonSV_chi2)").Define("MuonSV_mass_tightid", "muonSV_mass_tightid[0]").Define("MuonSV_pvalue", "muonSV_mass_tightid[1]")

#df_tightid1 = df_L1_7p0_HLT_7p0_IP4p0.Define("muonSV_mass_tightid", "SV_mass_tightid(MuonBPark_matched_dr, MuonBPark_pt, MuonBPark_isTriggering, MuonBPark_tightId, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, muonSV_x, muonSV_y, muonSV_z, muonSV_mass, muonSV_chi2)").Define("MuonSV_mass_tightid", "muonSV_mass_tightid[0]").Define("MuonSV_pvalue", "muonSV_mass_tightid[1]")




#df_dxy = df_muonBPark_SV.Define("muonSV_mass_dxy", "SV_mass_dxy(MuonBPark_isTriggering, MuonBPark_dxy, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, SV_x, SV_y, SV_z, SV_mass, SV_chi2)")
df_dxy = df_muonBPark_SV.Define("muonSV_mass_dxy", "SV_mass_dxy(MuonBPark_matched_dr, MuonBPark_pt, MuonBPark_isTriggering, MuonBPark_softId, MuonBPark_dxy, MuonBPark_dxyErr, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, muonSV_x, muonSV_y, muonSV_z, muonSV_mass, muonSV_chi2)").Define("MuonSV_mass_dxy", "muonSV_mass_dxy[0]")

df_dxysig = df_muonBPark_SV.Define("MuonBPark_dxysig", "MuonBPark_dxy/MuonBPark_dxyErr").Define("muonSV_mass_dxysig", "SV_mass_dxysig(MuonBPark_isTriggering, MuonBPark_dxysig, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, muonSV_x, muonSV_y, muonSV_z, muonSV_mass, muonSV_chi2)")

#h_pvalue = df_tightid.Histo1D(("hist_muonsv_pvalue_nano", "; p value; Events", 100, 0.0, 1.0), "MuonSV_pvalue")

bins = np.logspace(np.log10(0.1), np.log10(200), num = 10001)

nbins = 1999
xmin = 0.1
xmax = 200.0
logxmin = np.log10(xmin)
logxmax = np.log10(xmax)
binwidth = (logxmax-logxmin)/nbins
xbins = [xmin]
for i in range(1, nbins+1):
    xbins.append(xmin + 10**(logxmin + i*binwidth))

xbins = np.array(xbins)
print("xbins = ", xbins)

#h_muonsv_mass_tightid = df_tightid.Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.04 GeV", 1000, bins), "muonSV_mass_tightid")
h_muonsv_mass_tightid = df_tightid.Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.01 GeV", nbins, xbins), "MuonSV_mass_tightid")
h_muonsv_mass_dxy = df_dxy.Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.01 GeV", 20000, 0.0, 200.0), "MuonSV_mass_dxy")
h_muonsv_mass_dxysig = df_dxysig.Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.04 GeV", 5000, 0.0, 200.0), "muonSV_mass_dxysig")
#h_muonsv_mass = df_muonBPark_SV.Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.04 GeV", 5000, 0.0, 200), "SV_mass")

#h_muonsv_mass_tightid_logx_binning = df_tightid.Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.04 GeV", 10000, bins), "MuonSV_mass_tightid")
'''
h_muonsv_mass_tightid1 = df_L1_7p0_HLT_7p0_IP4p0.Define("muonSV_mass_tightid", "SV_mass_tightid(MuonBPark_matched_dr, MuonBPark_pt, MuonBPark_isTriggering, MuonBPark_dxy, MuonBPark_softId, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, muonSV_x, muonSV_y, muonSV_z, muonSV_mass, muonSV_chi2)").Define("MuonSV_mass_tightid", "muonSV_mass_tightid[0]").Define("MuonSV_pvalue", "muonSV_mass_tightid[1]").Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.004 GeV", 20000, 0.0, 200), "MuonSV_mass_tightid")
h_muonsv_mass_tightid2 = df_L1_7p0_HLT_8p0_IP3p0.Define("muonSV_mass_tightid", "SV_mass_tightid(MuonBPark_matched_dr, MuonBPark_pt, MuonBPark_isTriggering, MuonBPark_dxy, MuonBPark_softId, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, muonSV_x, muonSV_y, muonSV_z, muonSV_mass, muonSV_chi2)").Define("MuonSV_mass_tightid", "muonSV_mass_tightid[0]").Define("MuonSV_pvalue", "muonSV_mass_tightid[1]").Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.004 GeV", 20000, 0.0, 200), "MuonSV_mass_tightid")
h_muonsv_mass_tightid3 = df_L1_8p0_HLT_7p0_IP4p0.Define("muonSV_mass_tightid", "SV_mass_tightid(MuonBPark_matched_dr, MuonBPark_pt, MuonBPark_isTriggering, MuonBPark_dxy, MuonBPark_softId, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, muonSV_x, muonSV_y, muonSV_z, muonSV_mass, muonSV_chi2)").Define("MuonSV_mass_tightid", "muonSV_mass_tightid[0]").Define("MuonSV_pvalue", "muonSV_mass_tightid[1]").Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.004 GeV", 20000, 0.0, 200), "MuonSV_mass_tightid")
h_muonsv_mass_tightid4 = df_L1_8p0_HLT_8p0_IP5p0.Define("muonSV_mass_tightid", "SV_mass_tightid(MuonBPark_matched_dr, MuonBPark_pt, MuonBPark_isTriggering, MuonBPark_dxy, MuonBPark_softId, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, muonSV_x, muonSV_y, muonSV_z, muonSV_mass, muonSV_chi2)").Define("MuonSV_mass_tightid", "muonSV_mass_tightid[0]").Define("MuonSV_pvalue", "muonSV_mass_tightid[1]").Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.004 GeV", 20000, 0.0, 200), "MuonSV_mass_tightid")
h_muonsv_mass_tightid5 = df_L1_8p0_HLT_9p0_IP5p0.Define("muonSV_mass_tightid", "SV_mass_tightid(MuonBPark_matched_dr, MuonBPark_pt, MuonBPark_isTriggering, MuonBPark_dxy, MuonBPark_softId, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, muonSV_x, muonSV_y, muonSV_z, muonSV_mass, muonSV_chi2)").Define("MuonSV_mass_tightid", "muonSV_mass_tightid[0]").Define("MuonSV_pvalue", "muonSV_mass_tightid[1]").Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.004 GeV", 20000, 0.0, 200), "MuonSV_mass_tightid")
h_muonsv_mass_tightid6 = df_L1_8p0_HLT_9p0_IP6p0.Define("muonSV_mass_tightid", "SV_mass_tightid(MuonBPark_matched_dr, MuonBPark_pt, MuonBPark_isTriggering, MuonBPark_dxy, MuonBPark_softId, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, muonSV_x, muonSV_y, muonSV_z, muonSV_mass, muonSV_chi2)").Define("MuonSV_mass_tightid", "muonSV_mass_tightid[0]").Define("MuonSV_pvalue", "muonSV_mass_tightid[1]").Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.004 GeV", 20000, 0.0, 200), "MuonSV_mass_tightid")
h_muonsv_mass_tightid7 = df_L1_9p0_HLT_8p0_IP5p0.Define("muonSV_mass_tightid", "SV_mass_tightid(MuonBPark_matched_dr, MuonBPark_pt, MuonBPark_isTriggering, MuonBPark_dxy, MuonBPark_softId, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, muonSV_x, muonSV_y, muonSV_z, muonSV_mass, muonSV_chi2)").Define("MuonSV_mass_tightid", "muonSV_mass_tightid[0]").Define("MuonSV_pvalue", "muonSV_mass_tightid[1]").Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.004 GeV", 20000, 0.0, 200), "MuonSV_mass_tightid")
h_muonsv_mass_tightid8 = df_L1_9p0_HLT_9p0_IP5p0.Define("muonSV_mass_tightid", "SV_mass_tightid(MuonBPark_matched_dr, MuonBPark_pt, MuonBPark_isTriggering, MuonBPark_dxy, MuonBPark_softId, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, muonSV_x, muonSV_y, muonSV_z, muonSV_mass, muonSV_chi2)").Define("MuonSV_mass_tightid", "muonSV_mass_tightid[0]").Define("MuonSV_pvalue", "muonSV_mass_tightid[1]").Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.004 GeV", 20000, 0.0, 200), "MuonSV_mass_tightid")
h_muonsv_mass_tightid9 = df_L1_9p0_HLT_9p0_IP6p0.Define("muonSV_mass_tightid", "SV_mass_tightid(MuonBPark_matched_dr, MuonBPark_pt, MuonBPark_isTriggering, MuonBPark_dxy, MuonBPark_softId, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, muonSV_x, muonSV_y, muonSV_z, muonSV_mass, muonSV_chi2)").Define("MuonSV_mass_tightid", "muonSV_mass_tightid[0]").Define("MuonSV_pvalue", "muonSV_mass_tightid[1]").Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.004 GeV", 20000, 0.0, 200), "MuonSV_mass_tightid")
h_muonsv_mass_tightid10 = df_L1_10p0_HLT_9p0_IP5p0.Define("muonSV_mass_tightid", "SV_mass_tightid(MuonBPark_matched_dr, MuonBPark_pt, MuonBPark_isTriggering, MuonBPark_dxy, MuonBPark_softId, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, muonSV_x, muonSV_y, muonSV_z, muonSV_mass, muonSV_chi2)").Define("MuonSV_mass_tightid", "muonSV_mass_tightid[0]").Define("MuonSV_pvalue", "muonSV_mass_tightid[1]").Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.004 GeV", 20000, 0.0, 200), "MuonSV_mass_tightid")
h_muonsv_mass_tightid11 = df_L1_10p0_HLT_9p0_IP6p0.Define("muonSV_mass_tightid", "SV_mass_tightid(MuonBPark_matched_dr, MuonBPark_pt, MuonBPark_isTriggering, MuonBPark_dxy, MuonBPark_softId, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, muonSV_x, muonSV_y, muonSV_z, muonSV_mass, muonSV_chi2)").Define("MuonSV_mass_tightid", "muonSV_mass_tightid[0]").Define("MuonSV_pvalue", "muonSV_mass_tightid[1]").Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.004 GeV", 20000, 0.0, 200), "MuonSV_mass_tightid")
h_muonsv_mass_tightid12 = df_L1_12p0_HLT_12p0_IP6p0.Define("muonSV_mass_tightid", "SV_mass_tightid(MuonBPark_matched_dr, MuonBPark_pt, MuonBPark_isTriggering, MuonBPark_dxy, MuonBPark_softId, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, muonSV_x, muonSV_y, muonSV_z, muonSV_mass, muonSV_chi2)").Define("MuonSV_mass_tightid", "muonSV_mass_tightid[0]").Define("MuonSV_pvalue", "muonSV_mass_tightid[1]").Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.004 GeV", 20000, 0.0, 200), "MuonSV_mass_tightid")
'''
h_muonsv_mass_tightid1 = df_L1_7p0_HLT_7p0_IP4p0.Define("muonSV_mass_tightid", "SV_mass_tightid(MuonBPark_matched_dr, MuonBPark_pt, MuonBPark_isTriggering, MuonBPark_dxy, MuonBPark_softId, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, muonSV_x, muonSV_y, muonSV_z, muonSV_mass, muonSV_chi2)").Define("MuonSV_mass_tightid", "muonSV_mass_tightid[0]").Define("MuonSV_pvalue", "muonSV_mass_tightid[1]").Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.004 GeV", nbins, xbins), "MuonSV_mass_tightid")
h_muonsv_mass_tightid2 = df_L1_7p0_HLT_8p0_IP3p0.Define("muonSV_mass_tightid", "SV_mass_tightid(MuonBPark_matched_dr, MuonBPark_pt, MuonBPark_isTriggering, MuonBPark_dxy, MuonBPark_softId, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, muonSV_x, muonSV_y, muonSV_z, muonSV_mass, muonSV_chi2)").Define("MuonSV_mass_tightid", "muonSV_mass_tightid[0]").Define("MuonSV_pvalue", "muonSV_mass_tightid[1]").Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.004 GeV", nbins, xbins), "MuonSV_mass_tightid")
h_muonsv_mass_tightid3 = df_L1_8p0_HLT_7p0_IP4p0.Define("muonSV_mass_tightid", "SV_mass_tightid(MuonBPark_matched_dr, MuonBPark_pt, MuonBPark_isTriggering, MuonBPark_dxy, MuonBPark_softId, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, muonSV_x, muonSV_y, muonSV_z, muonSV_mass, muonSV_chi2)").Define("MuonSV_mass_tightid", "muonSV_mass_tightid[0]").Define("MuonSV_pvalue", "muonSV_mass_tightid[1]").Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.004 GeV", nbins, xbins), "MuonSV_mass_tightid")
h_muonsv_mass_tightid4 = df_L1_8p0_HLT_8p0_IP5p0.Define("muonSV_mass_tightid", "SV_mass_tightid(MuonBPark_matched_dr, MuonBPark_pt, MuonBPark_isTriggering, MuonBPark_dxy, MuonBPark_softId, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, muonSV_x, muonSV_y, muonSV_z, muonSV_mass, muonSV_chi2)").Define("MuonSV_mass_tightid", "muonSV_mass_tightid[0]").Define("MuonSV_pvalue", "muonSV_mass_tightid[1]").Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.004 GeV", nbins, xbins), "MuonSV_mass_tightid")
h_muonsv_mass_tightid5 = df_L1_8p0_HLT_9p0_IP5p0.Define("muonSV_mass_tightid", "SV_mass_tightid(MuonBPark_matched_dr, MuonBPark_pt, MuonBPark_isTriggering, MuonBPark_dxy, MuonBPark_softId, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, muonSV_x, muonSV_y, muonSV_z, muonSV_mass, muonSV_chi2)").Define("MuonSV_mass_tightid", "muonSV_mass_tightid[0]").Define("MuonSV_pvalue", "muonSV_mass_tightid[1]").Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.004 GeV", nbins, xbins), "MuonSV_mass_tightid")
h_muonsv_mass_tightid6 = df_L1_8p0_HLT_9p0_IP6p0.Define("muonSV_mass_tightid", "SV_mass_tightid(MuonBPark_matched_dr, MuonBPark_pt, MuonBPark_isTriggering, MuonBPark_dxy, MuonBPark_softId, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, muonSV_x, muonSV_y, muonSV_z, muonSV_mass, muonSV_chi2)").Define("MuonSV_mass_tightid", "muonSV_mass_tightid[0]").Define("MuonSV_pvalue", "muonSV_mass_tightid[1]").Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.004 GeV", nbins, xbins), "MuonSV_mass_tightid")
h_muonsv_mass_tightid7 = df_L1_9p0_HLT_8p0_IP5p0.Define("muonSV_mass_tightid", "SV_mass_tightid(MuonBPark_matched_dr, MuonBPark_pt, MuonBPark_isTriggering, MuonBPark_dxy, MuonBPark_softId, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, muonSV_x, muonSV_y, muonSV_z, muonSV_mass, muonSV_chi2)").Define("MuonSV_mass_tightid", "muonSV_mass_tightid[0]").Define("MuonSV_pvalue", "muonSV_mass_tightid[1]").Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.004 GeV", nbins, xbins), "MuonSV_mass_tightid")
h_muonsv_mass_tightid8 = df_L1_9p0_HLT_9p0_IP5p0.Define("muonSV_mass_tightid", "SV_mass_tightid(MuonBPark_matched_dr, MuonBPark_pt, MuonBPark_isTriggering, MuonBPark_dxy, MuonBPark_softId, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, muonSV_x, muonSV_y, muonSV_z, muonSV_mass, muonSV_chi2)").Define("MuonSV_mass_tightid", "muonSV_mass_tightid[0]").Define("MuonSV_pvalue", "muonSV_mass_tightid[1]").Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.004 GeV", nbins, xbins), "MuonSV_mass_tightid")
h_muonsv_mass_tightid9 = df_L1_9p0_HLT_9p0_IP6p0.Define("muonSV_mass_tightid", "SV_mass_tightid(MuonBPark_matched_dr, MuonBPark_pt, MuonBPark_isTriggering, MuonBPark_dxy, MuonBPark_softId, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, muonSV_x, muonSV_y, muonSV_z, muonSV_mass, muonSV_chi2)").Define("MuonSV_mass_tightid", "muonSV_mass_tightid[0]").Define("MuonSV_pvalue", "muonSV_mass_tightid[1]").Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.004 GeV", nbins, xbins), "MuonSV_mass_tightid")
h_muonsv_mass_tightid10 = df_L1_10p0_HLT_9p0_IP5p0.Define("muonSV_mass_tightid", "SV_mass_tightid(MuonBPark_matched_dr, MuonBPark_pt, MuonBPark_isTriggering, MuonBPark_dxy, MuonBPark_softId, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, muonSV_x, muonSV_y, muonSV_z, muonSV_mass, muonSV_chi2)").Define("MuonSV_mass_tightid", "muonSV_mass_tightid[0]").Define("MuonSV_pvalue", "muonSV_mass_tightid[1]").Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.004 GeV", nbins, xbins), "MuonSV_mass_tightid")
h_muonsv_mass_tightid11 = df_L1_10p0_HLT_9p0_IP6p0.Define("muonSV_mass_tightid", "SV_mass_tightid(MuonBPark_matched_dr, MuonBPark_pt, MuonBPark_isTriggering, MuonBPark_dxy, MuonBPark_softId, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, muonSV_x, muonSV_y, muonSV_z, muonSV_mass, muonSV_chi2)").Define("MuonSV_mass_tightid", "muonSV_mass_tightid[0]").Define("MuonSV_pvalue", "muonSV_mass_tightid[1]").Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.004 GeV", nbins, xbins), "MuonSV_mass_tightid")
h_muonsv_mass_tightid12 = df_L1_12p0_HLT_12p0_IP6p0.Define("muonSV_mass_tightid", "SV_mass_tightid(MuonBPark_matched_dr, MuonBPark_pt, MuonBPark_isTriggering, MuonBPark_dxy, MuonBPark_softId, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, muonSV_x, muonSV_y, muonSV_z, muonSV_mass, muonSV_chi2)").Define("MuonSV_mass_tightid", "muonSV_mass_tightid[0]").Define("MuonSV_pvalue", "muonSV_mass_tightid[1]").Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.004 GeV", nbins, xbins), "MuonSV_mass_tightid")


'''
df_dxy1 = df_L1_7p0_HLT_7p0_IP4p0.Define("muonSV_mass_dxy", "SV_mass_dxy(MuonBPark_matched_dr, MuonBPark_pt, MuonBPark_isTriggering, MuonBPark_softId, MuonBPark_dxy, MuonBPark_dxyErr, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, muonSV_x, muonSV_y, muonSV_z, muonSV_mass, muonSV_chi2)").Define("MuonSV_mass_dxy", "muonSV_mass_dxy[0]").Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.01 GeV", 20000, 0.0, 200.0), "MuonSV_mass_dxy")
df_dxy2 = df_L1_7p0_HLT_8p0_IP3p0.Define("muonSV_mass_dxy", "SV_mass_dxy(MuonBPark_matched_dr, MuonBPark_pt, MuonBPark_isTriggering, MuonBPark_softId, MuonBPark_dxy, MuonBPark_dxyErr, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, muonSV_x, muonSV_y, muonSV_z, muonSV_mass, muonSV_chi2)").Define("MuonSV_mass_dxy", "muonSV_mass_dxy[0]").Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.01 GeV", 20000, 0.0, 200.0), "MuonSV_mass_dxy")
df_dxy3 = df_L1_8p0_HLT_7p0_IP4p0.Define("muonSV_mass_dxy", "SV_mass_dxy(MuonBPark_matched_dr, MuonBPark_pt, MuonBPark_isTriggering, MuonBPark_softId, MuonBPark_dxy, MuonBPark_dxyErr, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, muonSV_x, muonSV_y, muonSV_z, muonSV_mass, muonSV_chi2)").Define("MuonSV_mass_dxy", "muonSV_mass_dxy[0]").Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.01 GeV", 20000, 0.0, 200.0), "MuonSV_mass_dxy")
df_dxy4 = df_L1_8p0_HLT_8p0_IP5p0.Define("muonSV_mass_dxy", "SV_mass_dxy(MuonBPark_matched_dr, MuonBPark_pt, MuonBPark_isTriggering, MuonBPark_softId, MuonBPark_dxy, MuonBPark_dxyErr, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, muonSV_x, muonSV_y, muonSV_z, muonSV_mass, muonSV_chi2)").Define("MuonSV_mass_dxy", "muonSV_mass_dxy[0]").Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.01 GeV", 20000, 0.0, 200.0), "MuonSV_mass_dxy")
df_dxy5 = df_L1_8p0_HLT_9p0_IP5p0.Define("muonSV_mass_dxy", "SV_mass_dxy(MuonBPark_matched_dr, MuonBPark_pt, MuonBPark_isTriggering, MuonBPark_softId, MuonBPark_dxy, MuonBPark_dxyErr, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, muonSV_x, muonSV_y, muonSV_z, muonSV_mass, muonSV_chi2)").Define("MuonSV_mass_dxy", "muonSV_mass_dxy[0]").Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.01 GeV", 20000, 0.0, 200.0), "MuonSV_mass_dxy")
df_dxy6 = df_L1_8p0_HLT_9p0_IP6p0.Define("muonSV_mass_dxy", "SV_mass_dxy(MuonBPark_matched_dr, MuonBPark_pt, MuonBPark_isTriggering, MuonBPark_softId, MuonBPark_dxy, MuonBPark_dxyErr, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, muonSV_x, muonSV_y, muonSV_z, muonSV_mass, muonSV_chi2)").Define("MuonSV_mass_dxy", "muonSV_mass_dxy[0]").Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.01 GeV", 20000, 0.0, 200.0), "MuonSV_mass_dxy")
df_dxy7 = df_L1_9p0_HLT_8p0_IP5p0.Define("muonSV_mass_dxy", "SV_mass_dxy(MuonBPark_matched_dr, MuonBPark_pt, MuonBPark_isTriggering, MuonBPark_softId, MuonBPark_dxy, MuonBPark_dxyErr, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, muonSV_x, muonSV_y, muonSV_z, muonSV_mass, muonSV_chi2)").Define("MuonSV_mass_dxy", "muonSV_mass_dxy[0]").Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.01 GeV", 20000, 0.0, 200.0), "MuonSV_mass_dxy")
df_dxy8 = df_L1_9p0_HLT_9p0_IP5p0.Define("muonSV_mass_dxy", "SV_mass_dxy(MuonBPark_matched_dr, MuonBPark_pt, MuonBPark_isTriggering, MuonBPark_softId, MuonBPark_dxy, MuonBPark_dxyErr, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, muonSV_x, muonSV_y, muonSV_z, muonSV_mass, muonSV_chi2)").Define("MuonSV_mass_dxy", "muonSV_mass_dxy[0]").Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.01 GeV", 20000, 0.0, 200.0), "MuonSV_mass_dxy")
df_dxy9 = df_L1_9p0_HLT_9p0_IP6p0.Define("muonSV_mass_dxy", "SV_mass_dxy(MuonBPark_matched_dr, MuonBPark_pt, MuonBPark_isTriggering, MuonBPark_softId, MuonBPark_dxy, MuonBPark_dxyErr, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, muonSV_x, muonSV_y, muonSV_z, muonSV_mass, muonSV_chi2)").Define("MuonSV_mass_dxy", "muonSV_mass_dxy[0]").Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.01 GeV", 20000, 0.0, 200.0), "MuonSV_mass_dxy")
df_dxy10 = df_L1_10p0_HLT_9p0_IP5p0.Define("muonSV_mass_dxy", "SV_mass_dxy(MuonBPark_matched_dr, MuonBPark_pt, MuonBPark_isTriggering, MuonBPark_softId, MuonBPark_dxy, MuonBPark_dxyErr, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, muonSV_x, muonSV_y, muonSV_z, muonSV_mass, muonSV_chi2)").Define("MuonSV_mass_dxy", "muonSV_mass_dxy[0]").Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.01 GeV", 20000, 0.0, 200.0), "MuonSV_mass_dxy")
df_dxy11 = df_L1_10p0_HLT_9p0_IP6p0.Define("muonSV_mass_dxy", "SV_mass_dxy(MuonBPark_matched_dr, MuonBPark_pt, MuonBPark_isTriggering, MuonBPark_softId, MuonBPark_dxy, MuonBPark_dxyErr, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, muonSV_x, muonSV_y, muonSV_z, muonSV_mass, muonSV_chi2)").Define("MuonSV_mass_dxy", "muonSV_mass_dxy[0]").Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.01 GeV", 20000, 0.0, 200.0), "MuonSV_mass_dxy")
df_dxy12 = df_L1_12p0_HLT_12p0_IP6p0.Define("muonSV_mass_dxy", "SV_mass_dxy(MuonBPark_matched_dr, MuonBPark_pt, MuonBPark_isTriggering, MuonBPark_softId, MuonBPark_dxy, MuonBPark_dxyErr, MuonBPark_vx, MuonBPark_vy, MuonBPark_vz, muonSV_x, muonSV_y, muonSV_z, muonSV_mass, muonSV_chi2)").Define("MuonSV_mass_dxy", "muonSV_mass_dxy[0]").Histo1D(("hist_muonsv_mass_nano", "; Dimuon mass (GeV); Events/0.01 GeV", 20000, 0.0, 200.0), "MuonSV_mass_dxy")
'''

h_muonsv_mass = h_muonsv_mass_tightid

print("TH1D done")

gStyle.SetOptStat(0)
gStyle.SetTextFont(42)

c1 = ROOT.TCanvas("", "", 800, 700)
c1.SetLogy()
c1.SetLogx()

h_muonsv_mass.SetLineColor(ROOT.kBlue+1)
h_muonsv_mass.GetXaxis().SetRangeUser(0.1, 200)
   
h_muonsv_mass.DrawClone("hist")

#c1.Print("/vols/cms/khl216/cmspaper_muonsvmass/twomuonswithdxysigreq_isTriggering_muonSV_new_eraD.pdf")  
#c1.Print("/vols/cms/khl216/cmspaper_muonsvmass/twotightmuons_isTriggering_muonSV_new_eraD_check3.pdf")

#c2 = ROOT.TCanvas("", "", 800, 700)
#c2.SetLogy()

#h_pvalue.DrawClone("hist")

#c1.Print("/vols/cms/khl216/cmspaper_muonsvmass/twomuonswithdxyreq_isTriggering_new_eraD_check4.pdf")
c1.Print("/vols/cms/khl216/cmspaper_muonsvmass/twotightmuons_isTriggering_muonSV_with_quality_cuts_new_eraD_logxbinning.pdf")  
#c2.Print("/vols/cms/khl216/cmspaper_muonsvmass/twotightmuons_isTriggering_muonSV_with_quality_cuts.pdf")


'''
myFile = TFile("/vols/cms/khl216/bparking_muon_vertex_mass_distributions_new3.root", "RECREATE")
print("TFile created")

myFile.WriteObject(h_muonsv_mass_tightid_logx_binning.GetPtr(), "tight_id_logx_binning")
myFile.WriteObject(h_muonsv_mass_tightid.GetPtr(), "tight_id")
#myFile.WriteObject(h_muonsv_mass_dxy.GetPtr(), "dxy_req_logx_binning")
'''

myFile = TFile("/vols/cms/khl216/bparking_muon_vertex_mass_distributions_with_triggers_fine_binning.root", "UPDATE")
print("TFile created")

myFile.WriteObject(h_muonsv_mass_tightid1.GetPtr(), "softid_L1_7p0_HLT_7p0_IP4p0_logxbinning")
myFile.WriteObject(h_muonsv_mass_tightid2.GetPtr(), "softid_L1_7p0_HLT_8p0_IP3p0_logxbinning")
myFile.WriteObject(h_muonsv_mass_tightid3.GetPtr(), "softid_L1_8p0_HLT_7p0_IP4p0_logxbinning")
myFile.WriteObject(h_muonsv_mass_tightid4.GetPtr(), "softid_L1_8p0_HLT_8p0_IP5p0_logxbinning")
myFile.WriteObject(h_muonsv_mass_tightid5.GetPtr(), "softid_L1_8p0_HLT_9p0_IP5p0_logxbinning")
myFile.WriteObject(h_muonsv_mass_tightid6.GetPtr(), "softid_L1_8p0_HLT_9p0_IP6p0_logxbinning")
myFile.WriteObject(h_muonsv_mass_tightid7.GetPtr(), "softid_L1_9p0_HLT_8p0_IP5p0_logxbinning")
myFile.WriteObject(h_muonsv_mass_tightid8.GetPtr(), "softid_L1_9p0_HLT_9p0_IP5p0_logxbinning")
myFile.WriteObject(h_muonsv_mass_tightid9.GetPtr(), "softid_L1_9p0_HLT_9p0_IP6p0_logxbinning")
myFile.WriteObject(h_muonsv_mass_tightid10.GetPtr(), "softid_L1_10p0_HLT_9p0_IP5p0_logxbinning")
myFile.WriteObject(h_muonsv_mass_tightid11.GetPtr(), "softid_L1_10p0_HLT_9p0_IP6p0_logxbinning")
myFile.WriteObject(h_muonsv_mass_tightid12.GetPtr(), "softid_L1_12p0_HLT_12p0_IP6p0_logxbinning")

'''
myFile = TFile("/vols/cms/khl216/bparking_muon_vertex_mass_distributions_with_triggers_fine_binning_dxy_new.root", "RECREATE")
print("TFile created")

myFile.WriteObject(df_dxy1.GetPtr(), "dxy_softid_L1_7p0_HLT_7p0_IP4p0")
myFile.WriteObject(df_dxy2.GetPtr(), "dxy_softid_L1_7p0_HLT_8p0_IP3p0")
myFile.WriteObject(df_dxy3.GetPtr(), "dxy_softid_L1_8p0_HLT_7p0_IP4p0")
myFile.WriteObject(df_dxy4.GetPtr(), "dxy_softid_L1_8p0_HLT_8p0_IP5p0")
myFile.WriteObject(df_dxy5.GetPtr(), "dxy_softid_L1_8p0_HLT_9p0_IP5p0")
myFile.WriteObject(df_dxy6.GetPtr(), "dxy_softid_L1_8p0_HLT_9p0_IP6p0")
myFile.WriteObject(df_dxy7.GetPtr(), "dxy_softid_L1_9p0_HLT_8p0_IP5p0")
myFile.WriteObject(df_dxy8.GetPtr(), "dxy_softid_L1_9p0_HLT_9p0_IP5p0")
myFile.WriteObject(df_dxy9.GetPtr(), "dxy_softid_L1_9p0_HLT_9p0_IP6p0")
myFile.WriteObject(df_dxy10.GetPtr(), "dxy_softid_L1_10p0_HLT_9p0_IP5p0")
myFile.WriteObject(df_dxy11.GetPtr(), "dxy_softid_L1_10p0_HLT_9p0_IP6p0")
myFile.WriteObject(df_dxy12.GetPtr(), "dxy_softid_L1_12p0_HLT_12p0_IP6p0")
'''






