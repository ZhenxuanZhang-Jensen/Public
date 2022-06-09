import uproot
from matplotlib import pyplot as plt
import numpy as np
import ROOT
from ROOT import gROOT 
import uproot
# create data for ROC curve
# XGboost input file name: XGboost_diphoMVA.root
# TMVA input file name: TMVA_diphoMVA.root
# JT TMVA input file name: DiphotonBDT_BDT.root
# XGboost_diphoMVA_t4.root
# test_sig, test_bkg_pfff, test_bkg_pp
event_signal = uproot.open("/eos/user/z/zhenxuan/Hgg_mass/MiniTree/XGboost_diphoMVA_t8.root:ggh_125_13TeV_UntaggedTag")
event_bkg1 = uproot.open("/eos/user/z/zhenxuan/Hgg_mass/MiniTree/XGboost_diphoMVA_t8.root:DataDriven_QCD")
event_bkg2 = uproot.open("/eos/user/z/zhenxuan/Hgg_mass/MiniTree/XGboost_diphoMVA_t8.root:pp")

# event_signal = uproot.open("/eos/user/z/zhenxuan/Hgg_mass/MiniTree/TMVA_diphoMVA.root:ggh_125_13TeV_UntaggedTag")
# event_bkg1 = uproot.open("/eos/user/z/zhenxuan/Hgg_mass/MiniTree/New_MCpp_DataDriven_QCD_SFs_sEoEWgt_2DpTWgt.root:DataDriven_QCD")
# event_bkg2 = uproot.open("/eos/user/z/zhenxuan/Hgg_mass/MiniTree/New_MCpp_DataDriven_QCD_SFs_sEoEWgt_2DpTWgt.root:pp")
# DiphotonMVA_self = 1. / ( 1. + exp( 0.5*log( 2./(DiphotonMVA_self+1.) - 1 ) ) )
diphoMVA_signal = event_signal['DiphotonMVA_self'].array(library = 'np')
diphoMVA_bkg1 = event_bkg1['DiphotonMVA_self'].array(library = 'np')
diphoMVA_bkg2 = event_bkg2['DiphotonMVA_self'].array(library = 'np')
# diphoMVA_signal = 1. / ( 1. + np.exp( 0.5*np.log( 2./(diphoMVA_signal+1.) - 1 ) ) )
# diphoMVA_bkg1 = 1. / ( 1. + np.exp( 0.5*np.log( 2./(diphoMVA_bkg1+1.) - 1 ) ) )
# diphoMVA_bkg2 = 1. / ( 1. + np.exp( 0.5*np.log( 2./(diphoMVA_bkg2+1.) - 1 ) ) )
# diphoMVA_signal = event_signal['DiphotonMVA_self'].array(library = 'np')
# diphoMVA_bkg1 = event_bkg1['DiphotonMVA_self'].array(library = 'np')
# diphoMVA_bkg2 = event_bkg2['DiphotonMVA_self'].array(library = 'np')

chunk_sig_df = event_signal.arrays(['weight','vtxprob','sigmarv','sigmawv'],library='pd')
weight_signal =(chunk_sig_df['weight']*(chunk_sig_df['vtxprob']*1./chunk_sig_df['sigmarv']+(1-chunk_sig_df['vtxprob'])*1./chunk_sig_df['sigmawv']))
weight_signal.clip(lower = 0)
weight_bkg1 = event_bkg1['weight'].array(library = 'np')
weight_bkg1[weight_bkg1<0] = 0
weight_bkg2 = event_bkg2['weight'].array(library = 'np')
weight_bkg2[weight_bkg2<0] = 0
Norm_SFs_bkg1 = event_bkg1['Norm_SFs'].array(library = 'np')
Norm_SFs_bkg2 = event_bkg2['Norm_SFs'].array(library = 'np')
diphoMVA_bkg = np.append(diphoMVA_bkg1,diphoMVA_bkg2)
weight_bkg = np.append(weight_bkg1,weight_bkg2)
reweight_sig_bkg = sum(weight_bkg) / sum(weight_signal)
weight_signal  = reweight_sig_bkg * weight_signal
Norm_SFs_bkg = np.append(Norm_SFs_bkg1,Norm_SFs_bkg2)
# pyroot histogram method
import ROOT
from ROOT import gROOT 
gROOT.GetListOfCanvases().Delete()

c = ROOT.TCanvas()
h_sig = ROOT.TH1F("sig","DiphoMVA_sig",100,-1,1)
h_bkg_pp = ROOT.TH1F("bkg_pp","DiphoMVA_pp",100,-1,1)
h_bkg_pfff = ROOT.TH1F("bkg_pfff","DiphoMVA_pfff",100,-1,1)
for i in range(len(diphoMVA_bkg1)):
    h_bkg_pfff.Fill(diphoMVA_bkg1[i],weight_bkg1[i]*Norm_SFs_bkg1[i])
    
    
for i in range(len(diphoMVA_bkg2)):
    h_bkg_pp.Fill(diphoMVA_bkg2[i],weight_bkg2[i]*Norm_SFs_bkg2[i])

h_bkg_pp.SetFillColor(38)
h_bkg_pfff.SetFillColor(30)
h_bkg_pfff.SetLineWidth(2)
h_bkg_pp.SetLineWidth(2)
h_bkg_pp.SetLineStyle(1)
h_bkg_pfff.SetLineStyle(1)
hs = ROOT.THStack('hs','DiphoMVA_xgboost')
# hs = ROOT.THStack('hs','DiphoMVA_TMVA')
hs.Add(h_bkg_pp)
hs.Add(h_bkg_pfff)
for i in range(len(diphoMVA_signal)):
    h_sig.Fill(diphoMVA_signal[i],weight_signal[i])

h_sig.SetLineWidth(2)
h_sig.SetLineColor(1)
h_sig.SetFillColorAlpha(1,0.35)
# bymax=h_sig.GetBinContent(100)+100000 # for TMVA max Y
sig_ymax=h_sig.GetBinContent(h_sig.GetMaximumBin()) # for XGboost max Y
# bkg_ymax=hs.GetBinContent(hs.GetMaximumBin())*1.01 # for XGboost max Y
# h_sig.SetMaximum(sig_ymax)
hs.SetMaximum(sig_ymax)
hs.Draw('histo')
h_sig.Draw('same,histo')
gROOT.GetListOfCanvases().Draw()
legend = ROOT.TLegend(0.2,0.2,0.4,0.4)
legend.AddEntry(h_sig,"sig","f")
legend.AddEntry(h_bkg_pp,"#gamma#gamma","f")
legend.AddEntry(h_bkg_pfff,"Data-driven QCD","f")
legend.SetTextFont(132)
legend.SetHeader("Legend")
legend.SetTextFont(42)
legend.Draw("same")
# gROOT.SetOptionStat(0)