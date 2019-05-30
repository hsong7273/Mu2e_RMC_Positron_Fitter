import ROOT as R
import numpy as np

f = R.TFile("/home/hsong/Mu2e/RMCTrkQual.root") # Weighted MDC2018 flatmugammamix retrained trkqual
tree = f.Get("NewTrkQual")


g = R.TFile("RMCSpectra.root","recreate")
rmc89 = R.TH1F("89","89",240,80,104)
rmc909 = R.TH1F("90.9","90.9",240,80,104)
rmc919 = R.TH1F("91.9","91.9",240,80,104)
rmc93 = R.TH1F("93","93",240,80,104)
rmc94 = R.TH1F("94","94",240,80,104)
rmc95 = R.TH1F("95","95",240,80,104)
rmc96 = R.TH1F("96","96",240,80,104)
rmc97 = R.TH1F("97","97",240,80,104)
rmc98 = R.TH1F("98","98",240,80,104)
rmc99 = R.TH1F("99","99",240,80,104)
rmc100 = R.TH1F("100","100",240,80,104)
rmc101 = R.TH1F("101","101",240,80,104)


for event in tree: #OLD DIO CUTS
    if (
        event.newtrkqual>.4 and
        event.t0>=700 and event.t0<=1695 and
        event.td>=.577 and event.td<=1 and
        event.d0>=-105 and event.d0<=80 and
        event.d0+2./event.om<=-450 and event.d0+2./event.om>=-680
    ):
        if (event.rmc89>0):   rmc89.Fill(event.FitMom,event.rmc89/89)
        if (event.rmc909>0):  rmc909.Fill(event.FitMom,event.rmc909/90.9)
        if (event.rmc919>0):  rmc919.Fill(event.FitMom,event.rmc919/91.9)
        if (event.rmc93>0):   rmc93.Fill(event.FitMom,event.rmc93/93)
        if (event.rmc94>0):   rmc94.Fill(event.FitMom,event.rmc94/94)
        if (event.rmc95>0):   rmc95.Fill(event.FitMom,event.rmc95/95)
        if (event.rmc96>0):   rmc96.Fill(event.FitMom,event.rmc96/96)
        if (event.rmc997>0):  rmc97.Fill(event.FitMom,event.rmc997/97)
        if (event.rmc98>0):   rmc98.Fill(event.FitMom,event.rmc98/98)
        if (event.rmc99>0):   rmc99.Fill(event.FitMom,event.rmc99/99)
        if (event.rmc100>0):  rmc100.Fill(event.FitMom,event.rmc100/100)
        if (event.rmc101>0):  rmc101.Fill(event.FitMom,event.rmc101/101)

NPOT = 3.6*10**20
pmustop = .00181
Alfucap = .61
AlfRMC = 1.43*10**(-5)
ngen = 9.58*10**8
ScaleFactor = 31*NPOT*pmustop*Alfucap*AlfRMC/ngen

g.Write()
