import numpy as np
import csv
import time
execfile("/home/hsong/Dropbox/Mu2e/RMCFitting/PhysFns.py")
execfile("/home/hsong/Dropbox/Mu2e/RMCFitting/PairParam.py")
c = R.TCanvas('c','c')

clos = R.TF1('ClosureApproximation',Clos,0,1) # Closure
dscb = R.TF1('dscb',fnc_dscb,-4.,4,7) # GENMOMRES
dscb.SetParName(0,"Norm")
dscb.SetParName(1,"x0")
dscb.SetParName(2,"sigma")
dscb.SetParName(3,"ANeg")
dscb.SetParName(4,"PNeg")
dscb.SetParName(5,"APos")
dscb.SetParName(6,"PPos")
dscb.SetParameters(1.0,-.633, .322, .458, 20., 1.94, 6.7)
#Mu2e Acceptance FlatePlus
acc = R.TF1('Efficiency',erf,50,106,9)
acc.SetParName(0,"Max")
acc.SetParName(1,"Threshold")
acc.SetParName(2,"Slope")
acc.SetParameters(.223,92.1,.072)

N = 50000000

def SampleGenMom(kmax,N,binwidth,hname): #N number of samples per kmax run
    hist = R.TH1F(hname,hname,int((105-80)/binwidth),80,105)
    for i in range(N):
        if i%100000.==0: print "kmax : ", kmax, float(i)/N
        gammamom = kmax*clos.GetRandom(80./kmax,1)
        posE = gammamom*(1-SampleEpsil(13,gammamom))
        Pmom = np.sqrt(posE**2-.511**2)
        hist.Fill(Pmom,acc(Pmom))
    return hist
start = time.time()
g = R.TFile("rmcspechists/rmc50M_pdf_89.root","recreate")

hist = SampleGenMom(89,N,.2,"hist")

hist.SetLineColor(2)
hist.Draw("hist")
c.SetLogy()
c.Draw()

g.Write()

end = time.time()
print end-start
