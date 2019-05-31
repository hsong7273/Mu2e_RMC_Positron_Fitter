import csv
import ROOT as R

c = R.TCanvas('c','c')
mom = R.RooRealVar("fitmom","fitmom",85,95)
errorlo = R.RooRealVar("errorlo","errorlo",-2,2)
errorhi = R.RooRealVar("errorhi","errorhi",-2,2)
pull = R.RooRealVar("pull","pull",-5,5)
pulldata = R.RooDataSet("pulldata","rmc90.9",R.RooArgSet(pull))

filename = "fits_kmax_90.9.ssv"

dh = R.RooDataSet.read("/home/hsong/Dropbox/Mu2e/RMCFitting/CDF_Interpolation/fitdata/"+filename,R.RooArgList(mom,errorlo,errorhi,pull))
#x = R.RooRealVar("x","x",-7,7)
mean = R.RooRealVar("mean","mean",0,-5,5)
sigma = R.RooRealVar("sigma","sigma",1,0,30)

gauss = R.RooGaussian("gauss","gauss",pull,mean,sigma)

gauss.fitTo(dh)
mean.Print()
sigma.Print()

frame = pull.frame()
dh.plotOn(frame)
gauss.plotOn(frame)
frame.SetTitle("KmaxFits kmax=90.9 MeV (N=2,000)")
#frame.SetXTitle("Pull")
frame.SetXTitle("Kmax (MeV/c)")

gauss.paramOn(frame,R.RooFit.Layout(0.1,0.45,0.9))

frame.Draw()
c.Update()
c.Draw()
