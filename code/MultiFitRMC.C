#include "TCanvas.h"
#include "RooPlot.h"
#include "TF1.h"
#include "TMath.h"
#include "Math/Math.h"
#include "RooTFnPdfBinding.h"
#include "Fit/DataRange.h"
#include "RooConstVar.h"
#include "RooRealVar.h"
#include "TRandom.h"
#include "/home/hsong/Dropbox/Mu2e/RMCFitting/dscb.h"
#include "/home/hsong/Dropbox/Mu2e/RMCFitting/dscb.cxx"

#include <iostream>
#include <fstream>
// #include <sstream>\

//ROOT Macro for TOYMC & Pull Studies

// #include "HistCDFInterp.h"
// #include "HistCDFInterp.cxx"
// #include "HistInterpolator.h"
// #include "HistInterpolator.cxx"

// #include "clsplit.h"
// #include "clsplit.cxx"
#include <chrono>
using namespace std::chrono;

using namespace RooFit;
using namespace std;

void HistFill( TTree* tree, int n, int maxN, string wtlabel, ofstream& myfile) {

  Double_t fitmom, newtrkqual, status, t0, d0, td, om;
  tree->SetBranchAddress("FitMom",&fitmom);
  tree->SetBranchAddress("newtrkqual",&newtrkqual);
  auto de = tree->GetBranch("de");
  auto demc = tree->GetBranch("demc");
  auto evtwt = tree->GetBranch("evtwt");
  auto nevent = tree->GetEntries();
  TRandom rand(n);

  for (Int_t i=0;i<nevent;i++){ //Loop through the tree
      tree->GetEvent(i);
      status = de->GetLeaf("status")->GetValue();
      t0 = de->GetLeaf("t0")->GetValue();
      om = de->GetLeaf("om")->GetValue();
      d0 = de->GetLeaf("d0")->GetValue();
      td = de->GetLeaf("td")->GetValue();
      if (
        n<maxN &&
        newtrkqual>.4 &&
        fitmom>82 &&
        demc->GetLeaf("proc")->GetValue()==13 &&
        t0>700 && t0<1695 &&
        td>.577 && td<1 &&
        d0>-105 && d0<80 &&
        d0+2./om<-450 && d0+2./om>-650 &&
        evtwt->GetLeaf(wtlabel.c_str())->GetValue()>rand.Uniform()*1.45
      ) {
        n+=1;
        ostringstream recomom;
        recomom << fitmom;
        myfile << recomom.str()+"\n";
      }
  }
  if (n<maxN) { HistFill(tree,n,maxN,wtlabel,myfile);}
}

void MultiFitRMC() {
  auto start = high_resolution_clock::now();
  // TCanvas *c1 = new TCanvas("c1","c1");
  RooRealVar x("RecoMom","RecoMom",82,102);
  //GENMOMRES FlatEPlus
  RooRealVar x0("x0","x0",-.633) ;
  RooRealVar sigma("sigma","sigma",0.322) ;
  RooRealVar ANeg("ANeg","ANeg",.458) ;
  RooRealVar PNeg("PNeg","PNeg",20.) ;
  RooRealVar APos("APos","APos",1.94) ;
  RooRealVar PPos("PPos","PPos",6.7) ;
  dscb res("momres","momres",x,x0,sigma,ANeg,PNeg,APos,PPos);

  RooRealVar kmax("kmax","kmax",93,85,105);
  //Mu2e Acceptance Curve
  RooRealVar Th("Th","Th",92.1); //
  RooRealVar Sl("Sl","Sl",.0722); //

  ClSplit clsplit("clsplit","clsplit",x,kmax,Th,Sl);

  RooFFTConvPdf rmcsmear("rmcsmear","Clos (X) Res",x,clsplit,res);
  rmcsmear.setBufferFraction(1.4);

  RooRealVar nrmc("nrmcreco","number of reconstructed rmc positron events",2700,0,10000);
  RooExtendPdf rmcsmearext("rmcsmearext","rmcsmearext",rmcsmear,nrmc);

  double kmaxtrue = 91.9;

  TFile* g = new TFile("/home/hsong/Mu2e/RMCTrkQual.root");
  TTree* tree = (TTree*) g->Get("NewTrkQual");

  TRandom randg;

  ofstream myfile; // final fitdata file
  const char *filename ="fitdata/fits_kmax_91.9.ssv";
  myfile.open(filename,ios::out);
  myfile << "fitted kmax LoError HiError Pull \n";

  double nExpected = 2500;
  double Nfits = 20;
  for (Int_t i=0; i<Nfits; i++){
    cout << i << "/" << Nfits << endl;

    double nSample = nExpected+randg.Gaus()*sqrt(nExpected); // Poisson error on number of expected events
    ofstream datfile;
    const char *filename ="temprecomom.ssv";
    datfile.open(filename,ios::out);
    HistFill(tree,0,nSample,"rmc919",datfile);
    datfile.close();
    RooDataSet* ubdata = RooDataSet::read("/home/hsong/Dropbox/Mu2e/RMCFitting/CDF_Interpolation/temprecomom.ssv",RooArgList(x));
    RooDataSet data = *ubdata;
    rmcsmearext.fitTo(data,Range(81,100),Minos(kTRUE),Extended(kTRUE));
    ostringstream fittedkmax;
    fittedkmax << kmax.getValV();
    ostringstream ErrorHi;
    ErrorHi << kmax.getAsymErrorHi();
    ostringstream ErrorLo;
    ErrorLo << kmax.getAsymErrorHi();
    //Calculate Pull
    double pull = 0;
    double err = 0;
    if (kmax.getValV()>kmaxtrue) { pull = (kmax.getValV()-kmaxtrue)/kmax.getAsymErrorHi(); }
    else { pull = -1*(kmax.getValV()-kmaxtrue)/kmax.getAsymErrorLo(); }
    ostringstream pullstr;
    pullstr << pull;
    myfile << fittedkmax.str()+" "+ErrorLo.str()+" "+ErrorHi.str()+" "+pullstr.str()+" \n";
  }
  myfile.close();
  auto stop = high_resolution_clock::now();

  auto duration = duration_cast<seconds>(stop - start);
  cout << duration.count() << endl;

}
