#include "TCanvas.h"
#include "RooPlot.h"
#include "TF1.h"
#include "TMath.h"
#include "Math/Math.h"
#include "RooTFnPdfBinding.h"
#include "Fit/DataRange.h"
#include "RooConstVar.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "TRandom.h"
#include "dscb.h"
#include "dscb.cxx"

#include <iostream>
#include <fstream>
#include <sstream>

// #include <chrono>
// using namespace std::chrono;

// #include "HistCDFInterp.h"
// #include "HistCDFInterp.cxx"
// #include "HistInterpolator.h"
// #include "HistInterpolator.cxx"

// #include "clsplit.h"
// #include "clsplit.cxx"

//root -l
//.L HistCDFInterp.cxx+
//.L HistInterpolator.cxx+
//.L clsplit.cxx+
//.x FitRMC.C

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
        ostringstream recomom; // For some reason RooDataSet.Add(fitmom) is not working
        recomom << fitmom;
        myfile << recomom.str()+"\n";
      }
  }
  if (n<maxN) { HistFill(tree,n,maxN,wtlabel,myfile);}
}

void FitRMC() {
  TCanvas *c1 = new TCanvas("c1","c1");
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

  // Weighted MDC2018 RMC data
  // TFile* f = new TFile("RMCSpectra.root");
  // TH1F* wthist = (TH1F*) f->Get("hist89");
  // RooDataHist data("data","RMC", x, wthist);
  double kmaxtrue = 90.9;
  // Random Sampling of RMC Positron spectrum
  TFile* g = new TFile("RMCTrkQual.root");
  TTree* tree = (TTree*) g->Get("NewTrkQual");
  // TH1F* hist = new TH1F("hist","hist",100,82,102);
  // tree->Print();

  TRandom randg;
  double nExpected = 1500; // average number of reco rmc positron events >82 MeV. (depends on kmax)
  double nSample = nExpected+randg.Gaus()*sqrt(nExpected); // Poisson error on number of expected events

  ofstream myfile;
  const char *filename ="temprecomom.ssv"; // to store sampled dataopints temporarily
  myfile.open(filename,ios::out); // Only because RooDataSet.add was not working 

  HistFill(tree,0,nSample,"rmc909",myfile);
  myfile.close();
  // RooDataHist data("data","rmc",x,wthist);

  RooDataSet* ubdata = RooDataSet::read("temprecomom.ssv",RooArgList(x));
  // ubdata.Print("v");
  RooDataSet nubdata = *ubdata;
  rmcsmearext.fitTo(nubdata,Range(82,100),Minos(kTRUE),Extended(kTRUE));

  RooPlot* frame = x.frame();
  nubdata.plotOn(frame);
  rmcsmearext.plotOn(frame);
  rmcsmearext.paramOn(frame,Layout(.5,.8,.7));

  ostringstream ss;
  ss << fixed << setprecision(1) << kmaxtrue;
  string title = "Weighted MDC2018 RMCext (true kmax = "+ss.str()+" MeV/c )";
  frame->SetTitle(title.c_str());
  frame->Draw();
  c1->SetLogy();
  c1->Draw();
}
