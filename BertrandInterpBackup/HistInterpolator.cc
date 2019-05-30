#include "HistInterpolator.hh"
#include "HistCDFInterp.hh"

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <TFile.h>
#include <TH1F.h>
#include <TRandom.h>
#include <TMath.h>


HistInterpolator::HistInterpolator(TString hname, TString suffix) : 
   _hname(hname), _mass(), _rebin(), _filenames(), _histos()
{   
    _mass.push_back(0.0170); _filenames.push_back("hist/sig_0.212_"+suffix+".root"); _rebin.push_back(1);
    _mass.push_back(0.0396); _filenames.push_back("hist/sig_0.215_"+suffix+".root"); _rebin.push_back(1);
    _mass.push_back(0.0612); _filenames.push_back("hist/sig_0.22_"+suffix+".root");  _rebin.push_back(1);
    _mass.push_back(0.1336); _filenames.push_back("hist/sig_0.25_"+suffix+".root");  _rebin.push_back(1);
    _mass.push_back(0.2129); _filenames.push_back("hist/sig_0.30_"+suffix+".root");  _rebin.push_back(1);
    _mass.push_back(0.3396); _filenames.push_back("hist/sig_0.40_"+suffix+".root");  _rebin.push_back(1);
    _mass.push_back(0.4531); _filenames.push_back("hist/sig_0.50_"+suffix+".root");  _rebin.push_back(1);
    _mass.push_back(0.7196); _filenames.push_back("hist/sig_0.75_"+suffix+".root");  _rebin.push_back(1);
    _mass.push_back(0.9774); _filenames.push_back("hist/sig_1.0_"+suffix+".root");   _rebin.push_back(1);
    _mass.push_back(1.4850); _filenames.push_back("hist/sig_1.5_"+suffix+".root");   _rebin.push_back(1);
    _mass.push_back(1.9888); _filenames.push_back("hist/sig_2.0_"+suffix+".root");   _rebin.push_back(2);
    _mass.push_back(2.4911); _filenames.push_back("hist/sig_2.5_"+suffix+".root");   _rebin.push_back(2);
    _mass.push_back(2.9926); _filenames.push_back("hist/sig_3.0_"+suffix+".root");   _rebin.push_back(3);
    _mass.push_back(3.4936); _filenames.push_back("hist/sig_3.5_"+suffix+".root");   _rebin.push_back(3);
    _mass.push_back(3.9944); _filenames.push_back("hist/sig_4.0_"+suffix+".root");   _rebin.push_back(4);
    _mass.push_back(4.4950); _filenames.push_back("hist/sig_4.5_"+suffix+".root");   _rebin.push_back(4);
    _mass.push_back(4.9955); _filenames.push_back("hist/sig_5.0_"+suffix+".root");   _rebin.push_back(4);
    _mass.push_back(5.4959); _filenames.push_back("hist/sig_5.5_"+suffix+".root");   _rebin.push_back(4);
    _mass.push_back(5.9963); _filenames.push_back("hist/sig_6.0_"+suffix+".root");   _rebin.push_back(4);
    _mass.push_back(6.2965); _filenames.push_back("hist/sig_6.3_"+suffix+".root");   _rebin.push_back(4);
    _mass.push_back(6.4956); _filenames.push_back("hist/sig_6.5_"+suffix+".root");   _rebin.push_back(4);
    _mass.push_back(6.6967); _filenames.push_back("hist/sig_6.7_"+suffix+".root");   _rebin.push_back(4);
    _mass.push_back(6.9968); _filenames.push_back("hist/sig_7.0_"+suffix+".root");   _rebin.push_back(4);


    //make a copy of the histos so you can change it without changing the original!
    for (unsigned int i=0;i<_mass.size();++i)
    {
        TH1F *hist = getHist(_filenames[i]);
        if (hist) _histos.push_back(*hist);
        delete hist;
    }  
    
}



//take the two neareast mass points
TH1F* HistInterpolator::interpolateHist(double mass) 
{    
    
    unsigned int imax = _mass.size()-1;
    if (mass < _mass[0])      return interp(0,1,_mass[0]);
    if (mass > _mass[imax])   return interp(imax-1,imax,_mass[imax]);
    
    
    unsigned int index(0);
    for (index=0;index<_mass.size()-1;++index)   
        if (mass >= _mass[index] && mass < _mass[index+1] ) break;	

    return interp(index,index+1,mass,0);
}


// clone histos before rebinning them, otherwise rebinning is permanent
// caller takes ownership of histogram
TH1F* HistInterpolator::interp(int index0, int index1, double mass, bool isRandom)
{
    double mass0 = _mass[index0];
    double mass1 = _mass[index1];
    int rebin    = std::min(_rebin[index0],_rebin[index1]);

    TH1F hist0(_histos[index0]);
    TH1F hist1(_histos[index1]);

    if (rebin >1) { hist0.Rebin(rebin); hist1.Rebin(rebin);}    
        
    HistCDFInterp hinterp(&hist0,&hist1,mass0,mass1);
    hinterp.interpolate(mass);

    return static_cast<TH1F*>(hinterp.hist()->Clone());
} 


TH1F* HistInterpolator::getHist(TString fname)
{
    TFile f0(fname);
      TH1F *hist =  static_cast<TH1F*>(f0.Get(_hname));
      if (!hist) {std::cout<<"No histo named "<<_hname<<std::endl; exit(1);}
      hist->SetDirectory(0);
    f0.Close();
    
    return hist;
}    
