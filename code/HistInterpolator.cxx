#include "HistInterpolator.h"
#include "HistCDFInterp.h"

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <TFile.h>
#include <TH1F.h>
#include <TRandom.h>
#include <TMath.h>

HistInterpolator::HistInterpolator(){ };

HistInterpolator::HistInterpolator(TString hname) :
   _hname(hname), _kmax(), _rebin(), _filenames(), _histos()
{
    _kmax.push_back(85.25); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_85.25.root"); _rebin.push_back(1);
    _kmax.push_back(85.5); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_85.5.root"); _rebin.push_back(1);
    _kmax.push_back(85.75); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_85.75.root"); _rebin.push_back(1);
    _kmax.push_back(86); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_86.0.root"); _rebin.push_back(1);
    _kmax.push_back(86.25); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_86.25.root"); _rebin.push_back(1);
    _kmax.push_back(86.5); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_86.5.root"); _rebin.push_back(1);
    _kmax.push_back(86.75); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_86.75.root"); _rebin.push_back(1);
    _kmax.push_back(87); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_87.0.root"); _rebin.push_back(1);
    _kmax.push_back(87.25); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_87.25.root"); _rebin.push_back(1);
    _kmax.push_back(87.5); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_87.5.root"); _rebin.push_back(1);
    _kmax.push_back(87.75); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_87.75.root"); _rebin.push_back(1);
    _kmax.push_back(88); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_88.0.root"); _rebin.push_back(1);
    _kmax.push_back(88.25); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_88.25.root"); _rebin.push_back(1);
    _kmax.push_back(88.5); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_88.5.root"); _rebin.push_back(1);
    _kmax.push_back(88.75); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_88.75.root"); _rebin.push_back(1);
    _kmax.push_back(89); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_89.0.root"); _rebin.push_back(1);
    _kmax.push_back(89.25); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_89.25.root"); _rebin.push_back(1);
    _kmax.push_back(89.5); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_89.5.root"); _rebin.push_back(1);
    _kmax.push_back(89.75); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_89.75.root"); _rebin.push_back(1);
    _kmax.push_back(90); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_90.0.root"); _rebin.push_back(1);
    _kmax.push_back(90.25); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_90.25.root"); _rebin.push_back(1);
    _kmax.push_back(90.5); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_90.5.root"); _rebin.push_back(1);
    _kmax.push_back(90.75); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_90.75.root"); _rebin.push_back(1);
    _kmax.push_back(91); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_91.0.root"); _rebin.push_back(1);
    _kmax.push_back(91.25); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_91.25.root"); _rebin.push_back(1);
    _kmax.push_back(91.5); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_91.5.root"); _rebin.push_back(1);
    _kmax.push_back(91.75); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_91.75.root"); _rebin.push_back(1);
    _kmax.push_back(92); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_92.0.root"); _rebin.push_back(1);
    _kmax.push_back(92.25); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_92.25.root"); _rebin.push_back(1);
    _kmax.push_back(92.5); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_92.5.root"); _rebin.push_back(1);
    _kmax.push_back(92.75); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_92.75.root"); _rebin.push_back(1);
    _kmax.push_back(93); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_93.0.root"); _rebin.push_back(1);
    _kmax.push_back(93.25); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_93.25.root"); _rebin.push_back(1);
    _kmax.push_back(93.5); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_93.5.root"); _rebin.push_back(1);
    _kmax.push_back(93.75); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_93.75.root"); _rebin.push_back(1);
    _kmax.push_back(94); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_94.0.root"); _rebin.push_back(1);
    _kmax.push_back(94.25); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_94.25.root"); _rebin.push_back(1);
    _kmax.push_back(94.5); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_94.5.root"); _rebin.push_back(1);
    _kmax.push_back(94.75); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_94.75.root"); _rebin.push_back(1);
    _kmax.push_back(95); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_95.0.root"); _rebin.push_back(1);
    _kmax.push_back(95.25); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_95.25.root"); _rebin.push_back(1);
    _kmax.push_back(95.5); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_95.5.root"); _rebin.push_back(1);
    _kmax.push_back(95.75); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_95.75.root"); _rebin.push_back(1);
    _kmax.push_back(96); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_96.0.root"); _rebin.push_back(1);
    _kmax.push_back(96.25); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_96.25.root"); _rebin.push_back(1);
    _kmax.push_back(96.5); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_96.5.root"); _rebin.push_back(1);
    _kmax.push_back(96.75); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_96.75.root"); _rebin.push_back(1);
    _kmax.push_back(97); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_97.0.root"); _rebin.push_back(1);
    _kmax.push_back(97.25); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_97.25.root"); _rebin.push_back(1);
    _kmax.push_back(97.5); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_97.5.root"); _rebin.push_back(1);
    _kmax.push_back(97.75); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_97.75.root"); _rebin.push_back(1);
    _kmax.push_back(98); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_98.0.root"); _rebin.push_back(1);
    _kmax.push_back(98.25); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_98.25.root"); _rebin.push_back(1);
    _kmax.push_back(98.5); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_98.5.root"); _rebin.push_back(1);
    _kmax.push_back(98.75); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_98.75.root"); _rebin.push_back(1);
    _kmax.push_back(99); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_99.0.root"); _rebin.push_back(1);
    _kmax.push_back(99.25); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_99.25.root"); _rebin.push_back(1);
    _kmax.push_back(99.5); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_99.5.root"); _rebin.push_back(1);
    _kmax.push_back(99.75); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_99.75.root"); _rebin.push_back(1);
    _kmax.push_back(100); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_100.0.root"); _rebin.push_back(1);
    _kmax.push_back(100.25); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_100.25.root"); _rebin.push_back(1);
    _kmax.push_back(100.5); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_100.5.root"); _rebin.push_back(1);
    _kmax.push_back(100.75); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_100.75.root"); _rebin.push_back(1);
    _kmax.push_back(101); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_101.0.root"); _rebin.push_back(1);
    _kmax.push_back(101.25); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_101.25.root"); _rebin.push_back(1);
    _kmax.push_back(101.5); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_101.5.root"); _rebin.push_back(1);
    _kmax.push_back(101.75); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_101.75.root"); _rebin.push_back(1);
    _kmax.push_back(102); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_102.0.root"); _rebin.push_back(1);
    _kmax.push_back(102.25); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_102.25.root"); _rebin.push_back(1);
    _kmax.push_back(102.5); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_102.5.root"); _rebin.push_back(1);
    _kmax.push_back(102.75); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_102.75.root"); _rebin.push_back(1);
    _kmax.push_back(103); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_103.0.root"); _rebin.push_back(1);
    _kmax.push_back(103.25); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_103.25.root"); _rebin.push_back(1);
    _kmax.push_back(103.5); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_103.5.root"); _rebin.push_back(1);
    _kmax.push_back(103.75); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_103.75.root"); _rebin.push_back(1);
    _kmax.push_back(104); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_104.0.root"); _rebin.push_back(1);
    _kmax.push_back(104.25); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_104.25.root"); _rebin.push_back(1);
    _kmax.push_back(104.5); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_104.5.root"); _rebin.push_back(1);
    _kmax.push_back(104.75); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_104.75.root"); _rebin.push_back(1);
    _kmax.push_back(105); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_105.0.root"); _rebin.push_back(1);
//    _kmax.push_back(); _filenames.push_back("BasisHists_RMC200M/rmc200M_pdf_.root"); _rebin.push_back(1);
    // _kmax.push_back(93.); _filenames.push_back("rmcspechists/rmc50M_pdf_93.root"); _rebin.push_back(1);
    // _kmax.push_back(96.); _filenames.push_back("rmcspechists/rmc50M_pdf_96.root"); _rebin.push_back(1);
    // _kmax.push_back(99.); _filenames.push_back("rmcspechists/rmc50M_pdf_99.root"); _rebin.push_back(1);

    //make a copy of the histos so you can change it without changing the original!
    for (unsigned int i=0;i<_kmax.size();++i)
    {
        TH1F *hist = getHist(_filenames[i]);
        if (hist) _histos.push_back(*hist);
        delete hist;
    }

}

//take the two neareast kmax points
TH1F* HistInterpolator::interpolateHist(double kmax)
{

    unsigned int imax = _kmax.size()-1;
    if (kmax < _kmax[0])      return interp(0,1,_kmax[0]);
    if (kmax > _kmax[imax])   return interp(imax-1,imax,_kmax[imax]);


    unsigned int index(0);
    for (index=0;index<_kmax.size()-1;++index)
        if (kmax >= _kmax[index] && kmax < _kmax[index+1] ) break;

    return interp(index,index+1,kmax,0);
}

// clone histos before rebinning them, otherwise rebinning is permanent
// caller takes ownership of histogram
TH1F* HistInterpolator::interp(int index0, int index1, double kmax, bool isRandom)
{
  double kmax0 = _kmax[index0];
  double kmax1 = _kmax[index1];
  int rebin    = std::min(_rebin[index0],_rebin[index1]);

  TH1F hist0(_histos[index0]);
  TH1F hist1(_histos[index1]);

  if (rebin >1) { hist0.Rebin(rebin); hist1.Rebin(rebin);}

  HistCDFInterp hinterp(&hist0,&hist1,kmax0,kmax1);
  hinterp.interpolate(kmax);
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
