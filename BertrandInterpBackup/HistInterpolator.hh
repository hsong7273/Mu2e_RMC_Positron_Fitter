#ifndef HistInterpolator_h
#define HistInterpolator_h

#include "HistCDFInterp.hh"
#include <TString.h>
#include <TRandom.h>
#include <TH1F.h>
#include <vector>


class HistInterpolator 
{

   public:
       
       HistInterpolator();
       HistInterpolator(TString hname, TString suffix); 
      ~HistInterpolator() {};

       TH1F* interpolateHist(double mass); 


       std::vector<int>  const& rebin() {return _rebin;}
       std::vector<TH1F> const& histo() {return _histos;}


   private:
            
      TH1F* interp(int index0, int index1, double mass, bool isRandom=0);
      TH1F* getHist(TString fname);

      TString               _hname;
      std::vector<double>   _mass;
      std::vector<int>      _rebin;
      std::vector<TString>  _filenames;
      std::vector<TH1F>     _histos;     
};

#endif
