#ifndef HISTINTERPOLATOR
#define HISTINTERPOLATOR

#include "HistCDFInterp.h"
#include <TString.h>
#include <TRandom.h>
#include <TH1F.h>
#include <vector>


class HistInterpolator
{

   public:

       HistInterpolator();
       HistInterpolator(TString hname);
       virtual ~HistInterpolator() {};

       TH1F* interpolateHist(double kmax);


       std::vector<int>  const& rebin() {return _rebin;}
       std::vector<TH1F> const& histo() {return _histos;}
       TString GetName() {return _hname;}

   private:
      double pastkmax = -1000;
      bool kmax_changed = true;
      TH1F* interp(int index0, int index1, double kmax, bool isRandom=0);
      TH1F* getHist(TString fname);

      TString               _hname;
      std::vector<double>   _kmax;
      std::vector<int>      _rebin;
      std::vector<TString>  _filenames;
      std::vector<TH1F>     _histos;
};

#endif
