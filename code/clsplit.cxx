// Implementation of RMC Closure x Pair Production x Mu2e Acceptance for kmax RMCFitting
// Subclass of RooAbsPdf
// H.Song (2019)

#include "RooFit.h"

#include "Riostream.h"
#include <math.h>
#include <iostream>
#include <vector>

#include "RooAbsReal.h"
#include "RooRealVar.h"
#include "RooAbsCategory.h"
#include <math.h>
#include "TF1.h"
#include "TMath.h"
#include "Math/Math.h"
#include "RooRandom.h"
#include "RooMath.h"

// #include "HistCDFInterp.h"
// #include "HistCDFInterp.cc"
// #include "HistInterpolator.h"
//#include "HistInterpolator.cxx"

#include "clsplit.h"

using namespace std;

ClassImp(ClSplit);

ClSplit::ClSplit(const char *name,
    const char *title,
    RooAbsReal& _x,
    RooAbsReal& _kmax,
    RooAbsReal& _Th,
    RooAbsReal& _Sl) :

  RooAbsPdf(name,title),
  x("x","Observable",this,_x),
  kmax("kmax","KinematicEndpoint",this,_kmax),
  Th("Th","Threshold",this,_Th),
  Sl("Sl","Slope",this,_Sl)
  { _hInterper = new HistInterpolator("hist"); }

ClSplit::ClSplit(const ClSplit& other, const char* name) :
  RooAbsPdf(other,name),
  x("x",this,other.x),
  kmax("kmax",this,other.kmax),
  Th("Th",this,other.Th),
  Sl("Sl",this,other.Sl),
  _hInterper(0)
  { _hInterper = new HistInterpolator(other._hInterper->GetName()); }

ClSplit::~ClSplit() { delete _hInterper; }


Double_t ClSplit::evaluate() const{
  //Check if kmax has changed and if so engage interpolation and differentiation
  if (fabs(pastkmax-kmax)<1e-6) {_kmaxchanged = false;}
  else {_kmaxchanged = true;}

  pastkmax = kmax;

  if (_kmaxchanged) {
    nhist = _hInterper->interpolateHist(kmax);
  }
  return nhist->Interpolate(x);
}
