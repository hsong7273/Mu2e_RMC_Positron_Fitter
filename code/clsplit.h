// Implementation of RMC Closure x Pair Production x Mu2e Acceptance for kmax RMCFitting
// Subclass of RooAbsPdf
// H.Song (2019)

#ifndef RMC_CLSPLIT
#define RMC_CLSPLIT

#include <iostream>
#include <vector>

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooTrace.h"

#include "HistInterpolator.h"


using namespace std;
class RooRealVar;

class ClSplit : public RooAbsPdf {
public :
  ClSplit() { };
  ClSplit(const char *name,
    const char *title,
    RooAbsReal& _x,
    RooAbsReal& _kmax,
    RooAbsReal& _Th,
    RooAbsReal& _Sl);

  ClSplit(const ClSplit& other, const char* name=0);
  virtual TObject* clone(const char* newname) const { return new ClSplit(*this, newname); }
  virtual ~ClSplit();

  mutable HistInterpolator* _hInterper;
  mutable TH1F* nhist;
  mutable bool _kmaxchanged = true;
  mutable double pastkmax = -1000;
protected:
  RooRealProxy x;
  RooRealProxy kmax;
  RooRealProxy Th;
  RooRealProxy Sl;

  Double_t evaluate() const;

private:
//ClosureApproximation convoluted with Pair Production multiplied by Mu2e e+ Acceptance
  ClassDef(ClSplit,1)
};

#endif
