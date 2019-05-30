/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            * 
 *****************************************************************************/

#ifndef ACC
#define ACC

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"
 
class Acc : public RooAbsPdf {
public:
  Acc() {} ; 
  Acc(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooAbsReal& _Th,
	      RooAbsReal& _Sl);
  Acc(const Acc& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new Acc(*this,newname); }
  inline virtual ~Acc() { }

protected:

  RooRealProxy x ;
  RooRealProxy Th ;
  RooRealProxy Sl ;
  
  Double_t evaluate() const ;

private:

  ClassDef(Acc,1) // Your description goes here...
};
 
#endif
