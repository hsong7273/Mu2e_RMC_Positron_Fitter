#ifndef HistCDFInterp_h
#define HistCDFInterp_h

#include <iomanip>
#include <map>
#include <vector>
#include <string>
#include <TH1.h>
#include <TString.h>


class HistCDFInterp {


   public:
   
      HistCDFInterp();
      HistCDFInterp(TH1F *hist1, TH1F *hist2, Double_t par1, Double_t par2);
      ~HistCDFInterp(); 
    
      void interpolate(double parinterp);
      TH1F* hist() const  {return _hmorph;}
      double par1() const {return _par1;}
      double par2() const {return _par2;}
     
      


   private:
      
      TH1F  _hist1;
      TH1F  _hist2;
      double _par1;
      double _par2;
      TH1F  *_hmorph;     
};

#endif
