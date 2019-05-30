#include "HistCDFInterp.hh"

#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <TString.h>
#include <TH2.h>
#include <TAxis.h>

using namespace std;


HistCDFInterp::HistCDFInterp() {};


HistCDFInterp::HistCDFInterp(TH1F *hist1, TH1F *hist2, Double_t par1, Double_t par2) :
  _hist1(*hist1),_hist2(*hist2),_par1(par1),_par2(par2),_hmorph(0)
{
  _hist1.Smooth(1);
  _hist2.Smooth(1);
}


HistCDFInterp::~HistCDFInterp()
{
   if (_hmorph) delete _hmorph;
}



void HistCDFInterp::interpolate(double parinterp)
{

  delete _hmorph;

  if (std::fabs(parinterp-_par1)<1e-5) _hmorph = new TH1F(_hist1);
  if (std::fabs(parinterp-_par2)<1e-5) _hmorph = new TH1F(_hist2);

  if (_hist1.GetSum() <= 0 || _hist2.GetSum() <=0 )    return;
  if (fabs(_par1-_par2)<1e-5)                          return;
//  if (parinterp < _par1 || parinterp > _par2)          return;


  Int_t nb1      = _hist1.GetNbinsX();
  Int_t nb2      = _hist2.GetNbinsX();
  Double_t xmin1 = _hist1.GetXaxis()->GetXmin();
  Double_t xmin2 = _hist2.GetXaxis()->GetXmin();
  Double_t xmax1 = _hist1.GetXaxis()->GetXmax();
  Double_t xmax2 = _hist2.GetXaxis()->GetXmax();


  Double_t wt1,wt2;
  if (_par2 != _par1) { //Weight by distance from parameters (kmaxes)
    wt1 = 1. - (parinterp-_par1)/(_par2-_par1);
    wt2 = 1. + (parinterp-_par2)/(_par2-_par1);
  }
  else { wt1 = wt2 = 0.5;}


  double xminn = (xmin1 == xmin2) ? xmin1 : wt1*xmin1 + wt2*xmin2; //Configure output histogram
  double xmaxn = (xmax1 == xmax2) ? xmax1 : wt1*xmax1 + wt2*xmax2;
  int      nbn = (nb1==nb2) ? nb1 : int(wt1*nb1 + wt2*nb2); //Just use histograms with identical binning
  //Don't want to round

  Double_t dx1 = (xmax1-xmin1)/double(nb1); //binwidths
  Double_t dx2 = (xmax2-xmin2)/double(nb2);
  Double_t dx  = (xmaxn-xminn)/double(nbn);

  if (wt1 == 0) {xminn = xmin2; xmaxn = xmax2; nbn = nb2;}
  if (wt2 == 0) {xminn = xmin1; xmaxn = xmax1; nbn = nb1;}


  Float_t *dist1    = _hist1.GetArray(); //convert histograms into ROOT arrays
  Float_t *dist2    = _hist2.GetArray();
  //initialize some vectors
  std::vector<double> sigdis1(1+nb1,0); //use as normalized cumulative distributions
  std::vector<double> sigdis2(1+nb2,0);
  std::vector<double> sigdisn(2+nb1+nb2,0);
  std::vector<double> xdisn(2+nb1+nb2,0);
  std::vector<double> sigdisf(1+nbn,0);

  for(Int_t i=1;i<nb1+1;i++) sigdis1[i] = dist1[i]; // start at 1 to avoid underflow bin
  for(Int_t i=1;i<nb2+1;i++) sigdis2[i] = dist2[i];


  Double_t total(0); //Calculate hist integrals for weighted normalization
  for(Int_t i=0;i<nb1+1;i++) total += sigdis1[i]; //start at 0 to nb1
  //"normalized" cumulative distribution
  for(Int_t i=1;i<nb1+1;i++) sigdis1[i] = sigdis1[i]/total + sigdis1[i-1];
  double norm1(total);

  total = 0.;
  for(Int_t i=0;i<nb2+1;i++) total += sigdis2[i];
  for(Int_t i=1;i<nb2+1;i++) sigdis2[i] = sigdis2[i]/total + sigdis2[i-1];
  double norm2(total);


  Int_t ix1l(nb1); // Check if Cumulative transform went wrong??
  Int_t ix2l(nb2); // Find where the cumulative dist starts to go down
  while(sigdis1[ix1l-1] >= sigdis1[ix1l]) ix1l = ix1l - 1;
  while(sigdis2[ix2l-1] >= sigdis2[ix2l]) ix2l = ix2l - 1;


  Int_t ix1(-1),ix2(-1); //Find where the cumulative dist starts to go up
  do {ix1 = ix1 + 1;} while(sigdis1[ix1+1] <= sigdis1[0]);
  do {ix2 = ix2 + 1;} while(sigdis2[ix2+1] <= sigdis2[0]);


  Int_t nx3 = 0;
  Double_t x1 = xmin1 + double(ix1)*dx1;
  Double_t x2 = xmin2 + double(ix2)*dx2;
  Double_t x = wt1*x1 + wt2*x2; // Calculate interpolated "turn on" position
  xdisn[nx3]   = x;
  sigdisn[nx3] = 0;



  Double_t yprev = -1;
  Double_t y,x20,x21,y20,y21;
  Double_t x10,x11,y10,y11;

  while (ix1 < ix1l || ix2 < ix2l) {

      Int_t i12type = -1;
      if ((sigdis1[ix1+1] <= sigdis2[ix2+1] || ix2 == ix2l) && ix1 < ix1l)
      {
	ix1 = ix1 + 1;
	while(sigdis1[ix1+1] < sigdis1[ix1] && ix1 < ix1l) ix1 = ix1 + 1;
	i12type = 1;

      } else if (ix2 < ix2l)
      {
	ix2 = ix2 + 1;
	while(sigdis2[ix2+1] < sigdis2[ix2] && ix2 < ix2l) ix2 = ix2 + 1;
	i12type = 2;
      }

      if (i12type == 1)
      {
	x1  = xmin1 + double(ix1)*dx1 ;
	y   = sigdis1[ix1];
	x20 = double(ix2)*dx2 + xmin2;
	x21 = x20 + dx2;
	y20 = sigdis2[ix2];
	y21 = sigdis2[ix2+1];

	if (y21 > y20) x2 = x20 + (x21-x20)*(y-y20)/(y21-y20);
	else x2 = x20;

      } else {
	x2  = xmin2 + double(ix2)*dx2 ;
	y   = sigdis2[ix2];
	x10 = double(ix1)*dx1 + xmin1;
	x11 = x10 + dx1;
	y10 = sigdis1[ix1];
	y11 = sigdis1[ix1+1];
	if (y11 > y10) x1 = x10 + (x11-x10)*(y-y10)/(y11-y10);
	else x1 = x10;
      }

      x = wt1*x1 + wt2*x2;
      if (y >= yprev) {
	nx3 = nx3+1;
	yprev = y;
	xdisn[nx3] = x;
	sigdisn[nx3] = y;
      }
  }

  x = xminn + double(nbn)*dx;
  Int_t ix = nbn;

  while(x >= xdisn[nx3])
  {
    sigdisf[ix] = sigdisn[nx3];
    ix = ix-1;
    x  = xminn + double(ix)*dx;
  }
  Int_t ixl = ix + 1;


  ix = 0;
  x  = xminn + double(ix+1)*dx;
  while(x <= xdisn[0])
  {
    sigdisf[ix] = sigdisn[0];
    ix = ix+1;
    x = xminn + double(ix+1)*dx;
  }
  Int_t ixf = ix;




  Int_t ix3 = 0;
  for (ix=ixf;ix<ixl;ix++)
  {
     x = xminn + double(ix)*dx;
     if (x < xdisn[0]) {
       y = 0;
     } else if (x > xdisn[nx3]) {
       y = 1.;
     } else {
       while(xdisn[ix3+1] <= x && ix3 < 2*nbn) ix3 = ix3 + 1;
       if      (xdisn[ix3+1]-x > 1.1*dx2)  y = sigdisn[ix3+1];
       else if (xdisn[ix3+1] > xdisn[ix3]) y = sigdisn[ix3] + (sigdisn[ix3+1]-sigdisn[ix3])*(x-xdisn[ix3])/(xdisn[ix3+1]-xdisn[ix3]);
       else     y = 0;
     }
     sigdisf[ix] = y;
  }

  _hmorph = new TH1F(Form("TH1F-interpolated%i",rand()),"Interpolated histogram",nbn,xminn,xmaxn);

  double  norm = (norm1 == norm2) ?  norm1 : wt1*norm1 + wt2*norm2;
  for(Int_t ixx=nbn-1;ixx>-1;ixx--)
  {
    x = xminn + double(ixx)*dx;
    y =  sigdisf[ixx+1]-sigdisf[ixx];
    _hmorph->SetBinContent(ixx+1,y*norm);
  }

}
