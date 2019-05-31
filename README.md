# Mu2e_RMC_Positron_Fitter
compile classes
- root -l
- .L code/HistCDFInterp.cxx+
- .L code/HistInterpolator.cxx+
- .L code/clsplit.cxx+
Run ROOT macros
- .x code/FitRMC.C
- .x code/MultiFitRMC.C


PairParam.py & PhysFns.py are used by CDF_Sampling.py to create the basis histograms

RooDataSet::read wants space separated values as inputs
 - The .ssv files are used by kmaxpullfit to create the pull/kmax distributions
