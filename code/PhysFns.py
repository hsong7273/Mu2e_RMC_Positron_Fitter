import ROOT as R

def PosP(x): #PairProduction in High Energy Limit
    result = 1-(4/3.)*x[0]*(1-x[0])
    return result
#Momentum Resolution double-sided crystall ball function
def fnc_dscb(xx,par):
    x = xx[0] #xx is indexed to move past double buffer to the double
    #Gaussian Core
    N = par[0]
    mu = par[1]
    sig = par[2]
    #transition parameters
    a1 = par[3]
    p1 = par[4]
    a2 = par[5]
    p2 = par[6]
    #Calculations
    u = (x-mu)/sig
    A1 = R.TMath.Power(p1/R.TMath.Abs(a1),p1)*R.TMath.Exp(-a1*a1/2)
    A2 = R.TMath.Power(p2/R.TMath.Abs(a2),p2)*R.TMath.Exp(-a2*a2/2)
    B1 = p1/R.TMath.Abs(a1) - R.TMath.Abs(a1)
    B2 = p2/R.TMath.Abs(a2) - R.TMath.Abs(a2)
    if u < -1*(a1):
        result = A1*R.TMath.Power(B1-u,-p1)
    elif u < a2:
        result = R.TMath.Exp(-u*u/2)
    else:
        result = A2*R.TMath.Power(B2+u,-p2)
    return result

def Clos(x): #Closure Approximation
    result = (1-2*x[0]+2*x[0]**2)*x[0]*(1-x[0])**2
    return result

def erf(xx,par): #Mu2e Tracker Acceptace
    x = xx[0] #xx is indexed to move past double buffer to the double

    M = par[0]
    Th = par[1]
    Sl = par[2]

    #Calculations
    result = M*(R.TMath.Erf(Sl*(x-Th))+1)/2.
    return result

def OldClosxSplit(x,par):
    result = 0
    k = par[0]
    Th = par[1]
    Sl = par[2]
    pmax = k-.511
    acc = (R.TMath.Erf(Sl*(x[0]-Th))+1)/2.
    z = x[0]/pmax
    result = acc*(z>0)*(z<1.)*((-28./45)*R.TMath.Power(z,5)+(17./6)*R.TMath.Power(z,4)-7*R.TMath.Power(z,3)+(16./3)*R.TMath.Power(z,2)*R.TMath.Log(z)+(14./9)*R.TMath.Power(z,2)+(4./3)*z*R.TMath.Log(z)+3*z+(7./30))
    return result
