import ROOT as R
import numpy as np

def F(z,gammaE): #Coulomb Correction function
    aZ = z/137.
    fc = (aZ**2)*(1/(1+(aZ)**2)+.020206-.0367*aZ**2+.0083*aZ**4-.002*aZ**6)
    if gammaE>50.:
        return (8./3)*np.log(z)+8*fc
    else:
        return (8./3)*np.log(z)
#screening variable & Screening functions
def ScreenFunction1(sval):
    if sval>1:
        sval = 42.24-8.368*np.log(sval+.952)
    else:
        sval = 42.392-sval*(7.796-1.961*sval)
    return sval
def ScreenFunction2(sval):
    if sval>1:
        sval = 42.24-8.368*np.log(sval+.952)
    else:
        sval = 41.734-sval*(6.484-1.25*sval)
    return sval

#Epsil = electronEnergy/photonEnergy
#Use photonE*Epsil output of this function for lower energy electron
#Use photonE*(1-Epsil) for higher energy electron
def SampleEpsil(Z,gammaE): #Sample Positron/Electron Energies
    greject = 0.
    epsil = 0.
    epsil0 = .511/gammaE
    if gammaE<2.:
        epsil = epsil0+(.5-epsil0)*np.random.rand()
    else:
        FZ = F(Z,gammaE) #Coulomb Correction Formula
        #limits of the screening variable
        screenfac = 136.*epsil0/(Z**(1./3))
        screenmax = np.exp((42.24-FZ)/8.368) - .952
        screenmin = min(4.*screenfac,screenmax)
        #limits of energy sampling
        epsil1 = .5-.5*np.sqrt(1.-screenmin/screenmax)
        epsilmin = max(epsil0,epsil1)
        epsilrange = .5 - epsilmin
        #sample energy rate of created electron
        F10 = ScreenFunction1(screenmin)-FZ
        F20 = ScreenFunction2(screenmin)-FZ
        NormF1 = max(F10*epsilrange*epsilrange,0.)
        NormF2 = max((1.5*F20),0.)
        while greject<np.random.rand():
            if NormF1/(NormF1+NormF2) > np.random.rand():
                epsil = .5 - epsilrange*(np.random.rand()**(1./3))
                screenvar = screenfac/(epsil*(1-epsil))
                greject = (ScreenFunction1(screenvar)-FZ)/F10
            else:
                epsil = epsilmin+epsilrange*np.random.rand()
                screenvar = screenfac/(epsil*(1-epsil))
                greject = (ScreenFunction2(screenvar)-FZ)/F20
    return epsil

def delta(epsil,epsil0,Z):
    return(136*epsil0)/(epsil*(1-epsil)*Z**(1./3))
def Phi1(d):
    result = 0
    if d<=1:
        result = 20.867 - 3.242*d + .625*d**2
    else:
        result = 21.12-4.184*np.log(d+.952)
    return result
def Phi2(d):
    result = 0
    if d<=1:
        result = 20.209 - 1.93*d + .086*d**2
    else:
        result = 21.12-4.184*np.log(d+.952)
    return result
#Code the Differential Cross Section from the required functions
#Differential cross section : PDF of positron energies for given photon energy
def PairCross(epsil,par):
    eps = epsil[0]
    Z = par[0]
    gammaE = par[1]
    epsil0 = .511/gammaE
    if eps<epsil0 or eps>1-epsil0:
        return 0
    else:
        FZ = F(Z,gammaE) #NOT a function of epsilon
        d = delta(eps,epsil0,Z)
        cr = (eps**2+(1-eps)**2)*(Phi1(d)-FZ/2.)+(2./3)*eps*(1-eps)*(Phi2(d)-FZ/2.)
        return cr
