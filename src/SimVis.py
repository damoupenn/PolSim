import numpy as np
from scipy.special import erfc,erfcinv

c = 0.3 # m/ns

class Dist:
    def __init__(self,
                bm,
                aL, bL, aH, bH,
                Smin, So, Smax,
                PImu, PIsig, PImax)
        #observation definitions
        self.bm = bm

        #params for flux distriubtions
        self.aL = aL
        self.bL = bL
        self.aH = aH

        #flux limits and turnover
        self.Smin = Smin
        self.So = So
        self.Smax = Smax

        #pol fraction distribution parameters
        self.PImu = PImu
        self.PIsig = PIsig
        self.PImax = PImax

        #useful things to know...
        self.FindNsrc()

    def FindNsrc(self):
        dNdS_L = (self.aL/(1-self.bL))*(self.So**(1-self.bL))*(1.-(self.Smin/self.So)**(1-self.bL))
        dNdS_H = (self.aH/(1-self.bH))*(self.Smax**(1-self.bH))*(1.-(self.So/self.Smax)**(1-self.bH))
        self.Nsrc = int((dNdS_L+dNdS_H)*self.bm)

    def DrawSrcFlux(self):
        Xo = (self.bm/self.Nsrc)*(self.aL/(1.-self.bL))*(self.So**(1.-self.bL))*(1.-(self.Smin/self.So)**(1.-self.bL))
        def lo(x):
            Si = self.Smin*(1.+(self.Smin**(self.bL-1.))*((1.-self.bL)/self.aL)*(self.Nsrc/self.bm)*x)**(1./(1.-self.bL))
            return Si
        def hi(x):
            Si = self.So*(1.+(self.So**(self.bH-1.))*((1.-self.bH)/self.aH)*(self.Nsrc/self.bm)*(x-Fo))**(1./(1.-self.bH))
            return Si
        X = np.random.uniform(0,1,self.Nsrc)
        S = np.where(FS >= Xo, hi(FS), lo(FS))
        return S

    def DrawPolFrac(self):
        Ex = np.log(mu)
        Sx = np.sqrt(np.log(1.+(self.PIsig/self.PImu)**2))
        X = np.random.uniform(0,1,self.Nsrc)
        X *= 0.5*erfc(-1.*(np.log(self.PImax)-Ex)/(Sx*np.sqrt(2.)))
        X = Ex - np.sqrt(2.)*Sx*erfcinv(2.*X)
        return np.exp(X)
