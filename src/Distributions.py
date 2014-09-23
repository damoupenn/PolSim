import numpy as np
from scipy.special import erfc, erfcinv
from scipy.interpolate import interp1d

from pylab import *

c = 0.3 # m/ns
twopi = 2.*np.pi

class Dist:
    def __init__(self, config_file, beam, RMfile):
        self.bm = beam
        self.RMfile = RMfile
        self.ReadConfig(config_file)
        self.FindNsrc()

    def ReadConfig(self, config_file):
        for line in open(config_file).readlines():
            if line.startswith('#'):
                continue
            if line.startswith('aL'):
                self.aL = float(line.split(':')[1])
            if line.startswith('bL'):
                self.bL = float(line.split(':')[1])
            if line.startswith('aH'):
                self.aH = float(line.split(':')[1])
            if line.startswith('bH'):
                self.bH = float(line.split(':')[1])
            if line.startswith('Smin'):
                self.Smin = float(line.split(':')[1])
            if line.startswith('Smax'):
                self.Smax = float(line.split(':')[1])
            if line.startswith('So'):
                self.So = float(line.split(':')[1])
            if line.startswith('Pmean'):
                self.PImu = float(line.split(':')[1])
            if line.startswith('Psig'):
                self.PIsig = float(line.split(':')[1])
            if line.startswith('Pmax'):
                self.PImax = float(line.split(':')[1])
            if line.startswith('Gmean'):
                self.Gmean = float(line.split(':')[1])
            if line.startswith('Gsig'):
                self.Gsig = float(line.split(':')[1])


    def FindNsrc(self):
        _s = np.linspace(self.Smin,self.Smax,500)
        def hi(x): return self.bm * self.aH * x**-self.bH
        def lo(x): return self.bm * self.aL * x**-self.bL
        dNdS = np.where(_s >= self.So, hi(_s), lo(_s))
        NgtS = np.cumsum(dNdS)*(_s[1]-_s[0])
        self.NgtS = (NgtS[-1] - NgtS)
        self.Nsrc = int(self.NgtS[0] - self.NgtS[-1])
        self.NgtS /= self.NgtS[0]
        self.Fs = interp1d(self.NgtS,_s)


    def DrawSrcFlux(self):
        return self.Fs(np.random.uniform(self.NgtS.min(),self.NgtS.max(),self.Nsrc))

    def DrawLogNormal(self):
        if self.PImu == 0:
            return np.zeros(self.Nsrc)
        else:
            Ex = np.log(self.PImu)
            Sx = np.sqrt(np.log(1.+(self.PIsig/self.PImu)**2))
            X = np.random.uniform(0,1,self.Nsrc)
            X *= 0.5*erfc(-1.*(np.log(self.PImax)-Ex)/(Sx*np.sqrt(2.)))
            X = Ex - np.sqrt(2.)*Sx*erfcinv(2.*X)
            return np.exp(X)

    def DrawSphere(self):
        cosTh = np.random.uniform(0, 1, self.Nsrc)
        phi = np.random.uniform(0, twopi, self.Nsrc)
        Th = np.arccos(cosTh)
        l = np.sin(Th)*np.cos(phi)
        m = np.sin(Th)*np.sin(phi)
        return l, m

    def DrawFromCDF(self):
        _load = np.load(self.RMfile)
        _rms = _load['arr_0']
        _Frm = _load['arr_1']
        return np.interp(np.random.uniform(0, 1, self.Nsrc), _Frm, _rms)

    def SimSkyParams(self):
        prms = {}
        prms['P'] = self.DrawLogNormal()
        prms['G'] = self.Gsig * np.random.standard_normal(self.Nsrc) + self.Gmean
        prms['F'] = self.DrawSrcFlux()
        prms['L'], prms['M'] = self.DrawSphere()
        prms['RM'] = self.DrawFromCDF()
        prms['X'] = np.random.uniform(0.,1.,self.Nsrc)
        return prms
