import numpy as np
from scipy.special import erfc,erfcinv

c = 0.3 # m/ns
twopi = 2.*np.pi

class Dist:
    def __init__(self, config_file, beam):
        self.bm = beam
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

    def DrawSphere(self):
        cosTh = np.random.uniform(0, 1, self.Nsrc)
        phi = np.random.uniform(0, twopi, self.Nsrc)
        Th = np.arccos(cosTh)
        l = np.sin(Th)*np.cos(phi)
        m = np.sin(Th)*np.sin(phi)
        return l, m

    def DrawFromCDF(self, fname):
        _load = np.load(fname)
        _rms = _load['arr_0']
        _Frm = _load['arr_1']
        return np.interp(np.random.uniform(0, 1, self.Nsrc), _Frm, _rms)

    def SimSkyParams(self, RmFname):
        prms = {}
        prms{'P'} = self.DrawLogNormal()
        prms{'G'} = self.Gsig * np.random.standard_normal(self.Nsrc) + self.Gmean
        prms{'F'} = self.DrawSrcFlux()
        prms{'L'}, prms{'M'} = self.DrawSphere()
        prms{'RM'} = DrawFromCDF(RmFname)
        return prms
