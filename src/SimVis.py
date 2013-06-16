import Distributions
import numpy as np

c = 0.3 # m/ns
twopi = 2.*np.pi

class SimVis:
    def __init__(self, DistObject, BeamObject, bl, fqs, mfreq=0.15):
        self.prms = DistObject.SimParams()
        self.beam = BeamObject
        self.bl = bl
        self.fqs = fqs
        self.l2 = (c/fqs)**2

    def spindex(self,G):
        return (self.mfreq/self.fqs)**G

    def RMspec(self, flx, rm, pa, L, G):
        phs = (rm/np.pi)*self.l2 + (self.fqs/c)*(self.bl*_L) + pa
        return flx * self.spindex(G) * np.exp(-2.j*np.pi * phs)

    def SimVis(self):
        prms = self.prms
        vis = None
        for F,L,G,P,RM,X in zip(prms['F'],prms['L'],prms['G'],prms['P'],prms['RM'],prms['X']):
            if vis is None:
                vis = P * self.RMspec(F, RM, X, L, G)
            else:
                vis += P * self.RMspec(F, RM, X, L, G)
        return vis
