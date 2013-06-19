import Distributions
import numpy as np

c = 0.3 # m/ns
twopi = 2.*np.pi

class SimVis:
    def __init__(self, DistObject, BeamObject, bl, fqs, mfreq=0.15):
        self.prms = DistObject.SimSkyParams()
        self.beam = BeamObject
        self.bl = bl
        self.fqs = fqs
        self.l2 = (c/fqs)**2
        self.mfreq = mfreq

    def spindex(self,G):
        return (self.mfreq/self.fqs)**G

    def RMspec(self, pf, flx, rm, pa, L, G):
        phs = (rm/np.pi)*self.l2 + self.fqs *(self.bl * L / c) + pa
        return pf * flx * self.spindex(G) * np.exp(-2.j*np.pi * phs)

    def SimVis(self, pol):
        prms = self.prms
        vis = None
        for F,L,M,G,P,RM,X in zip(prms['F'], prms['L'], prms['M'], prms['G'], prms['P'], prms['RM'], prms['X']):
            if vis is None:
                print self.beam
                vis = self.beam.Response(L, M, pol) * self.RMspec(P, F, RM, X, L, G)
            else:
                vis += self.beam.Response(L, M, pol) * self.RMspec(P, F, RM, X, L, G)
        return vis
