import Distributions
import numpy as np

c = 0.3 # m/ns
twopi = 2.*np.pi

class SimVis:
    def __init__(self, DistObject, BeamObject, bl, fqs):
        self.prms = DistObject.SimParams()
        self.beam = BeamObject
        self.bl = bl
        self.fqs = fqs
        self.l2 = (c/fqs)**2

    def RMspec(self, flx, nu, rm, pa, s, bl):
        phs = (rm/np.pi)*self.l2 + (nu/c)*(bl*s[0]) + pa
        return flx * np.exp(-2.j*np.pi * phs)
