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

    def srcspec(self, flx, L, G):
        phs = self.fqs * (self.bl * L / c)
        return flx * self.spindex(G) * np.exp(-2.j*np.pi * phs)

    def RMphs(self, rm, pa):
        phs = (rm/np.pi)*self.l2 + pa
        return np.exp(-2.j*np.pi * phs)

    def SimVis(self, pol, remove=None):
        prms = self.prms
        vis = None
        if remove:
            print 'Removing %s/%d sources'%(remove,len(prms['F']))
            flim = sorted(prms['F'])[-(remove+1)]
            print '\t limiting flux = %.2f mJy'%(flim*1000.)
        for F,L,M,G,P,RM,X in zip(prms['F'], prms['L'], prms['M'], prms['G'], prms['P'], prms['RM'], prms['X']):
            if remove:
                if F >= flim:
                    continue
            sgn = {'xx':1., 'yy':-1.}[pol]
            I = self.srcspec(F, L, G)
            Q = sgn * P * self.RMphs(RM, X)
            if vis is None:
                vis  = self.beam.Response(L, M, self.fqs, pol) * (1.+Q) * I
            else:
                vis += self.beam.Response(L, M, self.fqs, pol) * (1.+Q) * I
        return vis
