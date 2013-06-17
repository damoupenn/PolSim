import aipy as a
import numpy as np

class BeamError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class AipyBeam:
    def __init__(self,aa):
        self.aa = aa
    def Response(self, l, m, pol):
        z = np.sqrt(1.-l**2+m**2)
        aa.ants[0].set_active_pol(pol[0])
        BM = aa.ants[0].bm_response((L[i],M[i],z))[:,0]**2
        BM /= aa.ants[0].bm_response((0,0,1))[:,0]**2
        return BM

class Beam:
    def __init__(self, fromalm=False, fromaipy=False, hpx_alms=None, aa=None):
        if fromaipy:
            if aa is None:
                raise(BeamError("Pass an aipy AA object!"))
            else:
                self.beam = AipyBeam(aa)
        elif fromalm:
            if hpx_alms is None:
                raise(BeamError("Pass an array of healpix ALMs!"))
            else:
                self.beam = HPXBeam(aa)
        else:
            raise(BeamError("Beam needs to come from aipy or from a healpix alm!"))
