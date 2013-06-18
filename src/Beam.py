import aipy as a
import numpy as np
import healpy as hp

class BeamError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class AipyBeam:
    def __init__(self,aa):
        self.aa = aa
    def Response(self, l, m, pol):
        z = np.sqrt(1.-l**2-m**2)
        aa.ants[0].set_active_pol(pol[0])
        BM = aa.ants[0].bm_response((L[i],M[i],z))[:,0]**2
        BM /= aa.ants[0].bm_response((0,0,1))[:,0]**2
        return BM

class HPXBeam:
    def __init__(self, hpx_alms):
        print 'Generating beams (this may take a second...'
        self.alms = hpx_alms
        self.bmx = hp.alm2map(hpx_alms)
        self.bmy = self.rotate_xy(self.bmx)
    def Response(self, l, m, pol):
        z = np.sqrt(1.-l**2-m**2)
        if pol == 'xx':
            return self.bmx
        elif pol == 'yy':
            return self.bmy
        px = hp.vec2pix(hp.get_nside(self.map), l, m, z)
        return self.map[px]
    def rotate_xy(self):
        xyz = hp.pix2vec(hp.get_nside(self.bmx), range(len(self.bmx)))
        X2Y = a.coord.rot_m(-np.pi/2., np.array([0,0,1]))
        ycrd = np.dot(X2Y, xyz)
        ypx = self.vec2pix(get_nside(self.bmx), ycrd[0], ycrd[1], ycrd[2])
        self.ybm = self.bmx[ypx]

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
                self.beam = HPXBeam(hpx_alms)
        else:
            raise(BeamError("Beam needs to come from aipy or from a healpix alm!"))
