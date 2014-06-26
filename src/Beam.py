import aipy as a
import numpy as np
import healpy as hp
from scipy.interpolate import interp1d

class BeamError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class AipyBeam:
    def __init__(self,aa):
        self.aa = aa
    def Response(self, l, m, fq, pol):
        z = np.sqrt(1.-l**2-m**2)
        self.aa.ants[0].set_active_pol(pol[0])
        BM = self.aa.ants[0].bm_response((l,m,z))[:,0]**2
        BM /= self.aa.ants[0].bm_response((0,0,1))[:,0]**2
        return BM

def get_freq_from_filename(fname):
    return float(fname.split('_')[-2][:3])*1e-3

class HPXBeam:
    def __init__(self, fitsdir, nside=64):
        print 'Generating beams (this may take a second...)'
        self.fitsdir = fitsdir
        self.nside = nside
        self.load_beam()
    def load_beam(self):
        self.afreqs = []
        self.bm_maps = {'xx': [], 'yy': []}
        for ff in sorted(__import__('glob').glob('%s/*.fits'%self.fitsdir)):
            self.afreqs.append(get_freq_from_filename(ff))
            _bm = hp.read_map(ff)
            _bm = hp.ud_grade(_bm, self.nside)
            try:
                _bm /= _bm[self.zen]
            except(AttributeError):
                self.zen = hp.vec2pix(self.nside, 0, 0, 1)
                _bm /= _bm[self.zen]
            self.bm_maps['xx'].append(_bm)
            if not hasattr(self, 'npx'):
                self.npx = len(_bm)
            self.bm_maps['yy'].append(self.rotate_xy(_bm))
        self.ndeg = len(self.bm_maps['xx'])
    def rotate_xy(self, m):
        xyz = hp.pix2vec(self.nside, range(self.npx))
        X2Y = a.coord.rot_m(-np.pi/2., np.array([0,0,1]))
        ycrd = np.dot(X2Y, xyz)
        ypx = hp.vec2pix(self.nside, ycrd[0], ycrd[1], ycrd[2])
        return m[ypx]
    def Response(self, l, m, fq, pol):
        z = np.sqrt(1.-l**2-m**2)
        px = hp.vec2pix(self.nside, l, m, z)
        f = interp1d(self.afreqs,
                    [self.bm_maps[pol][i][px] for i in range(self.ndeg)],
                    kind='cubic')
        return f(fq)

class Beam:
    def __init__(self, fromalm=False, fromaipy=False, hpx_alms=None, aa=None):
        if fromaipy:
            if aa is None:
                raise(BeamError("Pass an aipy AA object!"))
            else:
                self.beam = AipyBeam(aa)
        elif fromalm:
            if hpx_alms is None:
                raise(BeamError("Pass an ALM Pickle!"))
            else:
                self.beam = HPXBeam(hpx_alms)
        else:
            raise(BeamError("Beam needs to come from aipy or from a healpix alm!"))
