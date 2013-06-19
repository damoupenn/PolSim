from _pspec import *
from capo import pfb
import SimVis as SV

class SimVis(SV.SimVis):
    def __init__(self, DistObject, BeamObject, bl, fqs, mfreq=0.15):
        SimVis.__init__(self, DistObject, BeamObject, bl, fqs, mfreq=mfreq)

    def GenCell(self, f0, df):
        return None

    def Vis2k3Pk(self, Vis, Nann=12):

        Nchan = len(self.fqs)
        Vis *= _pspec.jy2T(self.fqs)

        returnD, returnK = {},{}
        for i in range(Nann):

            if i == 0:
                continue
            if i+1 == Nann:
                continue

            dchan = Nchan/Nann
            ch1 = i*dchan
            ch2 = (i+1)*dchan

            ch_in = ch1 - dchan
            ch_out = ch2 + dchan

            _fq = self.freqs[ch1:ch2]
            _v = Vis[ch_in:ch_out]
            pfreqs = self.fqs[ch_in:ch_out]

            BW = _fq[-1]-_fq[0]
            Trms, _ks = _pspec.Trms_vs_fq(pfreqs, _v, B=BW, window='blackman-harris')
            fmid = _ks.keys()[0]
            _d = _pspec.k3pk_from_Trms([Trms[fmid]],[1.],
                    k=_ks[fmid][0], fq=np.median(_fqs), B=BW)[0]

            zmid = _pspec.f2z(fmid)
            returnD[zmid] = _d
            returnK[zmid] = _ks[fmid][0]

        return returnK, returnD
