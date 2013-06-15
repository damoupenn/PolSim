import numpy as np
import aipy as a
from capo import pspec,pfb
import RotMeasTools as RMT
from scipy.special import erfc,erfcinv
from pylab import *
import sys
import healpy
C_ell = False

dists = []
#dists.append('testI')
#dists.append('VLSS')
#dists.append('6C')
#dists.append('hiRM')
#dists.append('hiN')
#dists.append('RMsrc_1k')
#dists.append('RMsrc_2k')
#dists.append('RMsrc_3k')
#dists.append('RMsrc_5k')
#dists.append('RMsrc_8k')
#dists.append('RMsrc_10k')
#dists.append('RMbogeys')
#dists.append('Correlate')
dists.append('agree')

c = 0.3 #m/ns

fstart,fstop = 0.12,0.18 #GHz
Nchan = 1000
freqs = np.linspace(fstart,fstop,Nchan)
dnu = freqs[1]-freqs[0]
bm_poly = pspec.DEFAULT_BEAM_POLY #sr
bl = 16. #lambda
window = 'blackman-harris'

Niter = 25
Pmax = 0.3

Nann = 12

def FindNsrc(Smin,Smax,So,aL,bL,aH,bH,bm=0.76):
    dNdO_L = (aL/(1-bL))*(So**(1-bL))*(1.-(Smin/So)**(1-bL))
    dNdO_H = (aH/(1-bH))*(Smax**(1-bH))*(1.-(So/Smax)**(1-bH))
    N = int((dNdO_L+dNdO_H)*bm)
    return N

def DrawSrcFlx(Smin,Smax,So,aL,bL,aH,bH,bm=0.76):
    Nsrc = FindNsrc(Smin,Smax,So,aL,bL,aH,bH,bm)
    FS = np.random.uniform(0.,1.,Nsrc)
    Fo = (bm/Nsrc)*(aL/(1-bL))*(So**(1-bL))*(1.-(Smin/So)**(1-bL))
    def lo(x):
        Si = Smin*(1.+(Smin**(bL-1.))*((1-bL)/aL)*(Nsrc/bm)*x)**(1./(1.-bL))
        return Si
    def hi(x):
        Si = So*(1.+(So**(bH-1.))*((1.-bH)/aH)*(Nsrc/bm)*(x-Fo))**(1./(1.-bH))
        return Si
    S = np.where(FS >= Fo, hi(FS), lo(FS))
    return S

def DrawLogNormal(N,mu,sig,Xmax=0.1):
    Ex = np.log(mu)
    Vx = np.log(1.+(sig/mu)**2)
    mu,sig = Ex,np.sqrt(Vx)
    X = np.random.uniform(0.,1.,N)
    #set an upper limit of Xmax
    X *= 0.5*erfc(-1.*(np.log(Xmax)-mu)/(sig*np.sqrt(2.)))
    X = mu - np.sqrt(2.)*sig*erfcinv(2.*X)
    X = np.exp(X)
    return X

def DrawSphere(Nsrc):
    mu = np.random.uniform(0,1,Nsrc)
    phi = np.random.uniform(0,2.*np.pi,Nsrc)
    theta = np.arccos(mu)
    l = np.sin(theta)*np.cos(phi)
    m = np.sin(theta)*np.sin(phi)
    return l,m

def DrawFromCDF(N):
    _load = np.load('sim_scripts/Oppermann.copy.npz')
    _rms = _load['arr_0']
    _Frm = _load['arr_1']
    return np.interp(np.random.uniform(0,1,N),_Frm,_rms)

def SimSkyParams(Smin,Smax,So,aL,bL,aH,bH):
    Nsrc = FindNsrc(Smin,Smax,So,aL,bL,aH,bH)
    if Niter <= 25: print Nsrc
    #Polarizaiton Fraction
    if not dist == 'agree':
        Pmean = 0.0201
        Psig = 0.0384
    else:
        Pmean = 0.001
        Psig  = 0.001
    P = DrawLogNormal(Nsrc,Pmean,Psig)
    #Spectral Index
    GAMmean = 0.8
    GAMsig = 0.1
    GAM = GAMsig * np.random.standard_normal(Nsrc)+GAMmean
    #flux
    flx = DrawSrcFlx(Smin,Smax,So,aL,bL,aH,bH)
    #Position, uniform in l,m
    L,M = DrawSphere(Nsrc)
    if dist == 'Correlate':
        OpFile = 'sim_scripts/input_data/faraday.fits'
        RMmap = healpy.read_map(OpFile,hdu=3)
        Z = np.sqrt(1-L**2-M**2)
        px = healpy.vec2pix(healpy.get_nside(RMmap),L,M,Z)
        RM = RMmap[px]
        X = np.zeros(Nsrc)
    else:
        RM = DrawFromCDF(Nsrc)
        if dist == 'hiRM': RM *= 2.
        X= np.random.uniform(0.,1.,Nsrc)

    return flx,L,M,P,GAM,RM,X

def RMspec(flx,nu,rm,pa,s,bl):
    l2 = (0.3/nu)**2
    phs = (rm/(np.pi))*l2 + (nu/c)*(bl*s[0]) + pa
    return flx * np.exp(-2.j*np.pi * phs)

#Calculate the beam for x and y pol.
aa = a.cal.get_aa('psa_null',freqs)

for dist in dists:
    if 'RMsrc' in dist:
        Fmax = []
        Nrm = 1000*int((dist.split('_')[-1]).split('k')[0])
    if C_ell and dist != '6C': continue

    print '#'*30
    print 'Working on %s' %dist
    print '#'*30
    if dist in ('6C','hiRM','Oppermann','Correlate','agree') \
        or dist.startswith('hiN')   \
            or dist.startswith('RMsrc'):
        if dist.startswith('hiN'): Smin = 0.1
        else: Smin = 0.16
        Smax = 10.
        So=0.88
        aL,bL = 4000.*0.88**-0.76,1.75
        aH,bH = 4000.,2.81
    if dist in ('VLSS'):
        Smin = 0.5
        Smax = 100.
        So=1.
        aL,bL = 4865,2.3
        aH,bH = 4865,2.3
    if dist == 'testI':
        Smin=0.1
        Smax=100.
        So=1.
        aL,bL = 1000.,2.
        aH,bH = 1000.,2.

    Nsrc = FindNsrc(Smin,Smax,So,aL,bL,aH,bH)
    Pspec = {'Q':{},'I':{}}
    ks = {}
    print '%d Sources'%Nsrc

    if C_ell:
        POS,FLX,POL = [],[],[]
        filter = np.where(freqs <= 0.1503 - 0.0025/2., 0., 1.,)
        filter = np.where(freqs >= 0.1503 + 0.0025/2., 0., filter)
    for iter in range(Niter):
        print 'Working on iteration %d/%d' % (iter+1,Niter)
        flx,L,M,P,GAM,RM,X = SimSkyParams(Smin,Smax,So,aL,bL,aH,bH)
        if dist=='VLSS': flx *= (0.074/0.150)**0.79
        if 'RMsrc' in dist:
            inds = np.argsort(flx*P)[::-1]
            flx = np.take(flx,inds)
            L = np.take(L,inds)
            M = np.take(M,inds)
            P = np.take(P,inds)
            GAM = np.take(GAM,inds)
            RM = np.take(RM,inds)
            X = np.take(X,inds)
            number_of_rm_srcs = 0
        Vnu = {'Q':np.zeros(Nchan,dtype=np.complex),'I':np.zeros(Nchan,dtype=np.complex)}
        for i in range(Nsrc):
            #if dist == 'RMsrc' and P[i]*flx[i] >= Pmax: continue
            if 'RMsrc' in dist and i <= Nrm:
                _fmax = P[i]*flx[i]
                continue
            if dist == 'RMbogeys' and number_of_rm_srcs <= Nrm:
                if 10. >= np.abs(RM[i]):
                    _fmax = P[i]*flx[i]
                    number_of_rm_srcs += 1
                    continue
            z = np.sqrt(L[i]**2+M[i]**2)
            aa.ants[0].set_active_pol('x')
            BMx = aa.ants[0].bm_response((L[i],M[i],z))[:,0]**2
            BMx /= aa.ants[0].bm_response((0,0,1))[:,0]**2
            aa.ants[0].set_active_pol('y')
            BMy = aa.ants[0].bm_response((L[i],M[i],z))[:,0]**2
            BMy /= aa.ants[0].bm_response((0,0,1))[:,0]**2
            PF = flx[i]*(0.151/freqs)**(GAM[i])
            srcspec = RMspec(PF,freqs,RM[i],X[i],(L[i],M[i]),bl)
            if C_ell:
                _f = filter*srcspec*pspec.jy2T(freqs,bm_poly=bm_poly)
                FLX.append(np.abs(_f).sum()/filter.sum())
                POS.append((L[i],M[i]))
                POL.append(P[i])
            else:
                if dist=='testI':
                    Vnu['Q'] += (BMx+BMy)*P[i]*srcspec
                    Vnu['I'] += (BMx+BMy)*PF
                else:
                    Vnu['Q'] += (BMx+BMy)*srcspec*P[i] + (BMx-BMy)*PF
                    Vnu['I'] += (BMx-BMy)*srcspec*P[i] + (BMx+BMy)*PF

        if C_ell:
            nu0 = (np.sum(freqs*filter)/np.sum(filter))*1000.
            bw  = 6.0
            savename = '/home/damo/papers/faraday_leakage_2012/fig_data/PolSim_fq%3.1fMHz_bw%1.1fMHz_iter%d.npz'%(nu0,bw,iter)
            print 'Saving to: %s' % savename
            np.savez(savename,FLX=FLX,POS=POS,POL=POL)

        if iter == 0:
            np.savez('/home/damo/papers/faraday_leakage_2012/fig_data/params_sampled_'+dist+'.npz',
                RM=RM,flx=flx,Pflx=flx*P)
        for pol in Vnu:
            for i in range(Nann):
                if i == 0: continue
                if i+1 == Nann: continue
                ch1,ch2 = i*(Nchan/Nann),(i+1)*(Nchan/Nann)
                dchan = ch2-ch1
                ch_in,ch_out = ch1-dchan,ch2+dchan
                _freqs = freqs[ch1:ch2]
                Vnu_i = Vnu[pol][ch_in:ch_out]
                pfreqs = freqs[ch_in:ch_out]
                BW = _freqs[-1]-_freqs[0]
                Trms,_ks = pspec.Trms_vs_fq(pfreqs,Vnu_i, B=BW, window='blackman-harris')
                fmid = _ks.keys()[0]
                _d = pspec.k3pk_from_Trms([Trms[fmid]],[1.],k=_ks[fmid][0],fq=np.median(_freqs),B=BW)[0]
                zmid = str(pspec.f2z(fmid))
                if not zmid in Pspec[pol].keys():
                    Pspec[pol][zmid] = []
                    ks[zmid] = _ks[fmid][0]
                Pspec[pol][zmid].append(_d)
        if 'RMsrc' in dist:
            Fmax.append(_fmax)
            print '='*30
            print np.mean(Fmax)
            print '='*30

    if C_ell: sys.exit(0)
    np.savez('/home/damo/papers/faraday_leakage_2012/fig_data/ks_'+dist+'.npz',**ks)
    np.savez('/home/damo/papers/faraday_leakage_2012/fig_data/Pspec_Q_'+dist+'.npz',**Pspec['Q'])
    np.savez('/home/damo/papers/faraday_leakage_2012/fig_data/Pspec_I_'+dist+'.npz',**Pspec['I'])
