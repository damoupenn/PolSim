#! /usr/bin/env python

import numpy as np
import aipy as a 
import sys,os
from aipy.coord import *
from healpy import *
from pylab import *

Nside = 64
Pfrac = 0.001
Ntimes = 50

fq = np.linspace(0.12,0.18,1000)

r = Rotator(coord=['E','G'])
def hpm_gal2eq(N,_r):
    _nside = npix2nside(N)
    _px = _r(pix2vec(_nside,range(N)))
    return vec2pix(_nside,_px[0],_px[1],_px[2])

#Generate map of RM values
RMmap = read_map('sim_scripts/input_data/faraday.fits',hdu=3)
RMmap = alm2map(map2alm(RMmap,lmax=int(2.*np.pi*16)),get_nside(RMmap))
RMmap = RMmap[hpm_gal2eq(len(RMmap),r)]
RMmap = ud_grade(RMmap,Nside)

#Generate de Oliviera-Costa map at 150 MHz
Tgal = read_map('sim_scripts/input_data/dOC_gsm_150MHz.fits')
Tgal = read_map('/home/damo/coverage/lambda_haslam408_dsds.fits')*(408./150.)**2.76
Tgal = Tgal[hpm_gal2eq(len(Tgal),r)]
#Tgal *= Pfrac
Tgal = ud_grade(Tgal,Nside)

#Calculate polarization angle, aligned with the Galaxy.
X = np.load('sim_scripts/input_data/PolAngWmapK.npz')['arr_0']
X = X[hpm_gal2eq(len(X),r)]
X = ud_grade(X,Nside)

#Generate Beam map

from bm_prms import prms
bm = prms['beam'](np.array([0.15]),nside=Nside,lmax=20,mmax=20,deg=7)
bm.set_params(prms['bm_prms'])
px = range(nside2npix(Nside))
xyz = pix2vec(Nside,px)
BM = np.array([h.map[px] for h in bm.hmap])
    #Ax = np.polyval(poly,fq)
    #Ax = np.where(xyz[-1] > 0.001, Ax, 0.)
    #BM['x'].append(Ax**2)
    #ycrd = np.dot(X2Y,xyz)
    #ypx = vec2pix(get_nside(Ax),ycrd[0],ycrd[1],ycrd[2])
    #BM['y'].append(Ax[ypx])
X2Y = rot_m(np.pi/-2.,np.array([0,0,1]))
ycrd = np.dot(X2Y,xyz)
ypx = vec2pix(Nside,ycrd[0],ycrd[1],ycrd[2])

def rot_bm(hpm,ha,dec):
    _m = eq2top_m(-ha,dec)
    _xyz = pix2vec(get_nside(hpm),range(len(hpm)))
    ncrd = np.dot(_m,_xyz)
    n_px = vec2pix(get_nside(hpm),ncrd[0],ncrd[1],ncrd[2])
    return hpm[n_px]

#Setup coordinates. 
aa = a.cal.get_aa ('psa_null',np.array([fq]))
U_ns = np.array([32.,0.,0.])/0.3
Alt,Az = pix2ang(Nside,range(len(RMmap)))
Alt = np.pi/2. - Alt

x,y = np.cos(Alt)*np.cos(Az),np.cos(Alt)*np.sin(Az)
z = np.sqrt(1.-x**2 - y**2)

V = {}
for t in np.linspace(0,a.const.sidereal_day/a.const.s_per_day,Ntimes):
    aa.set_jultime(t+0.5)
    Ha = aa.sidereal_time()
    print Ha
    if not Ha in V.keys(): V[Ha] = {}
    for p in 'xy':
        if not p in V[Ha].keys(): V[Ha][p] = []
        for _fq in fq:
            _bm = np.polyval(BM,_fq)
            _bm = np.where(xyz[-1] > 0., _bm, 0.)
            if p == 'y': _bm = _bm[ypx]
            _bm = _bm**2
            Interferometry = _bm * np.exp(-2.j*np.pi * _fq *(U_ns[0]*x + U_ns[1]*y))
            lam2 = (0.3/_fq)**2
            if p == 'x':
                Integrand = rot_bm(Interferometry,Ha,aa.lat) * Tgal * (_fq/0.15)**-2.7 * ( 1. + Pfrac * np.exp(1.j*(X + 2.*RMmap*lam2)))
            elif p == 'y':
                Integrand = rot_bm(Interferometry,Ha,aa.lat) * Tgal * (_fq/0.15)**-2.7 * ( 1. - Pfrac * np.exp(1.j*(X + 2.*RMmap*lam2)))
            V[Ha][p].append(np.sum(Integrand))
        V[Ha][p] = np.array(V[Ha][p])
for t in V.keys():
    for p in V[t].keys():
        name = 'CorrSimData/CorrelatedVis%s_%s.npz'%(t,p)
        np.savez(name,V[t][p])
show()
