#! /usr/bin/env python

import numpy as np
import aipy as a
import PolSim

######################
# Startup Parameters #
######################

fstart = 0.12 #GHz
fstop = 0.18
Nchan = 1000
fqs = np.linspace(fstart, fstop, Nchan)
bl_EW = 16. #lambda
N_zbins = 12

CONFIG_FILE = 'config_files/6C.cf'
CAL_FILE = 'psa_null'
RM_DIST_FILE = 'input_data/Oppermann.copy.npz'

############################################
# Initialize Beam, SimVis and Dist Objects #
############################################

aa = a.cal.get_aa(CAL_FILE, fqs)
Bm = PolSim.Beam(fromaipy=True, aa=aa).beam
D = PolSim.Dist(CONFIG_FILE, 0.76, RM_DIST_FILE)
SV = PolSim.SimVis(D, Bm, bl_EW, fqs)

###########################
# Generate Power Spectra! #
###########################

visXX = SV.SimVis(pol = 'xx')
visYY = SV.SimVis(pol = 'yy')

K, PS = SV.Vis2k3Pk(visXX - visYY, Nann=N_zbins)
