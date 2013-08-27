#! /usr/bin/env python

from PolSim import *
import numpy as np
from scipy.interpolate import interp1d
import sys

BmPickle = sys.argv[1]

BM = Beam(fromalm=True, hpx_alms=BmPickle)

from pylab import *
fq = np.linspace(0.05,0.25,100)
plot(fq, BM.beam.Response(0,0,fq,'xx'))
plot(fq, BM.beam.Response(0,0,fq,'yy'))
show()
