import PolSim
import numpy as np
from pylab import *

config_file = 'config_files/6C.cf'

D = PolSim.Dist(config_file, 0.76, 'input_data/Oppermann.copy.npz')
prms = D.SimSkyParams()

SV = PolSim.SimVis(D, 'B', 0.15, np.linspace(0.12,0.18,100))

plot(SV.fqs, SV.SimVis().real)
show()
