import polsim

config_file = 'config_files/6C.cf'

D = PolSim.Distribution(config_file, 0.76)
prms = D.SimSkyParams('input_data/Oppermann.copy.npz')
