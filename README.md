PolSim
======

This module generates power spectra as in [my paper](http://arxiv.org/abs/1302.0876). 

## Use It!

Each file in `src` contains a class described by that filename:

- `Distributions.Dist` handles choosing source parameters  
- `SimVis.SimVis` handles simulating a visibility.
- `Pspec.SimVis` allows you to generate a power spectrum from a simulated visibility.

Observation parameters (frequencies, etc) are passed as arguments to SimVis. Beam can either be generated from an Aipy AntennaArray object or from a series of healpix alm object (not currently implemented).

Source parameters are specified in a .cf file. 

## Examples!!!

- Examples of how to use each object is in `scripts/SamplePSpec.py`.
- An example .cf file can be found in `config_files/6C.cf` These are the parameters I used in my paper.
