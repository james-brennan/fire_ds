#-*- coding: utf-8 -*-
"""

TODO:
#1. Vary the noise across wavelength currently its the same
    it should be larger on NIR etc for vegetation


python generate_test_data.py

-- script which generates fires using modeled reflectance ..

-- also generates noise artifacts

-- can model brf effects using brf models...

-- can generate sampling effects etc...

"""

# for each location need to add the fire in
#now need
import sys
sys.setrecursionlimit(10000)
import matplotlib.pyplot as plt


# import stuff we need
#from fire import aFire
#from utils import make_noises, generate_spectra
#from spectral_model import spectral_fire_model, spectral_temporal_response
#import fire_testing
from dataset import dataset
#import dataset

def main():
    ds = dataset()
    ds.model_background_state()
    ds.model_fires(30)
    ds.model_cloud_cover()
    ds._save_to_gif()
    #ds.BA_to_tif()
    #ds.surface_refl_to_tif()

if __name__ == "__main__":
    main()
