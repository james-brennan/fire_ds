#-*- coding: utf-8 -*-
"""
    Utility Functions
    =================

    Functions needed for various tasks

"""
import noise
from scipy.signal import convolve
import scipy.interpolate
import numpy as np
import scipy.ndimage.filters

def generate_spectra():
    veg_refl = np.genfromtxt('leafReflTrans')
    wv = veg_refl.T[0]
    veg_refl = veg_refl.T[1]

    # smooth out the higher frequency noise -- this is an idealised example...
    window = np.ones(3)/3.0
    veg_refl  = convolve(veg_refl, window, mode='same')
    # load a soil profile for the moment
    soil_ash_char = np.genfromtxt('soilRefl')
    # interpolate soil_ash_char to same wavelengths as vegetation
    wv1d = scipy.interpolate.interp1d(soil_ash_char.T[0], soil_ash_char.T[1])
    soil_ash_char2 = wv1d(wv)
    soil_ash_char2  = convolve(soil_ash_char2, window, mode='same')
    # decrease this a bit -- only using this for testing
    soil_ash_char2 *= 0.5
    #
    ## ok now try and make 1d example where fire happens on a date and recovers
    ## simple mixture model basically...
    return veg_refl, soil_ash_char2


#seems to only work in a loop atm!!!
# do same noise across bands but varying in space and time -- use 3d noise for now

def make_noises(timesteps, size, time_multiple=1):
    noises = np.zeros((timesteps, size,size))
    for t,z in enumerate(np.linspace(0,1,timesteps)):
        #print t
        for i,x in enumerate(np.linspace(0,1,size)):
            for j,y in enumerate(np.linspace(0,1,size)):
                #print x,y
                # increasign the number of octaves increases the noise...
                noises[t, j,i] = noise.snoise3(x,y,time_multiple*z, octaves=2) # not sure what octave to pick...
    return noises
