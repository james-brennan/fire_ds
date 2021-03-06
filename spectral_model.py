#-*- coding: utf-8 -*-
"""
    Spectral fire model
    ===================

    Model for spectral response following a fire

"""
import numpy as np

def spectral_temporal_response(num_days=20):
    """ Return of healthy vegetation length

        Computes the mixture of ash and recovered vegetation -->
        --> NOTE: obs very simple model of return of vegetation to health...
                -- no spectral change etc...
    """
    vegetation =  1.0 / (1.0 + np.exp(-np.linspace(-6,6,num_days)))
    char_ash = 1 - vegetation
    return vegetation, char_ash


def spectral_fire_model(surface_refl, timesteps, fire_i, ash_spectrum):
    """
    Parameters
    ----------
    surface_refl: array
        Array holding surface reflectance into which fire is modelled
    timesteps: int
        Number of timesteps
    fire_i: fire class
        Fire locations. Date of burn, x and y location
    ash_spectrum: Array
        spectral end member for ash
    """
    for i in xrange(len(fire_i.DOB)):
        #
        dob = int(fire_i.DOB[i])
        x =  fire_i.x[i]
        y =  fire_i.y[i]
        spectral_response_weights = spectral_temporal_response(num_days=50)
        # do mixture model
        if not (dob+50 > timesteps):
            surface_refl[dob:dob+50, :, x, y] = (
                            (spectral_response_weights[0]*surface_refl[dob:dob+50, :, x, y].T)+
                            spectral_response_weights[1] * ash_spectrum[:,np.newaxis]).T
    return surface_refl
