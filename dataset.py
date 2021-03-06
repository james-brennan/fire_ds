#-*- coding: utf-8 -*-
"""
    Dataset creation class
    =======================



"""
import noise
from scipy.signal import convolve
import scipy.interpolate
import sys
import numpy as np
sys.setrecursionlimit(10000)
import matplotlib.pyplot as plt
import scipy.ndimage.filters


from fire import fire
from utils import make_noises, generate_spectra, make_landscape
from spectral_model import spectral_fire_model, spectral_temporal_response
import fire_testing


class dataset(object):
    def __init__(self, timesteps=100,
                 bands=13, sizex=100, sizey=100, spatial_res=1 ):
        self.xSize = sizex
        self.ySize = sizey
        self.timesteps = timesteps
        self.bands = bands
        self.spatial_res = spatial_res
        self.surface_rho = None
        self.modelled_fire = False
        self.fires = []

    @classmethod
    def generate(*args):
        """
        Method which produces the dataset following specifications on
        -- fire
        -- brdf effects
        -- cloud cover
        """

    def make_stack(self):
        self.surface_rho = np.ones((self.timesteps, self.bands,
                                    self.xSize, self.ySize))

    def model_brdf_effects(self, vza, sza):
        """
        Add BRDF effects by simulating the kernels over the image

        -- 1 need to generate slight variations in vza and sza
            across the image... hmmm cehck the literature?
        """
        pass

    def model_cloud_cover(self):
        """
        Add clouds to remove some pixels

        run after fires has been modelled into scene...

        """
        """

        idea for now is to re-use perlin noise generator
        With a threshold for how many clouds there should be:
            -- could be semi-seasonal + random element

        """
        # model of cloud cover across timesteps
        #       maximum seasonal component 20%          random weather
        cloud = (0.2 * np.sin(np.linspace(0,np.pi,self.timesteps)) +
                    np.random.normal(0.2, 0.1, self.timesteps))
        # make sure max is one
        cloud[cloud > 1] = 1
        # and min is 0
        cloud[cloud < 0] = 0
        # generate some new perlin noise
        the_noise = make_noises(self.timesteps, self.xSize, time_multiple=200)
        """
         idea is to threshold the perlin noise by the percentile of cloud_cover
         this produces a binary mask that should produce realistic looking
         clouds
         """
        percentiles = np.array([np.percentile(no, cl*100) for no, cl
                                            in zip(the_noise, cloud)])
        # now threshold the noise below these values
        for day in xrange(self.timesteps):
            # mask in cloud
            the_noise[day][the_noise[day] < percentiles[day]] = True # is cloud
            # make into a binary
            the_noise[day][the_noise[day] != True] = False
        # convert all to a boolean
        the_noise = the_noise.astype(np.bool)
        # apply cloud mask to data
        self.cloud_mask = ~the_noise
        self.surface_rho = (self.surface_rho[:].T * self.cloud_mask[:,None].T).T
        # make it into a masked array for clarity
        #import pdb; pdb.set_trace()
        self.surface_rho = np.ma.array(self.surface_rho,
                                        mask=self.surface_rho==0)
        return 0


    def model_background_state(self):
        """
        Make some sort of realistic dynamics for the background
        state of reflectance
        """
        if self.surface_rho == None:
            self.make_stack()
        # load spectra
        vegetation_rho, ash_rho = generate_spectra()
        # select some bands...
        vegetation_rho = vegetation_rho[::25]
        ash_rho = ash_rho[::25]
        # keep these for later...
        self.ash_rho = ash_rho
        self.vegetation_rho = vegetation_rho

        # first fill the dataset with just normal leaf reflectance across all timesteps
        self.surface_rho = (self.surface_rho[:].T*vegetation_rho[None, :].T).T

        # now make a mean state...
        mean_state = make_noises(1, self.xSize)
        mean_state *= 0.3
        mean_state = np.tile(mean_state, (self.bands,self.timesteps,1,1),)
        mean_state = np.swapaxes(mean_state, 0,1)

        # alternative with more spatial features -- less like a cloud
        land_cover = make_landscape(self.xSize)

        # for each catergory in the land_cover assign a mean variation on
        # the spectra
        lc_spe_multiple = np.linspace(0.85, 1.15, len(np.unique(land_cover)))
        # associate these values with land_cover
        for i,j in enumerate(np.unique(land_cover)):
            land_cover[np.where(land_cover==j)] = lc_spe_multiple[i]

        #import pdb; pdb.set_trace()
        # multiply this with spectra
        lc_effect = np.tile(land_cover, (self.timesteps,self.bands,1,1),)
        #mean_state = np.swapaxes(mean_state, 0,1)
        self.surface_rho *= lc_effect

        return None
        # now generate temporal variability
        temporal = make_noises(self.timesteps, self.xSize)
        temporal *= 0.2

        # now add to the surface reflectance
        # need same dimensions
        temporal = np.tile(temporal, (self.bands,1,1,1),)
        # switch dimensions
        temporal = np.swapaxes(temporal, 0,1)
        self.surface_rho += temporal

    def model_fires(self, burned_pixels_percent=0):
        """

        if burned_pixels_percent > 0:
            inversion based approach to model
            a percentage of the pixels
            as burned. Uses scipy.optimize --
                NOTE: is quite slow...
        """
        # now generate some fires
        self.fire_locs = []
        if burned_pixels_percent == 0:
            # just model some fires...
            import pdb; pdb.set_trace()
            self.fires = burnIt_idea1(self, seeds=1, decay=5)
        elif  burned_pixels_percent > 0:
            # invert simulations to match with a pixel %...
            # 1. calculate total pixels...
            self.pixels = self.xSize * self.ySize
            # how many burnt?
            to_burn = ((burned_pixels_percent)/100.0) * self.pixels
            seeds = 5.0
            # per seed...
            to_burn /= seeds
            # call the fire routine
            # OLD ROUTINE
            #self.fires = burnIt_idea1(self, seeds=5, numPixels=to_burn)
            # try new method
            numFires=40
            i = 0
            while i <= numFires:
                # choose a date at random after initial possible
                # burn date
                spark_date = 15 + np.random.randint(0,15)
                # set size of fire
                fire_size = np.random.exponential(100)
                try:
                    f = fire_testing.make_fire(size=self.xSize, max_size=fire_size)
                    self.fires.append( f )
                    i+=1
                except fire_testing.FireIssue:
                    print  """ algorithm failed """
                    continue
                    # don't increment i count
                #import pdb; pdb.set_trace()
                #f = np.array(f)
                #if len(f) ==0:
                #    continue
                #burn_map = np.zeros((size, size))
                #burn_map[f.T[1].astype(int), f.T[2].astype(int)] = f.T[0]
                 #np.where(self.fires)
                #import pdb; pdb.set_trace()
                #self.fire_locs = np.array(self.fire_locs)
                #import pdb; pdb.set_trace()
                #import pdb; pdb.set_trace()
                # Now run mixture model for spectral response to fire
                #dob = 10
                #import pdb; pdb.set_trace()
            for fire_i in self.fires:
                self.surface_rho = spectral_fire_model(self.surface_rho, self.timesteps,
                                                    fire_i, self.ash_rho)

    def surface_refl_to_tif(self):
        """
        save the dataset to a file
        save the truth BA data as well...
        """
        import gdal
        filenames = ["s_refl.%03i.tif" % doy for doy in xrange(self.timesteps)]
        # first output
        outdriver = gdal.GetDriverByName("GTiff")
        # loop over dates...
        for doy in xrange(self.timesteps):
            outdata = outdriver.Create(str(filenames[doy]),
                                     self.xSize, self.ySize,
                                     self.bands, gdal.GDT_Float32,
                                     ['COMPRESS=LZW'])
            # write out each band of the surface refl
            for band in xrange(self.surface_rho[doy].shape[0]):
                outdata.GetRasterBand(band+1).WriteArray(self.surface_rho[doy][band])
        return None

    def BA_to_tif(self):
        """
        save the truth BA data
        """
        import gdal
        filenames = ["BA.%03i.tif" % doy for doy in xrange(self.timesteps)]
        # first output
        outdriver = gdal.GetDriverByName("GTiff")
        # loop over dates...
        for doy in xrange(self.timesteps):
            outdata   = outdriver.Create(str(filenames[doy]),
                                     self.xSize, self.ySize,
                                     1, gdal.GDT_Byte, ['COMPRESS=LZW'])
            outdata.GetRasterBand(1).WriteArray(self.fires[doy])
        return None

    def _save_to_gif(self):
        """
        make an animated gif of the surface refl
        """
        import matplotlib.pyplot as plt
        filenames = []
        for day in xrange(self.timesteps):
            #plt.figure()
            #plt.subplot(121)
            plt.imshow(self.surface_rho[day, 10], interpolation='nearest', cmap='Greys_r')
            plt.colorbar()
            fname = "rho_%03i.png" % day
            plt.title(fname)
            #plt.subplot(122)
            # plot by burndate
            #import pdb; pdb.set_trace()
            ##plt.imshow(self.fires[day], interpolation='nearest', cmap='Greys_r', vmin=0, vmax=100)
            #plt.colorbar()
            filenames.append(fname)
            plt.tight_layout()
            plt.savefig(fname)
            plt.close()
        # also run terminal command to make gif...
        import os
        os.system('convert -delay 20 -loop 0 *.png animation.gif')
