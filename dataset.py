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


from fire import aFire
from utils import make_noises, generate_spectra
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
        self.fires = None

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
        #       maximum seasonal component 80%          random weather
        cloud = (0.2 * np.sin(np.linspace(0,np.pi,self.timesteps)) +
                    np.random.normal(0.2, 0.1, self.timesteps))
        # make sure max is one
        cloud[cloud > 1] = 1
        # generate some new perlin noise
        the_noise = make_noises(self.timesteps, self.xSize, time_multiple=20)
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


    def model_fires(self, burned_pixels_percent=0):
        """

        if burned_pixels_percent > 0:
            inversion based approach to model
            a percentage of the pixels
            as burned. Uses scipy.optimize --
                NOTE: is quite slow...
        """
        if self.surface_rho == None:
            self.make_stack()
        # load spectra
        vegetation_rho, ash_rho = generate_spectra()
        # select some bands...
        vegetation_rho = vegetation_rho[::25]
        ash_rho = ash_rho[::25]

        # first fill the dataset with just normal leaf reflectance across all timesteps
        self.surface_rho = (self.surface_rho[:].T*vegetation_rho[None, :].T).T
        # now generate the noise
        noises = make_noises(self.timesteps, self.xSize)
        noises *= 0.1

        # now add to the surface reflectance
        #import pdb; pdb.set_trace()
        # need same dimensions
        noises = np.tile(noises, (self.bands,1,1,1),)
        # switch dimensions
        noises = np.swapaxes(noises, 0,1)
        self.surface_rho += noises

        # now generate some fires
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
            size=100
            f = fire_testing.fire(size=size)
            #import pdb; pdb.set_trace()
            f = np.array(f)
            #burn_map = np.zeros((size, size))
            #burn_map[f.T[1].astype(int), f.T[2].astype(int)] = f.T[0]
            self.fire_locs = f.T #np.where(self.fires)


        self.fire_locs = f.T #np.where(self.fires)
        #import pdb; pdb.set_trace()
        # Now run mixture model for spectral response to fire
        #dob = 10
        self.surface_rho = spectral_fire_model(self.surface_rho, self.timesteps,
                                                self.fire_locs, ash_rho)

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
            outdata   = outdriver.Create(str(filenames[doy]),
                                     self.xSize, self.ySize,
                                     self.bands, gdal.GDT_Float32,
                                     ['COMPRESS=LZW'])
            # write out each band of the surface refl
            for band in xrange(self.surface_rho[doy].shape[0]):
                outdata.GetRasterBand(band+1).WriteArray(self.surface_rho[doy][band])
        return None

    def _save_to_gif(self):
        """
        make an animated gif of the surface refl
        """
        import matplotlib.pyplot as plt
        filenames = []
        for day in xrange(self.timesteps):
            plt.imshow(self.surface_rho[day, 10], interpolation='nearest', cmap='jet')
            plt.colorbar()
            fname = "rho_%03i.png" % day
            plt.title(fname)
            filenames.append(fname)
            plt.tight_layout()
            plt.savefig(fname)
            plt.close()
        # also run terminal command to make gif...
        import os
        os.system('convert -delay 20 -loop 0 *.png animation.gif')

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


# lets just start of simple and generate a fire in the middle on a choosen dataset

def burnIt_idea1(dataset, seeds=1, numPixels=1000, temporal=True):
    """
    Generate burned areas that occur around the date of burn choosen
    easy way is to create a 3d boolean mask where only TRUE are where
    burn is (on DOB)


    This methods starts from a number of seeds and grows outwards randomly
    to produce burned areas...
    --- also has a temporal flag so that the burn spreads in time..
    """
    size=dataset.xSize
    timesteps = dataset.timesteps
    DOB = 10
    burns = []
    for i in xrange(seeds):
        burns.append((DOB, np.random.uniform(0, size, 2).astype(int)))
    # now for each of the seeds
    # run the fire spread algorithm
    fires = []
    for f in xrange(seeds):
        this_fire = aFire(burns[f][0], burns[f][1][0], burns[f][1][1], numPixels=numPixels)
        fires.append(this_fire)
    # now put the fires into a boolean mask...
    bools = np.zeros((timesteps, size, size)).astype(np.int)
    for fire in fires:
        # make sure placing in right place!
        for burnday in zip(fire.DOB, fire.x, fire.y):
            day = burnday[0]
            x = burnday[1]
            y = burnday[2]
            try:
                bools[day,x,y] = day
            except:
                pass
    # so bools provides the DOB for each pixel now...
    """
        NOTE: FOR NOW ALGORITHM REVISITS CELLS WHICH HAVE BEEN
        BURNT ALREADY. THIS IS UNREALISTIC. SO FOR NOW JUST TAKE FIRST
        OCCURENCE OF FIRE IN EACH PIXEL THROUGH TIME.
    """
    for x in xrange(size):
        for y in xrange(size):
            # take time stack...
            timestack = bools[:, x,y].copy()
            # if more than one fire occurs...
            idx = np.nonzero(timestack)
            #
            #import pdb; pdb.set_trace()
            if idx[0].shape[0] > 1:
                # more than one fire...
                # so set afters back to 0
                #import pdb; pdb.set_trace()
                timestack[idx][1:] =0
                # make a mask...
                mask = np.zeros_like(timestack)
                mask[idx[0][0]]=1
                # multiply mask by timestack
                timestack *= mask
                # make sure value is 1
                #timestack[idx[0][0]]=1
                # put the new timestack back
                #print timestack
                bools[:, x,y] = timestack
    #import pdb; pdb.set_trace()
    #print np.sum(bools, axis=0).max()
    return bools
