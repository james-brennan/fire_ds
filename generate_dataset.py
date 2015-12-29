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
#import noise
from scipy.signal import convolve
import scipy.interpolate
import sys
import numpy as np
sys.setrecursionlimit(1500)






class params(object):
    def __init__(self):
        """ params past to generate a dataset
        """
        self.size = 800,800
        self.repeat = 1

class experiment(params):
    """
    Vary the parameters to produce an experiment
    """
    pass


class dataset(object):
    def __init__(self, outname='nothing.tif',timesteps=100, bands=13, sizex=100, sizey=100, spatial_res=1 ):
        self.xSize = sizex
        self.ySize = sizey
        self.timesteps = timesteps
        self.bands = bands
        self.spatial_res = spatial_res
        self.surface_rho = None

    def make_stack(self):
        self.surface_rho = np.ones((self.timesteps, self.bands, self.xSize, self.ySize))

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
        pass

    def model_fires(self):
        """

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

        # ADD NOISE
        # created spatio-temporal noise patterns
        #noises = make_noises(timesteps,bands,size)

        # scale down from (-1,1) to a smaller range ...
        #noises *= 0.1

        # now add to the surface reflectance
        #self.surface_rho += noises

        # now generate some fires
        self.fires = burnIt_idea1(self, seeds=3)
        self.fire_locs = np.where(self.fires)
        #import pdb; pdb.set_trace()
        # Now run mixture model for spectral response to fire
        dob = 10
        self.surface_rho = spectral_fire_model(self.surface_rho, dob, self.timesteps,
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

def make_noises(timesteps, bands, size):
    noises = np.zeros((timesteps, bands, size,size))
    for t,z in enumerate(np.linspace(0,1,timesteps)):
        #print t
        for i,x in enumerate(np.linspace(0,1,size)):
            for j,y in enumerate(np.linspace(0,1,size)):
                #print x,y
                # increasign the number of octaves increases the noise...
                noises[t,:, j,i] = noise.snoise3(x,y,z, octaves=5) # not sure what octave to pick...
    return noises


class aFire():
    def __init__(self, DOB, x,y):
        self.x = [x]
        self.y = [y]
        self.DOB = [DOB]
        # run the fire
        self._spread(self.DOB[0], self.x[0], self.y[0], decay=1)

    def _spread(self, DOB, x,y, decay=20):
        """

        The key function...
        this spreads the fire with some randomness
        while also decaying so it doesn't grow forever...
        """
        if decay > 0:
            # 1. choose a direction to expand
            dx, dy = np.random.choice([-1,0,1],2)
            # verify it's not been burnt during this fire
            # so cant have both being true!
            #import pdb; pdb.set_trace()
            # think unfortunately need to che k
            if not (x+dx, y+dy) in zip(self.x, self.y):
                # decay a bit...
                decay -= 0.001
                # add this location to the store in the class..
                self.x.append(x+dx)
                self.y.append(y+dy)
                self.DOB.append(DOB+0.02) # spread more one day using a decimal increment for DOB
                # recurse from this location...
                #print decay
                self._spread(DOB+0.02, x+dx, y+dy, decay)
            else:
                # been burnt before
                # not sure what to do?
                # pick again
                self._spread(DOB+0.02, x, y, decay-0.001)
        else:
            return None


def burnIt_idea1(dataset, seeds=10, temporal=True):
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
        this_fire = aFire(burns[f][0], burns[f][1][0], burns[f][1][1])
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
    return bools


def spectral_temporal_response(num_days=20):
    """ Return of healthy vegetation length

        Computes the mixture of ash and recovered vegetation -->
        --> NOTE: obs very simple model of return of vegetation to health...
                -- no spectral change etc...
    """
    vegetation =  1.0 / (1.0 + np.exp(-np.linspace(-6,6,num_days)))
    char_ash = 1 - vegetation
    return vegetation, char_ash


def spectral_fire_model(surface_refl, dob, timesteps, fires_locations, ash_spectrum):
    for i in xrange(len(fires_locations[0])):
        #
        dob = int(fires_locations[0][i])
        x =  fires_locations[1][i]
        y =  fires_locations[2][i]
        spectral_response_weights = spectral_temporal_response(num_days=50)
        # do mixture model
        if not (dob+50 > timesteps):
            surface_refl[dob:dob+50, :, x, y] = (
                            (spectral_response_weights[0]*surface_refl[dob:dob+50, :, x, y].T)+
                            spectral_response_weights[1] * ash_spectrum[:,np.newaxis]).T
    return surface_refl


def main():
    ds = dataset()
    ds.model_fires()
    import pdb; pdb.set_trace()
    ds.BA_to_tif()
    ds.surface_refl_to_tif()

if __name__ == "__main__":
    main()    
