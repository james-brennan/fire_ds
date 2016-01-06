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
import noise
from scipy.signal import convolve
import scipy.interpolate
import sys
import numpy as np
sys.setrecursionlimit(10000)
import matplotlib.pyplot as plt
import scipy.ndimage.filters


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
        self.modelled_fire = False
        self.fires = None

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
        """

        idea for now is to re-use perlin noise generator
        With a threshold for how many clouds there should be:
            -- could be semi-seasonal + random element

        """
        # model of cloud cover across timesteps
        #       maximum seasonal component 80%          random weather
        cloud =0.4* np.sin(np.linspace(0,np.pi,self.timesteps)) + np.random.normal(0.2, 0.1, self.timesteps)
        # make sure max is one
        cloud[cloud>1]=1
        # generate some new perlin noise
        # for whole dataset!
        # try and apply day by day to save memory though

        the_noise = make_noises(self.timesteps, self.xSize, time_multiple=2)
        #import pdb; pdb.set_trace()
        # idea is to threshold the perlin noise by the percentile of cloud_cover
        # this produces a binary mask that should produce realistic looking
        # clouds
        percentiles = np.array([np.percentile(no, cl*100) for no, cl in zip(the_noise, cloud)])

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
        self.surface_rho = (self.surface_rho[:].T*self.cloud_mask[:,None].T).T
        # make it into a masked array for clarity
        #import pdb; pdb.set_trace()
        self.surface_rho = np.ma.array(self.surface_rho, mask=self.surface_rho==0)
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
            self.fires = burnIt_idea1(self, seeds=5, numPixels=to_burn)
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
                noises[t, j,i] = noise.snoise3(x,y,time_multiple*z, octaves=5) # not sure what octave to pick...
    return noises


class aFire(object):
    def __init__(self, DOB, x,y, numPixels=1000):
        self.x = [x]
        self.y = [y]
        self.DOB = [DOB]
        # run the fire
        self._spread3(self.DOB[0], self.x[0], self.y[0], numPixels=numPixels)

    def _spread3(self, DOB, x,y, numPixels=500, count=0):
        """

        The key function...
        this spreads the fire with some randomness
        to fill a certain number of pixels as supplied...
        """
        #import pdb; pdb.set_trace()
        if numPixels - count >= 0:
            # 1. choose a direction to expand
            """
            Key part of algorithm:
                Here we are choosing a surrounding box cell to move to
                    (or stay put).
                the various permuations of +x and +y
                are:
                 x    y
                 ======
                 0    0  -- current location
                 1    1
                -1   -1
                 1   -1
                -1    1
                -1    0
                 1    0
                 0    1
                 0    -1
            """
            movements = np.array([[0,0],[1,1], [-1,-1],
                                  [1,-1],[-1,1],[-1,0],
                                  [1,0],[0,1], [0,-1]] )

            choice = np.random.choice([0,1,2,3,4,5,6,7,8])
            if choice > 0:
                #import pdb; pdb.set_trace()
                dx, dy = movements[choice]
                count += 1
                # add this location to the store in the class..
                self.x.append(x+dx)
                self.y.append(y+dy)
                #print choice, x+dx, y+dy, decay
                #print decay
                self.DOB.append(DOB+0.02) # spread more one day using a decimal increment for DOB
                # recurse from this location...
                #print decay
                self._spread3(DOB+0.02, x+dx, y+dy, numPixels, count)
            else:
                # choice is 0...
                # don't want to recurse out of function
                # but want to try again?
                self._spread3(DOB+0.02, x, y, numPixels, count)
        else:
            return None

    def _spread4(self, DOB, x,y, numPixels=500, count=0):
        """

        The key function...
        this spreads the fire with some randomness
        to fill a certain number of pixels as supplied...

        THIS VERSION ALSO FOLLOWS A predominant wind-direction
        which varies across days!

        from

        http://jflevesque.com/2012/12/06/far-cry-how-the-fire-burns-and-spreads/

        """
        #import pdb; pdb.set_trace()
        if numPixels - count >= 0:
            # 1. choose a direction to expand
            """
            Key part of algorithm:
                Here we are choosing a surrounding box cell to move to
                    (or stay put).
                the various permuations of +x and +y
                are:
                 x    y
                 ======
                 0    0  -- current location
                 1    1
                -1   -1
                 1   -1
                -1    1
                -1    0
                 1    0
                 0    1
                 0    -1
            """
            movements = np.array([[1,1], [-1,-1],
                                  [1,-1],[-1,1],[-1,0],
                                  [1,0],[0,1], [0,-1]] )
            """

            To weight the directions more we derive the dot product between
            the wind vector and the direction vectors...

            we then sample a direction choice from this distribution

            """
            wind = np.array([-0.75+0.15*np.sin(-DOB),0.75+0.15*np.cos(-DOB)])
            print wind
            #
            wind_distance = np.array([np.dot(wind,mov) for mov in movements])
            idx = xrange(8)
            # need to do some annoying sorting atm
            to_sort = zip(idx, wind_distance)
            sorted_winds = sorted(to_sort, key=lambda x: x[1], reverse=True)
            sorted_dis = np.array([t[1] for t in sorted_winds])
            chs = np.array([t[0] for t in sorted_winds])
            # normalise between 0 and 1
            scaled_dist = 1 - (((1 - 0) * (np.max(sorted_dis) - sorted_dis)) / (np.max(sorted_dis) - np.min(sorted_dis)))
            # smooth it down a bit...
            #import pdb; pdb.set_trace()
            scaled_dist2 = scipy.ndimage.filters.uniform_filter1d(scaled_dist, 5)
            # weight as a distribution
            wind_dir_distr = scaled_dist2 / np.sum(scaled_dist2)
            #import pdb; pdb.set_trace()
            choice = np.random.choice([0,1,2,3,4,5,6,7],  p=wind_dir_distr)
            # now convert choice back to right index
            choice = chs[choice]
            #import pdb; pdb.set_trace()
            dx, dy = movements[choice]
            count += 1
            # add this location to the store in the class..
            self.x.append(x+dx)
            self.y.append(y+dy)
            #print choice, x+dx, y+dy, decay
            #print decay
            self.DOB.append(DOB+0.02) # spread more one day using a decimal increment for DOB
            # recurse from this location...
            #print decay
            self._spread4(DOB+0.02, x+dx, y+dy, numPixels, count)
        else:
            return None

    def spread_Drossel_and_Schwabl(self, DOB, x0,y0, xSize, ySize, numPixels=500,):
            """

            use the forest fire model of Drossel_and_Schwabl (1992) to
            simulate. Is a celluar automata model

            re-implemented from http://rosettacode.org/wiki/Forest_fire#Python

            """
            tree, burning, space = 1,2,0
            neighbourhood = np.array(((-1,-1), (-1,0), (-1,1),
                            (0,-1),          (0, 1),
                            (1,-1),  (1,0),  (1,1)))
            # do setup stuff
            #grid = np.zeros((xSize, ySize,)).astype(np.int)
            grid = np.random.randint(0,2,size=(xSize,ySize)).astype(np.int)
            # do some clumping
            import scipy.ndimage.filters
            #import pdb; pdb.set_trace()
            grid = scipy.ndimage.filters.median_filter(grid, 15).astype(np.int)
            grid = scipy.ndimage.filters.median_filter(grid, 15).astype(np.int)
            # initally all cells are trees
            #grid[:, :] = 1
            # choose a location for fire to start
            # must be where trees are
            idx = np.where(grid==1)[0]
            idy = np.where(grid==1)[1]
            import pdb; pdb.set_trace()
            if idx.shape[0] == 0:
                # no trees for some reason
                return None
            x0 = idx[np.random.randint(idx.shape[0])]
            y0 = idy[np.random.randint(idy.shape[0])]
            grid[x0,y0] = 2

            newgrid = grid.copy()
            count = 0

            # some constants
            # not sure what these do atm
            p = 0.5
            f = 0.001
            date = 0
            its = 3
            print numPixels
            while numPixels - count > 0:
                # still got pixels to burn...
                # loop over all pixels oh god...
                  for x in range(xSize-1):
                      for y in range(ySize-1):
                          if grid[(x,y)] == burning:
                              newgrid[(x,y)] = space
                          elif grid[(x,y)] == space:
                              newgrid[(x,y)] = tree if np.random.random()<= p else space
                          elif grid[(x,y)] == tree:
                              #import pdb; pdb.set_trace()
                              if np.any([grid[(x+dx,y+dy)]==burning for
                                            dx,dy in neighbourhood]):
                                           # a neighbour is burning..
                                           newgrid[(x,y)] = burning
                                           # ALSO add out to self.fires for this date...
                                           self.DOB.append(DOB+date) # assumed spread rate
                                           self.x.append(x)
                                           self.y.append(y)
                                           count += 15
                                           #print True

                              else:
                                        grid[(x,y)]=tree
                                        #print False
                  # after each iteration replace grid with the new grid
                  #import pdb; pdb.set_trace()
                  grid = newgrid

                  date += 0.25
                  count +=100
                  print date, self.x[-1],self.y[-1],  numPixels - count
            print '+++++======='
            #import pdb; pdb.set_trace()
            return None

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
    ds.model_fires(30)
    import pdb; pdb.set_trace()
    #ds.model_cloud_cover()
    ds._save_to_gif()
    #ds.BA_to_tif()
    #ds.surface_refl_to_tif()

if __name__ == "__main__":
    main()
