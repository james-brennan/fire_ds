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


def burnIt_idea1(surface, size=100, DOB, num_seeds=10, temporal=True):
    """
    Generate burned areas that occur around the date of burn choosen
    easy way is to create a 3d boolean mask where only TRUE are where
    burn is (on DOB)


    This methods starts from a number of seeds and grows outwards randomly
    to produce burned areas...
    --- also has a temporal flag so that the burn spreads in time..
    """
    # 1. randomly allocate starting locations for seeds
    if temporal:
        #need to take account of time
        #burns = [xy, DOB, for xy in np.random.uniform(0, size, 2)]

        x,y,DOBs = np.random.random()
        burns.append(x,y,DOB)
    if !temporal:
        # just flatten the DOB to be the supplied DOB




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
    def __init__(self, outname, sizex, sizey, daily=True, fire_percent=20 ):
        pass




####################################
####################################
####################################



# lets just start of simple and generate a fire in the middle on a choosen dataset



def run_fire(num_days=20):
    """ Return of healthy vegetation length

        Computes the mixture of ash and recovered vegetation -->
        --> NOTE: obs very simple model of return of vegetation to health...
                -- no spectral change etc...
    """
    vegetation =  1.0 / (1.0 + np.exp(-np.linspace(-6,6,num_days)))
    char_ash = 1 - vegetation
    return vegetation, char_ash

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



# seems to only work in a loop atm!!!
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

    

### ok that worked now lets generate a 3d dataset
bands = 13
timesteps = 100
size=100
data = np.ones((timesteps, bands, size,size))


# first let's create some realistic scenes...

# can use perlin noise to generate a more realistic variation in vegetation reflectance
# -- eventually can a have a few different classes of vegetation
# ---   -- and use in mixture model...

# select bands...
vbands = veg_refl[::25]

# first fill the dataset with just normal leaf reflectance across all timesteps
surface_refl = (data[:].T*vbands[None, :].T).T



# created spatio-temporal noise patterns
noises = make_noises(timesteps,bands,size)

# scale down from (-1,1) to a smaller range ...
noises *= 0.1

# now add to the surface reflectance
surface_refl += noises


# 1. so now got some sort of realistic pattern of vegetation dynamics in the
# background...




import sys
sys.setrecursionlimit(1500)

class aFire():
    def __init__(self, DOB, x,y):
        self.x = [x]
        self.y = [y]
        self.DOB = [DOB]
        # run the fire
        self._spread(self.DOB[0], self.x[0], self.y[0], decay=1)

    def _spread(self, DOB, x,y, decay):
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






size=100
seeds = 10
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
bools = np.zeros((timesteps, size, size)).astype(np.bool)
for fire in fires:
    # make sure placing in right place!
    for burnday in zip(fire.DOB, fire.x, fire.y):
        day = burnday[0]
        x = burnday[1]
        y = burnday[2]
        try:
            bools[day,x,y] = True
        except:
            pass


# so bools provides the DOB for each pixel now...

# So final step is to use the burn spectral model
# to change the surface_refl on DOB
# current thought is just a simple loop -- not sure on clever numpy vectorisation
# yet
ash_spectrum = soil_ash_char2[::25]
fires_locations = np.where(bools==True)
for i in xrange(len(fires_locations[0])):
    #
    dob = int(fires_locations[0][i])
    x =  fires_locations[1][i]
    y =  fires_locations[2][i]
    spectral_response_weights = run_fire(num_days=50)
    # do mixture model
    if not (dob+50 > timesteps):
        surface_refl[dob:dob+50, :, x, y] = (
                        (spectral_response_weights[0]*surface_refl[dob:dob+50, :, x, y].T)+
                        spectral_response_weights[1] * ash_spectrum[:,np.newaxis]).T
