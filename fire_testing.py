"""

  New fire spread model...


  Celluar automata approach

  implementation of Drossel_and_Schwabl which works!

"""
import scipy.ndimage.filters
from numba import jit

class SeedError(Exception):
    pass


def fire(size=100, steps=100, max_size=900):

    """
     1. setup

     1. Make initial grid...
    """
    DOBS = []
    xs = []
    ys =[]
    tree, burning, space = 1,2,0
    neighbourhood = np.array(((-1,-1), (-1,0), (-1,1),
                    (0,-1),          (0, 1),
                    (1,-1),  (1,0),  (1,1)))
    #neighbourhood = np.array(( (-1,0),
    #                    (0,-1),          (0, 1),
    #                            (1,0),  ))
    # do setup stuff
    #grid = np.zeros((xSize, ySize,)).astype(np.int)
    l0 = np.random.randint(0,2,size=(size,size)).astype(np.int)
    # do some smoothing on it to make it more realistic...
    #import pdb; pdb.set_trace()
    l0 = scipy.ndimage.filters.median_filter(l0, 5).astype(np.int)
    l0 = scipy.ndimage.filters.median_filter(l0, 15).astype(np.int)
    init = l0.copy()
    #import pdb; pdb.set_trace()
    #l0 = np.ones((size, size,)).astype(np.int)
    """
    2. Seed fire in a vegetation location
    """

    idx = np.where(l0==1)[0]
    idy = np.where(l0==1)[1]
    #import pdb; pdb.set_trace()
    if idx.shape[0] == 0:
        # no trees for some reason
        return 'what'
    x0 = idx[np.random.randint(1,idx.shape[0])]
    y0 = idy[np.random.randint(1,idy.shape[0])]
    l0[x0,y0] = 2

    """
    Some constants to add randomness
    """
    f = 0.99

    sto = []
    """

    3. Run simulation for timesteps
    """
    fails = 0
    burnt_count = 0
    doy = 0
    lt = l0.copy()
    for t in xrange(steps):
      if burnt_count < max_size:
          for x in range(1, size-1):
              for y in range(1, size-1):
                  # check if cell is burning
                  if l0[(x,y)] == burning:
                      lt[(x,y)] = space
                  elif l0[(x,y)] == space:
                      lt[(x,y)] = space # keep as space...
                  elif l0[(x,y)] == tree:
                      # check if a neighbour is burning
                       #import pdb; pdb.set_trace()
                       if np.any([l0[(x+dx,y+dy)]==burning for
                                     dx,dy in neighbourhood]):
                            # a neighbour is burning..
                            if np.random.random() < f:
                                lt[(x,y)] = burning
                                # store out to arrays...
                                burnt_count += 1
                                # vary up and down day a bit...

                                DOBS.append(doy + np.random.randint(-1,2))
                                xs.append(x)
                                ys.append(y)
                            else:
                                # what about burning after a day?
                                #delayed_burn.append((x,y))
                                lt[(x,y)] = space
          # sort out cells that are meant to have been
          # burnt in the last time step...
          # idea set to burning again?
          #if len(delayed_burn) > 0:
              #xx,yy = delayed_burn.pop()
              #lt[(xx,yy)] = burning
          # after each iteration re-assign grid...
          print t
          doy +=1
          #import pdb; pdb.set_trace()
          l0 = lt.copy()
          #sto.append(lt)
          # check seeding worked properly
          # burned area should increase each time right...
          if abs(l0.sum() - lt.sum()) == 0:
              fails += 1
          if fails > 20:
              # return what you'ev got
              print 'here'
              #return zip(DOBS, xs, ys), init
              #raise SeedError("Seeding has failed...")
    return zip(DOBS, xs, ys), init


size=200
f, init = fire(size=size)
f = np.array(f)

burn_map = np.zeros((size, size))
burn_map[f.T[1].astype(int), f.T[2].astype(int)] = f.T[0]

plt.figure()
plt.subplot(121)
plt.imshow(burn_map)
plt.subplot(122)
plt.imshow(init)
