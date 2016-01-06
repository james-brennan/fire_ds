#-*- coding: utf-8 -*-
"""

    Fire class implementation
    =========================

    generates burn scars which evolve in time+space

"""
from scipy.signal import convolve
import scipy.interpolate
import sys
import numpy as np
sys.setrecursionlimit(10000)
import scipy.ndimage.filters

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

        NEED TO THINK A BIT MORE ABOUT THIS...
        ESSENTIALLY WANT TO FILL IN A SHAPE better
        -- SO NEED TO INTERPOLATE BETWEEN DIRECTIONS?
        -- FIRE NEEDS TO SPREAD IN ALL DIRECTIONS MORE...
            -- DOWNWEIGHT IMPORTANCE OF BIAS?
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
            wind = np.array([-0.75, 0])
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
            #
            scaled_dist2 = scipy.ndimage.filters.uniform_filter1d(scaled_dist, 8, mode='constant')
            #import pdb; pdb.set_trace()
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
