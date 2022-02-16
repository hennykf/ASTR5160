from numpy.random import random
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

ra = 360.*(random(1000000))
dec = (180/np.pi)*np.arcsin(1.-random(1000000)*2.)

#KFH Convert random coords to healpix
pix = hp.ang2pix(1, ra, dec, lonlat=True)
print("Nside=1 pixel area: "+ str(hp.nside2pixarea(1)))

#KFH Find the number of points in each pixel
hist, bin_edges = np.histogram(pix, 11)
print("Number of points in each pixel, pixels 1-12: " + str(hist))

#KFH There are approximately equal number of points in every bin,
# showing that the pixels are equal-area.
