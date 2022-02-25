import pymangle
from  spherical_caps import cap
from astropy import units as u
import numpy as np
import os
import matplotlib.pyplot as plt

cap1 = np.asarray(cap(76*u.deg, 36*u.deg, 5))
cap2 = np.asarray(cap(75*u.deg, 35*u.deg, 5))
coord_list = np.asarray([[cap1, cap2]])

def create_mangle_file(array_of_polys, filename):
    """
    Function that takes an array consisting of 4-vector caps, and
    saves a .ply file in the correct format for pymangle to use.
    
    Inputs
    ------
    array_of_polys : ndarray
    - an array of first dimension n=number of polygons and second
    dimension n=number of caps for that polygon. For example, an
    array with 2 polygons with 3 caps each is; 
    np.asarray([[cap1, cap2, cap3],[cap4, cap5, cap6]])
    - the filename of the result, including extension .ply.
    
    Outputs
    -------
    A file saved to the location current_directory/filename

    """
    f = open(os.getcwd() + "/" + filename, "w")

    # KFH Set first line of ply file
    line1 = '{} polygons'.format(array_of_polys.shape[0])
    f.write(str(line1))
    
    # KFH Generate polygon entries for mangle object
    for i in np.arange(len(array_of_polys)):
        poly = i+1
        f.write('\n')
        f.write(str('polygon {} ({} caps, 1 weight, 0 pixel, 0 str):'
                    .format(poly, len(array_of_polys[:][0]))))

        # KFH For each polygon make a line in ply file for each cap
        for x in np.arange(len(array_of_polys[:][:])):
            print(x)
            cap = array_of_polys[:][x]
            f.write('\n')
            for i in cap:
                f.write(str(i) + ' ')
            
    f.close()

user = os.getenv('USER')
formatter = '/d/www/{}/public_html/week6'
webdir = formatter.format(user)

# KFH Generate caps
cap1 = np.asarray(cap(76*u.deg, 36*u.deg, 5))
cap2 = np.asarray(cap(75*u.deg, 35*u.deg, 5))
coord_list = np.asarray([[cap1, cap2]])

# KFH Generate mangle formatted files
create_mangle_file(coord_list, 'intersection.ply')
create_mangle_file(np.asarray([[cap1], [cap2]]), 'bothcaps.ply')

minter = pymangle.Mangle("intersection.ply")
mboth = pymangle.Mangle("bothcaps.ply")

ra_inter, dec_inter = minter.genrand(10000)
ra_both, dec_both = mboth.genrand(10000)


plt.clf()

fig = plt.figure(figsize=[10, 10])
ax = fig.add_subplot(111)
ax.scatter(ra_inter, dec_inter, s=0.7, alpha=0.5, color="green", label="Intersection")
ax.scatter(ra_both, dec_both, s=0.7, alpha=0.5, color="red", label="Both caps")
ax.set_xlabel("RA (deg)")
ax.set_ylabel("Dec (deg)")
ax.legend()

plt.savefig(webdir + '/intersection.png')

