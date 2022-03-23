import numpy as np
import os
import pandas as pd
from astropy.coordinates import Angle


ra = "16h35m26s"
dec = "+09d47m53s"

color = [14.397, -0.192, -0.974, -0.093, -0.116]
v = color[0]
b = v + color[1]
u = color[2] + b
r = -(color[3] - v)
i = -(color[4] - r)

def one_pt_sdss_data(ra, dec, textfile):
    angles = Angle([ra, dec]).degree
    ra_deg = angles[0]
    dec_deg = angles[1]
    print(angles)
    if (os.path.exists(textfile)):
        os.remove(textfile)

    # KFH Query sdss database and append results of each object query to textfile
    [os.system('python sdssDR9query.py {} {} >> {}'.format(ra_deg, dec_deg, textfile))]

def ubvri_to_ugriz(list_of_filters, ra, dec, textfile):
    """
    list_of_filters should be of the format [V, B-V, U-B, V-R, R-I]
    """

    if list_of_filters[4] < 1.15 and list_of_filters[2] < 0:
        # KFH Calculate ugriz filters from uvbri
        ug = 1.28*list_of_filters[2] + 1.14
        gr = 1.09*list_of_filters[1] - 0.23
        ri = 0.98*list_of_filters[4] - 0.22
        rz = 1.69*list_of_filters[4] - 0.42
        g = list_of_filters[0] + 0.64*list_of_filters[1] - 0.13
        rr = list_of_filters[0] - 0.46*list_of_filters[1] + 0.11
        z = -(rz - rr)
        
    else:
        print("Out of color range")

    # KFH query sdss with point
    one_pt_sdss_data(ra, dec, textfile)

    # KFH read in resulting file
    data = pd.read_csv(textfile, sep=',', header=None, index_col=False)
    data.columns = ['ra', 'dec', 'u', 'g', 'r', 'i', 'z', 'sep']
    
    print("Calculated g and z:" + str(g) + ", " + str(z))
    print("SDSS  g and z:" + data['g'].to_string(index=False) + ", " + data['z'].to_string(index=False))
    print(data['ra'].to_string(index=False), data['dec'].to_string(index=False))
    
ubvri_to_ugriz(color, ra, dec, "pg1633_099a.txt")
