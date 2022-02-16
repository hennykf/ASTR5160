import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u

def ra_vect_array(ra):
    """
    Function to convert an RA bound into a spherical cap vector 4-array
    in cartesian coordinates
    Inputs:
    ra: the RA bound of your cap in xxhxxmxxs format
    """
    # KFH Find the center of the bounded area
    new_ra = str(int(ra[0:2])+6) + ra[2:]

    # KFH Calculate the cartesian coordinates of spherical cap
    area = SkyCoord(ra=new_ra, dec = '+00d00m00s', frame='icrs').cartesian
    vec4array=[area.x.value, area.y.value, area.z.value, 1]
    print("RA Vector 4-array: " + str(vec4array))
    return(vec4array)

def dec_vect_array(dec):
    """                                                                          
    Function to convert a dec bound into a spherical cap vector 4-array          
    in cartesian coordinates                                                     
    Inputs:                                                                      
    dec: the declination bound of your cap in xxdxxmxxs format                             
    """

    # KFH Calculate the cartesian coordinates of spherical cap                   
    area = SkyCoord(ra=0*u.degree, dec = 90*u.degree, frame='icrs').cartesian
    h = 1 - np.sin(int(dec[1:3])*np.pi/180)
    vec4array=[area.x.value, area.y.value, area.z.value, h]
    print("Dec Vector 4-array: " + str(vec4array))
    return(vec4array)

ra_vect_array('05h00m00s')
dec_vect_array('+36d00m00s')
