import numpy as np
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
from astropy import units as u

def ra_vect_array(ra, neg=False):
    """
    Function to convert an RA bound into a spherical cap vector 4-array
    in cartesian coordinates
    Inputs:
    ra: the RA bound of your cap in xxhxxmxxs format
    """
    # KFH Find the center of the bounded area
    new_ra = str(int(ra[0:2])+6) + ra[2:]
    h = np.where(neg, -1, 1)
    # KFH Calculate the cartesian coordinates of spherical cap
    area = SkyCoord(ra=new_ra, dec = '+00d00m00s', frame='icrs').cartesian
    vec4array=[area.x.value, area.y.value, area.z.value, h]
    return(vec4array)

def dec_vect_array(dec, neg=False):
    """                                                                          
    Function to convert a dec bound into a spherical cap vector 4-array          
    in cartesian coordinates                                                     
    Inputs:                                                                      
    dec: the declination bound of your cap in xxdxxmxxs format                             
    """

    # KFH Calculate the cartesian coordinates of spherical cap                   
    area = SkyCoord(ra=0*u.degree, dec = 90*u.degree, frame='icrs').cartesian
    h = np.where(neg, -(1 - np.sin(int(dec[1:3])*np.pi/180)), 1 - np.sin(int(dec[1:3])*np.pi/180))
    vec4array=[area.x.value, area.y.value, area.z.value, h]
    return(vec4array)

def cap(ra, dec, radius, neg=False):
    """
    Function to convert an RA and dec into a spherical cap vector 4-array
    in cartesian coordinates, centered around the input RA,dec.
    """
    area = SkyCoord(ra, dec, frame='icrs').cartesian
    h = np.where(neg, -(1-np.cos(radius*np.pi/180)),(1-np.cos(radius*np.pi/180)))
    vec4array=[area.x.value, area.y.value, area.z.value, h]
    return(vec4array)

if __name__ == '__main__':
    ra_vect_array('05h00m00s')
    dec_vect_array('+36d00m00s')
    cap('05h00m00s', '+36d00m00s', 1)
