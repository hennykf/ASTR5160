import pymangle
from spherical_caps import cap, ra_vect_array, dec_vect_array
from mangle import create_mangle_file
from astropy import units as u
import numpy as np
import os
import  matplotlib.pyplot as plt

# KFH Generate coordinate boundaries and create mangle file
racap1 = ra_vect_array('05h00m00s')
racap2 = ra_vect_array('06h00m00s')

deccap1 = dec_vect_array('+30d00m00s')
deccap2 = dec_vect_array('+40d00m00s')

area_str = (75*np.pi/180 - 90*np.pi/180)*(np.sin(30*np.pi/180) - np.sin(40*np.pi/180))
coord_array = np.asarray([[racap1, racap2, deccap1, deccap2]])
create_mangle_file(coord_array, 'rectangle.ply', weight=.9, area=[area_str])

racap3 = ra_vect_array('10h00m00s')
racap4 = ra_vect_array('12h00m00s')

deccap3 = dec_vect_array('+60d00m00s')
deccap4 = dec_vect_array('+70d00m00s')

area_str2 = (150*np.pi/180 - 180*np.pi/180)*(np.sin(60*np.pi/180) - np.sin(70*np.pi/180))
                         
coord_array2 = np.asarray([[recap1, recap2, deccap1, deccap2],
                           [racap3, racap4, deccap3, deccap4]])
create_mangle_file(coord_array2, 'two_rectangles.ply', weight=[0.9, 0.2], area=[area_str, area_str2])
