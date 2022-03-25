import numpy as np
from tasks.week6.mangle import create_mangle_file
from tasks.week6.spherical_caps import  ra_vect_array, dec_vect_array, cap
from astropy.coordinates import Angle
import pymangle
import astropy.units as u
from hmwk.hmwk2 import rect_area, populate_rect, compare_area
from numpy.random import random

def populate_sphere(num_points):
    ra = 360.*(random(num_points))-180
    dec = (180/np.pi)*np.arcsin(1.-random(num_points)*2.)
    return(ra, dec)
    
# KFH Use functions from week6 to generate caps for whole area  and create a mangle file
ra_min_cap = ra_vect_array('10h15m00s')
ra_max_cap = ra_vect_array('11h15m00s')
dec_min_cap = dec_vect_array('+30d00m00s')
dec_max_cap = dec_vect_array('+40d00m00s')
coords = ['10h15m00s', '11h15m00s', '30d00m00s', '40d00m00s']
coords_deg = [Angle(d).value for d in coords]
ar = rect_area(coords_deg[0], coords_deg[1], coords_deg[2], coords_deg[3])

create_mangle_file(np.asarray([[ra_min_cap, ra_max_cap, dec_min_cap, dec_max_cap]]),
                   'survey_area.ply', area = [ar])

coords_plates = [[155,34], [159,36],[163, 34], [167, 36]]
plate_cap1 = cap(155*u.deg, 34*u.deg, 2)
plate_cap2 = cap(159*u.deg, 36*u.deg, 2)
plate_cap3 = cap(163*u.deg, 34*u.deg, 2)
plate_cap4 = cap(167*u.deg, 36*u.deg, 2)

create_mangle_file(np.asarray([[plate_cap1, ra_min_cap, ra_max_cap, dec_min_cap, dec_max_cap],
                               [plate_cap2, ra_min_cap, ra_max_cap, dec_min_cap, dec_max_cap],
                               [plate_cap3, ra_min_cap, ra_max_cap, dec_min_cap, dec_max_cap],
                               [plate_cap4, ra_min_cap, ra_max_cap, dec_min_cap, dec_max_cap]]),
                   "plate_survey_overlap.ply")


# KFH Generate random points on a sphere and create mangle mask
points = 10000000
ra, dec = populate_sphere(points)
m = pymangle.Mangle('survey_area.ply')

in_mask = m.contains(ra, dec)
fract = np.count_nonzero(in_mask==1) / points

num_degrees_in_sphere = 4*np.pi*(180/np.pi)**2
print(num_degrees_in_sphere)
print(fract)




