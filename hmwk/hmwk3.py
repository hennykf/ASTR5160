import numpy as np
from tasks.week6.mangle import create_mangle_file
from tasks.week6.spherical_caps import  ra_vect_array, dec_vect_array, cap
from astropy.coordinates import Angle, SkyCoord
import pymangle
import astropy.units as u
from hmwk.hmwk2 import rect_area
from astropy.table import Table
import matplotlib.pyplot as plt

"""
This script constructs mangle masks for a rectangular area and several plates mostly
within the rectangle, and finds the area where the plates and rectangle overlap.
It then plots quasars, indicating those quasars that fall within the area of the plates.
"""
    
# KFH Use functions from week6 to generate caps
# KFH Generate caps for rectangular survey area
ra_min_cap = ra_vect_array('10h15m00s')
ra_max_cap = ra_vect_array('11h15m00s', neg=True)
dec_min_cap = dec_vect_array('+30d00m00s')
dec_max_cap = dec_vect_array('+40d00m00s', neg=True)

# KFH Generate caps for plates
plate_cap1 = cap(155*u.deg, 34*u.deg, 2)
plate_cap2 = cap(159*u.deg, 36*u.deg, 2)
plate_cap3 = cap(163*u.deg, 34*u.deg, 2)
plate_cap4 = cap(167*u.deg, 36*u.deg, 2)

# KFH Create mangle of rectangle only
create_mangle_file(np.array([[ra_min_cap, ra_max_cap, dec_min_cap, dec_max_cap]]), 'survey_area.ply')

# KFH Create array of polys that each include the survey rectangle and one plate
poly_array = np.asarray([[plate_cap1, ra_min_cap, ra_max_cap, dec_min_cap, dec_max_cap],
                    [plate_cap2, ra_min_cap, ra_max_cap, dec_min_cap, dec_max_cap],
                    [plate_cap3, ra_min_cap, ra_max_cap, dec_min_cap, dec_max_cap],
                    [plate_cap4, ra_min_cap, ra_max_cap, dec_min_cap, dec_max_cap]])
create_mangle_file(poly_array, "plate_survey_overlap.ply", weight=[1,1,1,1], area=[0,0,0,0])

# KFH Create mangle objects
rect = pymangle.Mangle('survey_area.ply')
m = pymangle.Mangle('plate_survey_overlap.ply')

# KFH Randomly generate equal-distributed points in rectangle
points = 100000
ra, dec = rect.genrand(points)

# KFH Create mask of points in plate areas
in_mask = m.contains(ra, dec)
fract = np.count_nonzero(in_mask==1) / points

# KFH Find area of rectangular survey
full_survey_area = rect_area(Angle('10h15m00s').degree, Angle('11h15m00s').degree,
                             Angle('+30d00m00s').degree, Angle('+40d00m00s').degree)

# KFH Calculate number of square degrees in plates-survey overlap
num_degrees_in_mask = full_survey_area  * fract
print(str(num_degrees_in_mask) + " square degrees in plate+survey overlapping area")

# KFH Import quasar file
qsos = Table.read("/d/scratch/ASTR5160/week8/HW3quasarfile.dat", format='ascii.no_header',
                  names=['position'])

# KFH Convert to more conventional skycoord format
ra = ['{}h{}m{}s'.format(q[0:2], q[2:4], q[4:9]) for q in qsos['position']]
dec = ['{}d{}m{}s'.format(q[9:12], q[12:14], q[14:19]) for q in qsos['position']]
coords = SkyCoord(ra, dec, frame='icrs')
qsos.add_columns([coords.ra, coords.dec],  names=['ra', 'dec'])

# KFH Create subsets of quasars in and not in plate areas
in_plate = m.contains(qsos['ra'], qsos['dec'])
qsos_in_plate = qsos[in_plate==1]
qsos_nin_plate = qsos[in_plate==0]

# KFH Calculate number density per square degree
num_density = len(qsos_in_plate) / num_degrees_in_mask

# KFH Plot figure
plt.clf()
fig = plt.figure(figsize=[10,10])
ax = fig.add_subplot(111)

ax.scatter(qsos_in_plate['ra'], qsos_in_plate['dec'], label = "Points in plate area", s=1)
ax.scatter(qsos_nin_plate['ra'], qsos_nin_plate['dec'], label = "Points outside of plates", s=1)
ax.set_xlabel("RA (deg)")
ax.set_ylabel("Dec (deg)")
ax.legend()
ax.set_title("Quasars, mask area={}, mask quasar number density={} qsos/degree".format(str("%.2f" % num_degrees_in_mask),
                                                                                    str("%.2f" % num_density)))
fig.savefig("/d/www/kianah/public_html/hmwk/quasars_in_survey.png")
print("Figure saved to: /d/www/kianah/public_html/hmwk/quasars_in_survey.png")






