import matplotlib.pyplot as plt
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import search_around_sky
import numpy as np

coo1 = SkyCoord(ra = 263.75, dec = -17.9, frame='icrs', unit='deg')
coo2 = SkyCoord('20h24m59.9s', '+10d06m00s', frame = 'icrs', unit=(u.hourangle, u.deg))

coord1 = coo1.cartesian
coord2 = coo2.cartesian

c1 = np.asarray([coord1.x.value, coord1.y.value, coord1.z.value])
c2 = np.asarray([coord2.x.value, coord2.y.value, coord2.z.value])

dot_prod = np.dot(c1, c2)
z_ang = np.arccos(dot_prod/(np.linalg.norm(c1)*np.linalg.norm(c2)))*360/(2*np.pi)

#KFH Now we can check the results against the separation function
sep = coo1.separation(coo2).degree
print("Angular separation, manual and skycoords calculated;")
print(z_ang, sep)

# 2.
# KFH First lets generate the random ra and dec values
rand_ra1 = np.random.random(100)+2
rand_ra2 = np.random.random(100)+2

rand_dec1 = 4*np.random.random(100)-2
rand_dec2 = 4*np.random.random(100)-2

objs1 = SkyCoord(ra=rand_ra1*u.hour, dec=rand_dec1*u.degree, frame='icrs')
objs2 = SkyCoord(ra=rand_ra2*u.hour, dec=rand_dec2*u.degree, frame='icrs')


idx1, idx2, d2d, d3d = objs2.search_around_sky(objs1, seplimit = 0.167*u.degree)

# KFH Plot the datasets and separation datasets
plt.clf()
plt.scatter(objs1[idx1].ra.hour, objs1[idx1].dec.degree, color='m',
            marker='o', label='points within 10 arcmin of other dataset')
plt.scatter(objs2[idx2].ra.hour, objs2[idx2].dec.degree, color='m',
            marker='o')
plt.scatter(rand_ra1, rand_dec1, color='b', marker='+', label='rand points1',
            alpha=0.5)
plt.scatter(rand_ra2, rand_dec2, color='r', marker='v', label = 'rand points2',
            alpha=0.5)
plt.xlabel('RA (hours)')
plt.ylabel('Dec (degrees)')
plt.legend(loc='upper right')

plt.savefig('/d/www/kianah/public_html/week4_sphere_dist.png')
print('Figure found at /d/www/kianah/public_html/week4_sphere_dist.png')

