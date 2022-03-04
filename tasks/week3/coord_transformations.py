import astropy.coordinates as ap
import math
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.time import Time
import numpy as np
import pandas as pd

def radec_to_cartesian(ra_in, dec_in):
    sc = ap.SkyCoord(ra=ra_in, dec=dec_in)
    print("Skycoord result:" + str(sc.cartesian))
    
    x_test = math.cos(sc.ra.radian)*math.cos(sc.dec.radian)
    y_test = math.sin(sc.ra.radian)*math.cos(sc.dec.radian)
    z_test = math.sin(sc.dec.radian)
    print("Adam's equations result:" + str(x_test) +", " +  str( y_test) + ", " + str( z_test))

def galactic_to_radec(l_in, b_in):
    sc = ap.SkyCoord(l_in*u.degree, b_in*u.degree, frame='galactic')
    print(sc.fk5)
    const = ap.get_constellation(sc.fk5)
    print("The constellation the GC is in is:" + str(const))

#radec_to_cartesian('00h45m30s', '+24d15m05s')    
galactic_to_radec(0, 0)
# The GC is near the edge of Sagittarius.

def galactic_coords_zenith():
    # First lets get our list of ra's every night
    dates = pd.date_range('2022-01-01', '2022-12-31')
    dates_str = np.asarray(dates.strftime('%Y-%m-%d'))
    midnights = []
    for i in dates_str:
        i = i+'T00:00:00'
        midnights.append(i)

    loc = ap.EarthLocation.of_address('Laramie, WY')
    full_data = []
    for i in midnights:
        coords = ap.SkyCoord(az=0*u.degree, alt=90*u.degree, frame='altaz', obstime=i, location=loc).galactic
        full_data.append([coords.l.degree, coords.b.degree])

    full_data = np.asarray(full_data)
    print(full_data)
    plt.clf()
    plt.scatter(full_data[:,0], full_data[:,1], s=1, color = "blue")
    plt.xlabel("Galactic coord l (degrees)")
    plt.ylabel("Galactic coord b (degrees)")
    plt.title("Galactic coords of zenith in Laramie")

    plt.savefig('/d/www/kianah/public_html/week3_galactic_coords_zenith2.png')

galactic_coords_zenith()
radec_to_cartesian('00h45m30s', '+24d15m05s')                                                               
galactic_to_radec(0, 0)
