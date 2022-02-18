import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from  astropy.time import Time
from calendar import monthrange
import argparse

def best_airmass(month):
    """
    Accepts a month, and finds the lowest-airmass quasar out of a set file of quasars
    at Kitt Peak National Observatory at 11 PM every night that month, in 2022.

    Input:
    month : int
    - An integer from 1-12 corresponding to a month of the year

    Output:
    An astropy Table containing the date, quasar position, and airmasses of the lowest-airmass
    quasar each night of that month.
    """
    month = int(month)
    quasars = Table.read('/d/scratch/ASTR5160/week4/HW1quasarfile.txt', format='ascii.no_header',
                         names = ['position'])
    
    # KFH Separate RA and dec and convert into degrees
    ra = ['{}h{}m{}s'.format(q[0:2], q[2:4], q[4:9]) for q in quasars['position']]
    dec = ['{}d{}m{}s'.format(q[9:12], q[12:14], q[14:19]) for q in quasars['position']]
    
    coords = SkyCoord(ra, dec, frame='icrs')
    
    quasars.add_columns([coords.ra, coords.dec],  names=['ra', 'dec'])

    # KFH Generate list of days in input month and create times each night
    month_num, day_num = monthrange(2022, month)
    days = np.arange(1, day_num+1)
    
    times = []
    for i in days:
        time = '2022-' + str('%02d' % month) + '-' + str('%02d' % i) + 'T23:00:00.000'
        times.append(Time(time))

    # KFH Find alt and az of quasars at 11 PM MST on each day of the month
    kpno = EarthLocation.of_site('kpno')
    loc_altaz = [coords.transform_to(AltAz(location=kpno, obstime=t)) for t in times]

    # KFH Find airmass of all objects
    airmass = np.asarray([1/np.cos(90 -t.alt.degree) for t in loc_altaz])


    low_am = np.argmin(np.abs(airmass),1)
    low_airmasses =  quasars[low_am]

    # KFH Put all vals into final table of the lowest-airmass objs each night in a month
    final_tab = Table([times, low_airmasses['position'], low_airmasses['ra'], low_airmasses['dec'],np.min(np.abs(airmass),1)],
                  names=['Date', 'Quasar Coordinates', 'RA', 'Dec', 'Airmass'])
    return(final_tab)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Enter an integer corresponding to a month of the year')
    parser.add_argument('month')
    args = parser.parse_args()
    print(best_airmass(args.month))
