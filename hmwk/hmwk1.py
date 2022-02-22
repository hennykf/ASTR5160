import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
import astropy.units as u
from astropy.time import Time
from calendar import monthrange
import argparse

# ADM just hard-coding the offset between UTC and MST...
utc_off = 7*u.hour
# ADM ...and the filename.
filename = "/d/scratch/ASTR5160/week4/HW1quasarfile.txt"

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
    quasars = Table.read(filename, format='ascii.no_header', names=['position'])
    
    # KFH Separate RA and dec and convert into degrees
    ra = ['{}h{}m{}s'.format(q[0:2], q[2:4], q[4:9]) for q in quasars['position']]
    dec = ['{}d{}m{}s'.format(q[9:12], q[12:14], q[14:19]) for q in quasars['position']]
    
    coords = SkyCoord(ra, dec, frame='icrs')
    
    quasars.add_columns([coords.ra, coords.dec],  names=['ra', 'dec'])

    # KFH Generate list of days in input month and create times each night
    month_num, day_num = monthrange(2022, month)
    days = np.arange(1, day_num+1)
    
    # ADM More compact way of constructing times.
    t = ['2022-' + str(month) + '-' + str(i) + ' 23:00:00' for i in days]
    times = Time(t, format='iso') + utc_off

    # KFH Find alt and az of quasars at 11 PM MST on each day of the month
    kpno = EarthLocation.of_site('kpno')

    # ADM find the airmasses at the given times in just one step.
    airmass = np.vstack(
        [coords.transform_to(AltAz(location=kpno, obstime=t)).secz.value
         for t in times])

    # ADM find the index of the lowest airmass, first ensuring that 
    # ADM negative numbers are converted to a very-high airmass.
    ii = np.argmin(np.where(airmass > 0, airmass, 1e16), 1)
    low_am = np.min(np.where(airmass > 0, airmass, 1e16), 1)
    low_qs = quasars[ii]

    # KFH Put all vals into final table of the lowest-airmass objs each night in a month
    # ADM just to be illustrative, I changed the times back into MST instead of UTC.
    final_tab = Table(
        [times-utc_off, low_qs['position'], low_qs['ra'], low_qs['dec'], low_am],
        names=['Date (MST)', 'Quasar Coordinates', 'RA', 'Dec', 'Airmass'])

    return(final_tab)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Enter an integer corresponding to a month of the year')
    parser.add_argument('month', type=int)
    args = parser.parse_args()
    print(best_airmass(args.month))
