import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from  astropy.time import Time
from calendar import monthrange

month = 9

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

# KFH Find alt and az of quasars at 11 PM MST
kpno = EarthLocation.of_site('kpno')
loc_altaz = AltAz(location=kpno, obstime = times)
coords.transform_to(loc_altaz)
