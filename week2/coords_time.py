import  astropy.coordinates
import astropy.time
import numpy as np

# Calc some coordinates
coords = astropy.coordinates.SkyCoord('01h35m30s', '+15d15m22s', frame='icrs')

# Calc same coords manually
dec = 1*(15+(15/60) + (22/3600))
ra = 15*(1+(35/60)+(30/3600))

print(coords)
print(ra, dec)

# Find MJD and JD
now_jd = astropy.time.Time.now().jd
now_mjd = astropy.time.Time.now().mjd

# Check MJD manually
test_mjd = now_jd - 2400000.5

print(now_mjd, test_mjd)

# Generate a bunch of days around the current MJD
some_days = np.arange(now_mjd-5, now_mjd+5)
print(some_days)
