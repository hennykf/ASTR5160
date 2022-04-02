import numpy as np
from tasks.week9.mag_systems import  sweep_fluxes
from astropy.coordinates import Angle, search_around_sky, SkyCoord
from astropy.table import Table
from astropy import units as u
from astropy.visualization.wcsaxes import SphericalCircle
import matplotlib.pyplot as plt
from sklearn import neighbors
import math

def sweep_obj(ra, dec, rad=.05, objtype=False):
    """
    Find sweep file that contains the object at the ra, dec, and calculate
    magnitude in given filters.

    Inputs
    ------
    ra: str
    - RA of object in xxhxxmxxs format
    dec: str
    - Declination of object in xxdxxmxxs
    rad: float
    - the radius of the circle to search in arcseconds
    objtype: bool
    - if true, subsets to only 'PSF' and r<20 magnitude objects

    Returns
    -------
    obj: ndarray
    - an array of all objects found to be within the rad of the ra, dec.
    """

    # KFH Handle both hms format and degree format
    if type(ra)==str:
        ra_deg = Angle(ra).degree
        dec_deg = Angle(dec).degree

    else:
        ra_deg = ra
        dec_deg = dec

    # KFH Generate points in circle to get all sweeps files needed
    circle = SphericalCircle([SkyCoord(ra_deg, dec_deg, frame='icrs', unit='degree').ra,
                              SkyCoord(ra_deg, dec_deg, frame='icrs', unit='degree').dec], (rad/3600)*u.deg)

    circle_coords = circle.get_xy()
    ra_circle, dec_circle = [list(x) for x in zip(*circle_coords)]
    coords = Table([ra_circle, dec_circle], names=('ra', 'dec'))
        
    # KFH Find filename where object is present
    filename = set(['sweep-{}{}{}-{}{}{}.fits'.format(f"{int(q['ra'] - (q['ra'] % 10)):03d}",
                                                      'm' if (q['dec'] - (q['dec'] %5)) < 0 else 'p',
                                                      f"{int(q['dec'] - (q['dec'] % 5)):03d}",
                                                      f"{int(q['ra'] + (10-q['ra'] % 10)):03d}",
                                                      'm' if (q['dec'] - (q['dec'] %5)) < 0 else 'p',
                                                      f"{int(q['dec'] + (5-q['dec'] % 5)):03d}") for q in coords])

    
    path = [('/d/scratch/ASTR5160/data/legacysurvey/dr9/south/sweep/9.0/' + f) for f in filename]
    
    # KFH If multiple sweeps files, combine data from all sweeps files
    full_dat = []
    for i in path:
        dat = Table.read(i, format='fits')

        # KFH Subset to only psf and mag <20 objs
        if objtype:
            dat = dat[dat['TYPE']== b'PSF']
            dat = dat[(22.5-2.5*np.log10(dat['FLUX_R'])) < 20]
        full_dat.append(dat)

    full_dat = np.concatenate(full_dat)

    # KFH Get index of closest object from sweep file
    obj1 = SkyCoord(np.array([ra_deg]), np.array([dec_deg]), frame='icrs', unit='degree')
    file_objs = SkyCoord(full_dat["RA"], full_dat["DEC"], frame='icrs', unit='degree')
    idx1, idx2, sep2d, dist3d = file_objs.search_around_sky(obj1, seplimit=rad*u.arcsec)

    obj = full_dat[idx2]
    return(obj)

def coord_match(datatable, fits_file):
    # KFH Finds matches of fits file in datatable, returns the datatable rows that match
    tab = SkyCoord(datatable["RA"], datatable["DEC"], frame='icrs', unit='degree')
    fits = Table.read(fits_file, format='fits')
    fits2 = SkyCoord(fits['RA'], fits['DEC'], frame='icrs', unit='degree')
    idx1, idx2, sep2d, dist3d = tab.search_around_sky(fits2, seplimit=.5*u.arcsec)

    return(datatable[idx2])

def quasar_class(stars, qsos, unknown):
    """
    Uses the k-NN algorithm to classify quasars based on g-z and r-W1 colors

    Inputs
    ------
    stars: ndarray
    - A random subset of a larger dataset, will be treated as training data
    qsos: ndarray
    - Array of confirmed quasars
    unknown: ndarray
    - A random subset of same larger dataset, to be classified

    Returns
    -------
    
    """
    # KFH First calculate g-z and r-W1
    star_color = [flux_to_mag(stars['FLUX_G']) - flux_to_mag(stars['FLUX_Z']),
                   flux_to_mag(stars['FLUX_R']) - flux_to_mag(stars['FLUX_W1'])]
    qsos_color = [flux_to_mag(qsos['FLUX_G']) - flux_to_mag(qsos['FLUX_Z']),
                   flux_to_mag(qsos['FLUX_R']) - flux_to_mag(qsos['FLUX_W1'])]
    unknown_color = [flux_to_mag(unknown['FLUX_G']) - flux_to_mag(unknown['FLUX_Z']),
                      flux_to_mag(unknown['FLUX_R']) - flux_to_mag(unknown['FLUX_W1'])]

    # KFH reshape and remove nans
    star_clean = [list(a) for a in zip(star_color[0], star_color[1]) if np.isfinite(a).any()==True]
    qsos_clean = [list(a) for a in zip(qsos_color[0], qsos_color[1]) if np.isfinite(a).any()==True]
    unknown_clean = [list(a) for a in zip(unknown_color[0], unknown_color[1]) if np.isfinite(a).any()==True]
    
    data = np.concatenate([star_color, qsos_color])
    data_class = np.concatenate([np.zeros(len(stars), dtype='i'),
                                 np.ones(len(qsos), dtype='i')])
    target_names = np.array(["star", "qso"])
    
    knn = neighbors.KNeighborsClassifier(n_neighbors=1)
    knn.fit(data, data_class)

    predicted_class = knn.predict(unknown)

    
def flux_to_mag(flux):
    # KFH Convert a flux to magnitude
    mag = 22.5-2.5*np.log10(flux)
    return(mag)

if __name__=="__main__":

    # KFH The type is exponential, indicative of a non-stellar source
    obj = sweep_obj(188.53667, 21.04572, rad=.05)

    # KFH All values are 2, so they are saturated
    print("The object is saturated, with mask G, R, Z band values of: ", obj['ALLMASK_G'], obj['ALLMASK_R'], obj['ALLMASK_Z'])

    # KFH This object on the sky viewer is definitely not a galaxy
    # KFH According to simbad it's a possible blazar. This would explain why
    # KFH it's so blue, since blazars are pointing at us. So it's not a galaxy itself
    # KFH but represents a galaxy.

    psfobjs = sweep_obj(180,30, 10800, objtype = 'PSF')
    print(len(psfobjs))

    qsos = coord_match(psfobjs, '/d/scratch/ASTR5160/week10/qsos-ra180-dec30-rad3.fits')
    print(len(qsos))

    # KFH To avoid copying code, subsequent code is part of day 2 tasks.
    print("All subsequent info is now part of day 2 tasks")
    training_data = psfobjs[np.random.randint(psfobjs.shape[0], size=275)]
    test_data = psfobjs[np.random.randint(psfobjs.shape[0], size=275)]
    quasar_class(training_data, qsos, test_data)
    
    
