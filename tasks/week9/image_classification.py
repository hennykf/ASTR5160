from tasks.week8.cross_matching import sweep_subset
from tasks.week9.mag_systems import sweep_fluxes
import numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord

def read_fits_files(fits_files, sweep_folder, filter_list):

    # KFH Find sweep files where objs might be
    sweep_filenames = sweep_subset(fits_files[0])
    sweep_paths = [(sweep_folder + b) for b in sweep_filenames]

    # KFH Read in the tables
    tab = [Table.read(f, format='fits') for f in fits_files]
    database = [Table.read(b) for b in sweep_paths]

    # KFH Get index of closest object from sweep file
    objects = [[SkyCoord(i['RA'], i['DEC'], frame='icrs', unit='degree') for i in tb] for tb in tab]
    sweep_objs = [[SkyCoord(i["RA"], i["DEC"], frame='icrs', unit='degree') for i in db] for db in database]
    idx1, idx2, sep2d, dist3d = sweep_objs.search_around_sky(objects, seplimit=0.5*u.arcsec)

    flux = [sweep_objs[idx2][c] for c in filter_list]

    print(flux[0:10])

    flux2  = [(flux[f] / flux['MW_TRANSMISSION_' + f.split('_')[1]])for f in filter_list]
    mags = [(22.5 - 2.5*np.log10(f)) for f in flux2]

    print(mags[0:10]
              


obj_files = ["/d/scratch/ASTR5160/week10/stars-ra180-dec30-rad3.fits",
             "/d/scratch/ASTR5160/week10/qsos-ra180-dec30-rad3.fits"]
sweep_folder = '/d/scratch/ASTR5160/data/legacysurvey/dr9/south/sweep/9.0/'
filter_list = ['FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_W2',
               'MW_TRANSMISSION_G', 'MW_TRANSMISSION_R', 'MW_TRANSMISSION_Z',
               'MW_TRANSMISSION_W1', 'MW_TRANSMISSION_W2']
read_fits_files(obj_files, sweep_folder, filter_list)
