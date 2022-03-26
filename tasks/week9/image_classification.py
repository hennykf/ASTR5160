from tasks.week8.cross_matching import sweep_subset
from tasks.week9.mag_systems import sweep_fluxes
import numpy as np
from astropy.table import Table
import astropy.units as u
from astropy.coordinates import SkyCoord
import pandas as pd
import matplotlib.pyplot as plt

def read_fits_files(fits_files, sweep_folder, filter_list):


    # KFH Read in the fits files
    all_tab = []
    for i in fits_files:
        tab = Table.read(i, format='fits')
        all_tab.append(tab['NAME', 'RA', 'DEC'])
    all_tab = np.concatenate(all_tab)

    # KFH Find sweep files where fits files objs might be
    sweep_filenames =  set(['sweep-{}{}{}-{}{}{}.fits'.format(f"{int(q['RA'] - (q['RA']%10)):03d}",
                                                  'm' if (q['DEC'] - (q['DEC']%5)) < 0 else 'p',
                                                   f"{int(q['DEC'] - (q['DEC']%5)):03d}",
                                                   f"{int(q['RA'] + (10-q['RA']%10)):03d}",
                                                   'm' if (q['DEC'] - (q['DEC']%5)) < 0 else 'p',
                                                              f"{int(q['DEC'] + (5-q['DEC']%5)):03d}")
                            for q in all_tab])
    
    sweep_paths = [(sweep_folder + b) for b in sweep_filenames]

    # KFH read in sweeps files
    full_database = []
    for i in sweep_paths:
        filedata = Table.read(i, format='fits')
        full_database.append(filedata)   
    full_database = np.concatenate(full_database)
    print(full_database[0]['RA'])
    
    print("1")
    # KFH Get index of closest object from sweep file
    objects = SkyCoord(all_tab['RA'], all_tab['DEC'], frame='icrs', unit='degree')
    print("2")
    sweep_objs = SkyCoord(full_database["RA"], full_database["DEC"], frame='icrs', unit='degree')
    
    idx1, idx2, sep2d, dist3d = sweep_objs.search_around_sky(objects, seplimit=0.5*u.arcsec)

    # KFH Get columns from database corresponding to filters and indices
    transmission_filters = [('MW_TRANSMISSION_' + f.split('_')[1]) for f in filter_list]
    flux = full_database[idx2][filter_list + transmission_filters]
    
    print(flux.shape)

    # KFH For each object, find transmission-corrected flux in selected filters and calc magnitude
    mag_all = []
    for i in flux:
        flux_filters = [(i[f] / i['MW_TRANSMISSION_' + f.split('_')[1]])for f in filter_list]
        mags = [(22.5 - 2.5*np.log10(f)) for f in flux_filters]
        mag_all.append(mags)
        
    df = pd.DataFrame(mag_all, columns = filter_list)
    return(df)
              


obj_files = ["/d/scratch/ASTR5160/week10/stars-ra180-dec30-rad3.fits",
             "/d/scratch/ASTR5160/week10/qsos-ra180-dec30-rad3.fits"]
sweep_folder = '/d/scratch/ASTR5160/data/legacysurvey/dr9/south/sweep/9.0/'
filter_list = ['FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_W2']
mags = read_fits_files(obj_files, sweep_folder, filter_list)
print(mags['FLUX_G'])

plt.clf()
plt.scatter(mags['FLUX_G'] - mags['FLUX_Z'], mags['FLUX_R'] - mags['FLUX_W1'])
plt.savefig('/d/www/kianah/public_html/week9/stars_quasars.png')
