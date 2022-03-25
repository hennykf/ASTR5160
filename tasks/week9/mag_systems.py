import numpy as np
import os
import pandas as pd
from astropy import units as u
from astropy.coordinates import Angle, search_around_sky, SkyCoord
from astropy.table import Table

def one_pt_sdss_data(ra, dec, textfile):
    """
    Query the sdssDR9 database with one coordinate set and return the results in a txt file

    Inputs
    ------
    ra: str                                                                                                                          
    - RA of object in xxhxxmxxs format                                                                                               
    dec: str                                                                                                                         
    - Declination of object in xxdxxmxxs                                                                                             
    textfile: str                                                                                                                    
    - file name of textfile containing results of sdss query  

    Products
    --------
    textfile saved in current directory containing results of query
    """
    if type(ra)==str:
        ra_deg = Angles(ra).degree
        dec_deg = Angles(dec).degree
    else:
        ra_deg = ra
        dec_deg = dec
        
    # KFH Make sure to get rid of existing file if present
    if (os.path.exists(textfile)):
        os.remove(textfile)

    # KFH Query sdss database and append results of each object query to textfile
    [os.system('python sdssDR9query.py {} {} >> {}'.format(ra_deg, dec_deg, textfile))]

def ubvri_to_ugriz(list_of_filters, ra, dec, textfile):
    """
    Accepts list of magnitudes of filters in the format [V, B-V, U-B, V-R, R-I], and position
    of a star, and returns calculated g and z filter magnitudes for that star.

    Inputs
    ------
    list_of_filters: str
    - List of the magnitudes of the star at the following filters: [V, B-V, U-B, V-R, R-I]
    ra: str                                                                                                                      
    - RA of object in xxhxxmxxs format                                                                                               
    dec: str                                                                                                                         
    - Declination of object in xxdxxmxxs
    textfile: str
    - file name of textfile containing sdss query

    Returns
    -------
    list  containing:
    - g filter magnitude, calculated 
    - z filter magnitude calculated
    - g filter magnitude from sdss query
    - z filter magnitude from sdss query
    """

    if list_of_filters[4] < 1.15:
        # KFH Calculate ugriz filters from uvbri
        ug = 1.28*list_of_filters[2] + 1.13
        gr = 1.02*list_of_filters[1] - 0.22
        ri = 0.91*list_of_filters[4] - 0.20
        rz = 1.72*list_of_filters[4] - 0.41
        g = list_of_filters[0] + 0.60*list_of_filters[1] - 0.12
        rr = list_of_filters[0] - 0.42*list_of_filters[1] + 0.11
        z = -(rz - rr)
        
    else:
        print("Out of color range")

    # KFH query sdss with point
    one_pt_sdss_data(ra, dec, textfile)

    # KFH read in resulting file
    data = pd.read_csv(textfile, sep=',', header=None, index_col=False)
    data.columns = ['ra', 'dec', 'u', 'g', 'r', 'i', 'z', 'sep']

    return([g,z,data['g'][0], data['z'][0]])


def sweep_fluxes(ra, dec, list_of_cols, type="mag"):
    """
    Find sweep file that contains the object at the ra, dec, and calculate
    magnitude in given filters.

    Inputs
    ------
    ra: str
    - RA of object in xxhxxmxxs format
    dec: str
    - Declination of object in xxdxxmxxs
    list_of_cols: list
    - A list of the column names in the sweep file.
      Info on columns can be found at https://www.legacysurvey.org/dr7/files/

    Returns
    -------
    mag: list
    - A list of the magnitudes of the object at the filters given.
    """

    ra_deg = Angle(ra).degree
    dec_deg = Angle(dec).degree
        
    # KFH Find filename where object is present
    filename = 'sweep-{}{}{}-{}{}{}.fits'.format(f"{int(ra_deg - (ra_deg%10)):03d}",
                                                 'm' if (dec_deg - (dec_deg%5)) < 0 else 'p',
                                                 f"{int(dec_deg - (dec_deg%5)):03d}",
                                                 f"{int(ra_deg + (10-ra_deg%10)):03d}",
                                                 'm' if (dec_deg - (dec_deg%5)) < 0 else 'p',
                                                  f"{int(dec_deg + (5-dec_deg%5)):03d}")
    path = '/d/scratch/ASTR5160/data/legacysurvey/dr9/south/sweep/9.0/' + str(filename)
    dat = Table.read(path)

    # KFH Get index of closest object from sweep file
    obj1 = SkyCoord(np.array([ra_deg]), np.array([dec_deg]), frame='icrs', unit='degree')
    file_objs = SkyCoord(dat["RA"], dat["DEC"], frame='icrs', unit='degree')
    idx1, idx2, sep2d, dist3d = file_objs.search_around_sky(obj1, seplimit=1*u.arcsec)

    if type == "mag":
        # KFH Compute magnitudes of all filters given
        mag = [(22.5 - 2.5*np.log10(dat[idx2][c]))[0] for c in list_of_cols]

    elif type=="flux":
        # KFH Return raw fluxes for listed columns
        mag = [dat[idx2][c] for c in list_of_cols]

    else:
        print("Type must be 'mag' or 'flux'")

        pass
    
    return(mag)

if __name__=="__main__":
    # KFH convert from ubvri to ugriz
    conversions = ubvri_to_ugriz([15.256, 0.873, 0.320, 0.505, 0.511],
                                 ra = "16h35m26s",
                                 dec = "+09d47m53s",
                                 textfile = "pg1633_099a.txt")

    # Get magnitudes of objs from sweeps
    sweeps = sweep_fluxes("16h35m26s", "+09d47m53s", ['FLUX_G', 'FLUX_R', 'FLUX_Z'])
    sweeps_wise = sweep_fluxes("16h35m26s", "+09d47m53s",
                               ['FLUX_W1', 'FLUX_W2', 'FLUX_W3', 'FLUX_W4'])
    
    # KFH These roughly correspond to the navigate tool values
    print("g and z calculated, g and z from sdss: " + str(conversions))
    print("g, r, and z magnitudes from sweep: " + str(sweeps))

    # KFH We can see that this star was not detected in the w4 band.
    print("w1, w2, w3, w4 magnitudes from sweep: " + str(sweeps_wise))
