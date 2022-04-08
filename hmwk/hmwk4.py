import  numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
import time
from tasks.week8.sdssDR9query.py import sdssDR9query


def within_radius(ra_center, dec_center, rad, fits_file):
	"""
	Function to find what objects in fits_file are within
	rad of ra_center and dec_center
	
	Inputs
	------
	ra_center: float
	- The RA of the center of the spherical cap
	dec_center: float
	- The declination of the center of the spherical cap
	rad: float
	- The radius of the spherical cap in degrees
	fits_file: str
	- full path and filename of fits file to be filtered
	
	Returns
	-------
	objs_in_rad: ndarray
	- The table of objects from fits_file that fall within the
	  radius of the spherical cap
	"""
	
	# KFH Read in fits file and create cartesian coords column
	tab = Table.read(fits_file)
	tab_locs = SkyCoord(tab['RA'], tab['DEC'], frame='icrs', unit='degree')
	
	# KFH Find objects within radius of ra_center, dec_center
	ra_center = np.asarray([ra_center])
	dec_center = np.asarray([dec_center])
	center = SkyCoord(ra_center, dec_center, frame='icrs', unit='degree')
	idx1, idx2, sep2d, dist3d = tab_locs.search_around_sky(center, seplimit=rad*u.deg)
	
	objs_in_rad = tab[idx2]

	return(objs_in_rad)
	


def crossmatch(sources, sweep_folder, rad):
	"""
	Crossmatch sources array with sweep_path sources
	
	"""
	# KFH Find sweep files where fits files objs might be
	sweep_filenames =  set(['sweep-{}{}{}-{}{}{}.fits'.format(f"{int(q['RA'] - (q['RA']%10)):03d}",
	 'm' if (q['DEC'] - (q['DEC']%5)) < 0 else 'p',
	 f"{int(q['DEC'] - (q['DEC']%5)):03d}", 
	 f"{int(q['RA'] + (10-q['RA']%10)):03d}", 
	 'm' if (q['DEC'] - (q['DEC']%5)) < 0 else 'p', 
	 f"{int(q['DEC'] + (5-q['DEC']%5)):03d}") for q in sources])
    
	sweep_paths = [(sweep_folder + b) for b in sweep_filenames]

	# KFH read in sweeps files
	full_database = []
	for i in sweep_paths:
		filedata = Table.read(i, format='fits')
		full_database.append(filedata)   
	full_database = np.concatenate(full_database)
	
	# KFH Crossmatch objs with sweep database
	objs = SkyCoord(sources['RA'], sources['DEC'], frame='icrs', unit='degree')
	sweep_objs = SkyCoord(full_database["RA"], full_database["DEC"], frame='icrs', unit='degree')
    
	idx1, idx2, sep2d, dist3d = sweep_objs.search_around_sky(objs, seplimit=rad*u.arcsec)
	matches = full_database[idx2]
	
	# KFH Filter to PSF objects with r<22 and W1-W2>0.5
	filtered = matches[matches['TYPE'] ==b'PSF']
	filtered = filtered[22.5-2.5*np.log10(filtered['FLUX_R']) < 22]
	filtered = filtered[(22.5-2.5*np.log10(filtered['FLUX_W1'])) -
	(22.5-2.5*np.log10(filtered['FLUX_W2'])) > 0.5]

	return(filtered)
	
def retrieve_sdss(datatable):
	"""
	Retrieve data from sdss for datatable objects
	"""

	url='http://skyserver.sdss3.org/dr9/en/tools/search/x_sql.asp'

	sdssQuery(datatable['RA'][0], datatable['DEC'][1])

	# KFH Create query for sdss
	all_queries = []
	for i in datatable:
		query = """SELECT top 1 ra,dec,u,g,r,i,z,GNOE.distance*60 FROM PhotoObj as PT
    	JOIN dbo.fGetNearbyObjEq(""" + str(i['RA']) + """,""" + str(i['DEC']) + """,0.02) as GNOE
    	on PT.objID = GNOE.objID ORDER BY GNOE.distance"""
		all_queries.append(query)

	
	print(all_queries)

if __name__=="__main__":

	t1 = time.time()

	# KFH Create spherical cap and find objs in cap
	objs_in_rad = within_radius(163., 50., 3.,
	 fits_file = '/d/scratch/ASTR5160/data/first/first_08jul16.fits')
	
	# KFH Crossmatch objs with sweep objs
	filtered_objs = crossmatch(objs_in_rad, '/d/scratch/ASTR5160/data/legacysurvey/dr9/north/sweep/9.0/', 1)
	print("Number of filtered objects:", len(filtered_objs))
	
	# KFH Retrieve sdss data
	retrieve_sdss(filtered_objs)
	
	t2=time.time()
	print('Time', t2-t1)
