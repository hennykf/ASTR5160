import  numpy as np
from tasks.week6.spherical_caps import cap
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u


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
	
	# KFH get 4vect of spherical cap
	cap1 = cap(ra_center*u.deg, dec_center*u.deg, rad*u.deg)
	
	# KFH Read in fits file and create cartesian coords column
	tab = Table.read(fits_file)
	print("Starting num objs ", len(tab))
	tab.add_column(SkyCoord(tab['RA'], tab['DEC'], frame='icrs', unit='degree').cartesian,
			name='cart_coords')
	
	# KFH Subset to only objects within rad of ra_center, dec_center
	objs_in_rad = tab[tab['cart_coords'].x.value**2 + \
				tab['cart_coords'].y.value**2 + \
				tab['cart_coords'].z.value <= (-1-cap1[3])**2]
	print("FIRST sources within survey ", len(objs_in_rad))
	
	return(objs_in_rad)
	


def crossmatch(sources, sweep_folder, rad, filt):
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
	

if __name__=="__main__":

	# KFH Create spherical cap
	within_radius(163., 50., 3., fits_file = '/d/scratch/ASTR5160/data/first/first_08jul16.fits')
	

