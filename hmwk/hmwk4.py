import  numpy as np
from astropy.table import Table
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import astropy.units as u
import time
import sdssDR9query
from tasks.week11.data_flagging import flux_to_mag
from argparse import ArgumentParser


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
	
	# KFH Subset to only objects within radius
	objs_in_rad = tab[idx2]

	return(objs_in_rad)
	


def crossmatch(sources, sweep_folder, rad):
	"""
	Crossmatch sources array with sweep_folder sources within given radius

	Inputs
	------
	sources: ndarray
	- Table of objects to be matched, eg. output array from 'within_radius' function
	sweep_folder: str
	- path of folder containing sweep files
	rad: float
	- Radius to crossmatch sources in 'sources' and in sweep files, in arcseconds
	
	"""
	# KFH Find sweep files where fits files objs might be
	sweep_filenames =  set(['sweep-{}{}{}-{}{}{}.fits'.format(f"{int(q['RA'] - (q['RA']%10)):03d}",
	 'm' if (q['DEC'] - (q['DEC']%5)) < 0 else 'p',
	 f"{int(q['DEC'] - (q['DEC']%5)):03d}", 
	 f"{int(q['RA'] + (10-q['RA']%10)):03d}", 
	 'm' if (q['DEC'] - (q['DEC']%5)) < 0 else 'p', 
	 f"{int(q['DEC'] + (5-q['DEC']%5)):03d}") for q in sources])
    
	sweep_paths = [(sweep_folder + b) for b in sweep_filenames]

	# KFH Ignore runtime errors caused by logrithms
	np.seterr(all='ignore')

	# KFH read in sweeps files and subset
	full_database = []
	for i in sweep_paths:
		filedata = Table.read(i, format='fits')
		
		# KFH Filter to PSF objects with r<22 and W1-W2>0.5
		filtered = filedata[filedata['TYPE'] ==b'PSF']
		filtered = filtered[22.5-2.5*np.log10(filtered['FLUX_R']) < 22]
		filtered = filtered[(22.5-2.5*np.log10(filtered['FLUX_W1'])) -
		(22.5-2.5*np.log10(filtered['FLUX_W2'])) > 0.5]
		full_database.append(filtered)   
	full_database = np.concatenate(full_database)

	# KFH turn runtime errors back on
	np.seterr(all='warn')
	
	# KFH Crossmatch objs with sweep database
	objs = SkyCoord(sources['RA'], sources['DEC'], frame='icrs', unit='degree')
	sweep_objs = SkyCoord(full_database["RA"], full_database["DEC"], frame='icrs', unit='degree')
    
	idx1, idx2, sep2d, dist3d = sweep_objs.search_around_sky(objs, seplimit=rad*u.arcsec)
	matches = full_database[idx2]
	return(matches)
	
def retrieve_sdss(obj):
	"""
	Retrieve data from sdss for datatable objects. Function based off Adam Myers' code

	Input
	-----
	obj: ndarray
	- Object with columns 'RA' and 'DEC' to be queried in sdss database

	Returns
	-------
	res: list
	- List of ra, dec, u, i from sdss database
	"""

	# KFH Get ra and dec
	ra = obj['RA']
	dec = obj['DEC']

	# KFH Initialize class object
	qry = sdssDR9query.sdssQuery()
	
	# KFH Set query string
	query = """SELECT top 1 ra, dec, u,i FROM PhotoObj as PT \
		JOIN dbo.fGetNearbyObjEq(""" + str(ra) + """,""" + str(dec) + """,0.02) \
		as GNOE on PT.objID = GNOE.objID ORDER BY GNOE.distance"""

	qry.query = query
	for line in qry.executeQuery():
		result = line.strip()

	# KFH Don't overwhelm the servers
	time.sleep(1)

	# KFH Split query results into separate rows
	res = result.decode().split(',')
	return(res)

def mag_to_flux(mag):
	"""
	Converts a magnitude value, mag, to a flux in nanomaggies.
	Returns flux in nanomaggies
	"""
	flux = 10**((22.5 - mag) / 2.5)
	return(flux)

if __name__=="__main__":
	
	# KFH Begin timer
	t1 = time.time()

	ap = ArgumentParser(description='Find SDSS u-band properties of Legacy Surveys \
				point sources that have an IR excess and are detected in \
				the radio by the FIRST survey')
	ap.add_argument('ra', help='Right Ascension (degrees)', type=float)
	ap.add_argument('dec', help='Declination (degrees)', type=float)
	ap.add_argument('radius', help='Radius of spherical cap (degrees)', type=float)

	ns = ap.parse_args()
	# KFH Create spherical cap and find objs in cap
	objs_in_rad = within_radius(ns.ra, ns.dec, ns.radius, \
	 fits_file = '/d/scratch/ASTR5160/data/first/first_08jul16.fits')
	
	# KFH Crossmatch objs with sweep objs
	filtered_objs = crossmatch(objs_in_rad, '/d/scratch/ASTR5160/data/legacysurvey/dr9/north/sweep/9.0/', 1)
	print("Number of filtered objects:", len(filtered_objs))
	
	# KFH Retrieve sdss data
	sdss_objs = []
	for i in filtered_objs:
		data = retrieve_sdss(i)
		sdss_objs.append(data)

	# KFH Remove any entries where no objects were matched in sdss
	sdss_objs.remove(['No objects have been found'])
	print("Number of sdss matches:", len(sdss_objs))

	# KFH Turn ra, dec, u and i into columns to be matched	
	sdss_objs = np.array(sdss_objs).T.astype('float')
	sdss_tab = Table([sdss_objs[0], sdss_objs[1], sdss_objs[2], sdss_objs[3]],
			 names=('RA', 'DEC', 'u', 'i'))
	
	# KFH Find the brightest source in u band and calc flux of u and i
	ubrite1 = sdss_tab[sdss_tab['u'] == min(sdss_tab['u'])]
	ubrite1.add_columns([[mag_to_flux(ubrite1['u'])], [mag_to_flux(ubrite1['i'])]],
				 names=['u_flux', 'i_flux'])

	# KFH Match back with the filtered_obj that matches this ra and dec
	ub = SkyCoord(ubrite1['RA'], ubrite1['DEC'], frame='icrs', unit='degree')
	fo = SkyCoord(filtered_objs['RA'], filtered_objs['DEC'], frame='icrs', unit='degree')
	idx1, idx2, sep2d, dist3d = fo.search_around_sky(ub, seplimit=(.02/60)*u.deg)

	# KFH Match is the filtered_objs object with all columns of original fits file
	# KFH SDSS info about this obj can be found here:
	# https://skyserver.sdss.org/dr9/en/tools/explore/obj.asp?id=1237658206115201095
	match = filtered_objs[idx2]

	print("Info about ubrite object")
	print("This is a quasar, at a z=1.035. Quasars are very blue, so it's unsurprising \
 		that the brightest object in the survey at u band would be a quasar. Given that \
		quasars are powered by AGN, it is also not surprisng that it is at a high-ish \
		redshift, as there are a lot more AGN at z~1 than in the local universe")
	
	# KFH Plot wavelength vs flux for ubrite source
	fluxes = [ubrite1['u_flux'], match['FLUX_G'], match['FLUX_R'], ubrite1['i_flux'],
		 match['FLUX_Z'], match['FLUX_W1'], match['FLUX_W2'], match['FLUX_W3'], match['FLUX_W4']]
	freqs = [3543e-10, 4770e-10, 6231e-10, 7625e-10, 9134e-10, 3.4e-6, 4.6e-6, 12e-6, 22e-6] 

	plt.plot(freqs, fluxes, 'bo', markersize=3, linestyle='dashed')
	plt.xlabel('Wavelength (m)')
	plt.ylabel('Flux (nanomaggies)')
	plt.title('Flux vs wavelength, source at RA={}, dec={}'.format(str(ubrite1['RA'].value),
								str(ubrite1['DEC'].value)))

	# KFH End timer
	t2=time.time()
	print('Time', t2-t1)
	# KFH Because the timer ends after you close the plot.show() window, ending timer before showing
	plt.show()
