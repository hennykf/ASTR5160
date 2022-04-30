import numpy as np
import matplotlib.pyplot as plt
import emcee
from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
from tasks.week11.data_flagging import flux_to_mag
from sklearn import neighbors
import time
from argparse import ArgumentParser

def get_qsos(sources, sweep_folder, rad, saveloc):
	"""
	This function finds the matches to sources in all relevant sweeps files,
	and saves them to a fits file.
	
	Inputs
	------
	sources: Table
	- An astropy table containing the sources to find matches for
	sweep_folder: str
	- The folder where the sweeps files are found
	rad: int
	- The radius in arcseconds around which the matching should be performed
	saveloc: str
	- The full fulename where this file should be saved.
	
	Outputs
	-------
	fits file in location specified by saveloc, containing the sweeps data for
	input sources
	"""
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
		
		# KFH Filter to PSF objects with r<19
		filtered = filedata[filedata['TYPE'] ==b'PSF']
		filtered = filtered[22.5-2.5*np.log10(filtered['FLUX_R']) < 19]
		full_database.append(filtered)   
	full_database = np.concatenate(full_database)

	# KFH turn runtime errors back on
	np.seterr(all='warn')
	
	# KFH Crossmatch objs with sweep database
	objs = SkyCoord(sources['RA'], sources['DEC'], frame='icrs', unit='degree')
	sweep_objs = SkyCoord(full_database["RA"], full_database["DEC"], frame='icrs', unit='degree')
    
	idx1, idx2, sep2d, dist3d = sweep_objs.search_around_sky(objs, seplimit=rad*u.arcsec)
	matches = Table(full_database[idx2])
	
	# KFH Save out sweep data of matching quasars
	matches.write(saveloc, format='fits')


def quasar_class(stars, qsos, unknown):
	"""
	Uses the k-NN algorithm to classify quasars based on g-z and r-W1 colors
	Inputs
	------
	stars: ndarray
	- A set of data with minimal number of quasars, will be treated as training data
	qsos: ndarray
	- Array of confirmed quasars
	unknown: ndarray
	- A random subset of same larger dataset, to be classified
	Returns
	-------
	bools: bool
	- Boolean array of number of objects in 'unknown', True means object was
	identified as a quasar, False means object was ruled out or not identified as a quasar
	star_clean: list
	- List of objects used as stars in knn model
	qsos_clean: list
	- List of confirmed quasars used in knn model
	unknown_clean_eligible: list
	- List of possible quasars categorized by knn
	"""
	# KFH temporarily turn off runtime warnings
	np.seterr(divide='ignore', invalid='ignore')
	
	
	# KFH First calculate g-z and r-W1
	star_color = [flux_to_mag(stars['FLUX_G']) - flux_to_mag(stars['FLUX_Z']),
	 				flux_to_mag(stars['FLUX_R']) - flux_to_mag(stars['FLUX_W1'])]
	qsos_color = [flux_to_mag(qsos['FLUX_G']) - flux_to_mag(qsos['FLUX_Z']),
	 				flux_to_mag(qsos['FLUX_R']) - flux_to_mag(qsos['FLUX_W1'])]
	unknown_color = [flux_to_mag(unknown['FLUX_G']) - flux_to_mag(unknown['FLUX_Z']),
	 				flux_to_mag(unknown['FLUX_R']) - flux_to_mag(unknown['FLUX_W1'])]
	 				
	np.seterr(divide='warn', invalid='warn')

	# KFH reshape and remove nans
	star_clean = [list(a) for a in zip(star_color[0], star_color[1]) if np.isfinite(a).all()==True]
	qsos_clean = [list(a) for a in zip(qsos_color[0], qsos_color[1]) if np.isfinite(a).all()==True]
	
	# KFH Reshape and pull indices without nans
	unknown_clean = [list(a) for a in zip(unknown_color[0], unknown_color[1])]
	unknown_clean_idx = [idx for idx, element in enumerate(unknown_clean) if np.isfinite(element).all() == True]
	unknown_clean_eligible = [unknown_clean[i] for i in unknown_clean_idx]

	# KFH Combine known stars and quasars
	data = np.concatenate([star_clean, qsos_clean])
	data_class = np.concatenate([np.zeros(len(star_clean), dtype='bool'), np.ones(len(qsos_clean), dtype='bool')])
	target_names = np.array(["star", "qso"])
    
    # KFH Use knn to create classifier
	knn = neighbors.KNeighborsClassifier(n_neighbors=1)
	knn.fit(data, data_class)
	
	# KFH Create boolean array of False of length unknown_clean
	bools = np.zeros(len(unknown_clean), dtype='bool')

	# KFH Apply classifier to unknown objects
	predicted_class = knn.predict(unknown_clean_eligible)
	
	# KFH Apply predicted class to bools that could be quasars (doesn't work with a list comprehension btw)
	for (i, b) in zip(unknown_clean_idx, predicted_class):
		bools[i] = b

	return(bools, star_clean, qsos_clean, unknown_clean_eligible)
	
def cut_criteria(table, stars=True):
	"""
	Subsets an astropy Table applying various criteria, to get a subset that are more 
	likely to be stars or quasars than an unsorted sample
	
	Inputs
	------
	table: astropy Table
	- An astropy table of sources
	stars: bool
	- If true, reduces number of objects based on starlike criteria
	- If false, reduces number of objects based on quasarlike criteria
	
	If qsos=True, subsets using criteria that are most likely to apply to quasars
	"""
	# KFH Turn off log warnings
	np.seterr(divide='ignore', invalid='ignore')
	
	# KFH Subset to r<19 psf-type objects
	table = table[22.5-2.5*np.log10(table['FLUX_R']) < 19]
	
	# KFH turn on warnings
	np.seterr(divide='warn', invalid='warn')
	
	table = table[table['TYPE']==b'PSF']
	
	if stars:
		# KFH Subset to only things that are definitely not quasars
		table = table[table['PARALLAX'] > 1]
		table = table[np.abs(table['PMRA']) > 5]
		table = table[np.abs(table['PMDEC']) > 5]
	
	else:
		# KFH Subset to only things that are more likely to be
		table = table[table['PARALLAX'] < 1]
		table = table[np.abs(table['PMRA']) < 5]
		table = table[np.abs(table['PMDEC']) < 5]
	
	return(table)

def splendid_function(table):
	"""
	Takes an astropy Table and classifies the objects, based on a known table
	of quasars and a table of mostly non-quasars. Starlike table and unknown
	table are filtered to remove unlikely candidates for stars and quasars respectively
	to improve categorization.
	"""
	# KFH Read in known quasar data
	quasars = Table.read('quasars_sweep_data.fits', format='fits')

	# KFH Read in some other sources to be used as test "star" data
	stars = Table.read('/d/scratch/ASTR5160/data/legacysurvey/dr9/south/sweep/9.0/sweep-170p030-180p035.fits',
						 format='fits')
	# KFH Clean up test star data to remove likely quasars, and use a subset as test data
	stars = cut_criteria(stars)
	stars = stars[0:200]
	
	# KFH Eliminate things that are definitely NOT quasars
	test_table = cut_criteria(table, stars=False)
	
	# KFH Perform knn fit on criteria-cleared test_table
	binary, star_clean, qsos_clean, unknown_clean = quasar_class(stars, quasars, test_table)
	yes_quasars = test_table[binary]

	# KFH If object is determined to be a quasar, assign True to new column
	table['IS_QSO'] = [True if (i['RA'] in list(yes_quasars['RA'].value)) and 
								(i['DEC'] in list(yes_quasars['DEC'].value)) else
						 False for i in table]
	
	return(np.array(table['IS_QSO']))
	
if __name__=="__main__":

	# KFH The code used to generate quasar file
	# test_data = Table.read('/d/scratch/ASTR5160/week10/qsos-ra180-dec30-rad3.fits', format='fits')
	# get_qsos(test_data, '/d/scratch/ASTR5160/data/legacysurvey/dr9/south/sweep/9.0/', 1,
	#		 'quasars_sweep_data.fits')
	
	#mixed = Table.read('/d/scratch/ASTR5160/data/legacysurvey/dr9/south/sweep/9.0\
	#/sweep-170p025-180p030.fits', format='fits')
	#mixed = mixed[0:20000]

	start  = time.time()
	
	ap = ArgumentParser(description='Find how many quasars are in the given Table')
	ap.add_argument('datatable', help='Fits file name, formated like the sweeps files')
	ns = ap.parse_args()
	
	# KFH Read in file as astropy Table
	tab = Table.read(ns.datatable, format='fits')

	# KFH Run fit on Table and print number of identified quasars
	binary = splendid_function(tab)
	print("Number of quasars", list(binary).count(True))

	end = time.time()
	print("Time elapsed", end-start)


