from astropy.table import Table
import  matplotlib.pyplot as plt
import os

def retrieve_sdss_data(fits_file, textfile):
    """
    Query the sdssDR9 database with the coordinates of the first 100
    objects in the fits_file.

    Inputs
    ------
    fits_file: str
    - The full path name of the file with object coordinates. Must have
    an "RA" column and a "DEC" column
    textfile: str
    - The name of the textfile where the results of the querys will be
    saved. The file will be generated in the current working directory.

    Outputs
    -------
    A file containing the results of the sdss query for each of the first
    100 objects in the fits_file. Some objects may not appear in the sdss
    database, and if so, the line will read "No objects have been found"
    """

    # KFH Read in the fits file and take only the first 100 objs
    tab = Table.read(fits_file)
    first_100 = tab[0:101]

    # KFH Remove file if it exists so it doesn't append an existing file
    if (os.path.exists(textfile)):
        os.remove(textfile)

    # KFH Query sdss database and append results of each object query to textfile
    [os.system('python sdssDR9query.py {} {} >> {}'.format(q['RA'], q['DEC'],
                                                           textfile)) for q in first_100]


def sweep_subset(fits_file):
    """
    Function to find the sweep files in the directory 
    '/d/scratch/ASTR5160/data/legacysurvey/dr9/north/sweep/9.0'
    encapsulate the coordinates of each of the first 100 items in the fits_file.

    Inputs
    ------
    fits_file: str
    - Location of a fits file, must contain colums "RA" and "DEC".

    Returns
    -------
    filenames: list
    - list of all the sweep file names that  cover the
    10 deg RA x 5 deg DEC box a given object would be found in.
    """

    # KFH Get the first 100 objects of the fits file
    tab = Table.read(fits_file)
    first_100 = tab[0:101]
    
    # KFH Find the names of all the files where coordinates would be
    filenames = ['sweep-{}p{}-{}p{}.fits'.format(int(q['RA'] - (q['RA']%10)),
                                                 int(q['DEC'] - (q['DEC']%5)),
                                                 int(q['RA'] + (10-q['RA']%10)),
                                                 int(q['DEC'] + (5-q['DEC']%5))) for q in first_100]

    # KFH Return the set of unique filenames where we can find these objects
    return(set(filenames))


if __name__=='__main__':

    # KFH Plot all objects in the first file
    all_dat = Table.read('/d/scratch/ASTR5160/data/first/first_08jul16.fits')

    plt.clf()
    plt.scatter(all_dat['RA'], all_dat['DEC'], s=.01)
    plt.savefig('/d/www/kianah/public_html/week8/vla_first_pts.png')
    plt.xlabel('RA')
    plt.ylabel('Dec')

    # KFH Query sdss database and save results in first_100.txt
    retrieve_sdss_data('/d/scratch/ASTR5160/data/first/first_08jul16.fits', 'first_100.txt')

    # KFH Find the names of all the sweep files where we can find the first 100 objs in first
    files = sweep_subset('/d/scratch/ASTR5160/data/first/first_08jul16.fits')
    print(files)
