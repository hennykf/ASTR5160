import pandas as pd
import matplotlib.pyplot as plt

def plot_sdss(filename, savepath, size_var=False):
    """
    Function to plot RA and dec of some sdss sources
    Inputs
    ------
    filename: str
    - path to csv of source points
    savepath: str
    - path to desired save location of plot
    size_var: bool
    - If False, all points will be size 1.
    - If True, points will be size 5/source mag at g band
    """
    
    # KFH Import datatable
    filename = 'sdss_query.csv'
    objs = pd.read_csv(filename)

    plt.clf()
    if size_var:
        # KFH Make larger mags smaller points
        sz=20000/(objs['g']**3)
    else:
        # KFH All points same size
        sz=5
    
    # KFH Plot scatter with same size points
    plt.scatter(objs['ra'], objs['dec'], marker='o', s=sz)
    plt.xlabel('RA')
    plt.ylabel('Dec')
    plt.title('Positions of points within 2 arcmin of [300,-1]')
    plt.ticklabel_format(useOffset=False)
    plt.savefig(savepath)

if __name__=='__main__':
    # KFH Plot  same size points and differing size points 
    plot_sdss('sdss_query.csv', '/d/www/kianah/public_html/week8/sdss_ra_dec.png')

    plot_sdss('sdss_query.csv', '/d/www/kianah/public_html/week8/sdss_ra_dec_varsize.png', size_var=True)

          

