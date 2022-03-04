import sfdmap
import numpy as np
from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt

def ugriz_extinction(ra, dec):
    dustdir = '/d/scratch/ASTR5160/data/dust/v0_1/maps'
    m = sfdmap.SFDMap(dustdir, scaling=1)

    c = SkyCoord(ra, dec)
    ebv = m.ebv(c.ra.value, c.dec.value)
    ugriz = np.array([4.239 ,3.303 , 2.285 , 1.698 , 1.263])
    A = ebv*ugriz
    return A

def quasar_color(correction = False):
    quasar1 = [18.82, 18.81, 18.73, 18.82, 18.90]
    quasar2 = [19.37, 19.10, 18.79, 18.73, 18.63]

    plotloc = '/d/www/kianah/public_html/week3_quasars.png'

    if correction:
        a1 = ugriz_extinction('16h27m43.9s', '40d03m10.8s')
        a2 = ugriz_extinction('15h46m14.8s', '2d01m45.5s')
        quasar1 = quasar1 - a1
        quasar2 = quasar2 - a2
        
        plotloc = '/d/www/kianah/public_html/week3_quasars_ext.png'

    plt.clf()
    plt.scatter(quasar1[1]-quasar1[2], quasar1[2]-quasar1[3], label='246.933, 40.795')
    plt.scatter(quasar2[1]-quasar2[2], quasar2[2]-quasar2[3], label='236.562, 2.44')
    plt.legend()
    plt.xlim([.05,.35])
    plt.ylim([-.1,.1])
    plt.xlabel('g-r')
    plt.ylabel('r-i')
    plt.title('Quasar colors')
    plt.savefig(plotloc)

quasar_color()
quasar_color(correction=True)
# The quasars do not have similar colors, which makes sense because they peak in different wavelengths. Quasar 1 is brightest in i, and quasar 2 is brighter in z band.
# After a magnitude correction, the quasars are much closer in color.

def dustviz():
    ra1 = np.linspace(231.6, 241.6, num=100, endpoint=False)
    dec1 = np.linspace(-2.6, 7.4, num=100, endpoint=False)

    ra2 = np.linspace(240.4, 253.4, num=100, endpoint=False)
    dec2 = np.linspace(35.8, 45.8, num=100, endpoint=False)

    grid1 = np.meshgrid(ra1, dec1)
    grid2 = np.meshgrid(ra2, dec2)

    print(grid1, grid2)
dustviz()



ugriz_extinction('12h42m30s', '+41d12m00s')
