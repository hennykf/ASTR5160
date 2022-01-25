import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

objs = Table.read("/d/scratch/ASTR5160/week2/struc.fits")
plt.scatter(objs["RA"], objs["DEC"], s=2, marker = "o", color="blue")
ii = objs["EXTINCTION"][:,0]>=.22
plt.scatter(objs[ii]["RA"], objs[ii]["DEC"], color="purple", s=1, marker = "o")
saveloc = "/d/www/kianah/public_html/week2_pos_scatter2.png"
plt.savefig(saveloc)
plt.clf()

print("Figure saved at:" + saveloc)
