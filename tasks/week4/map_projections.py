import matplotlib.pyplot as plt
import numpy as np
from numpy.random import random

ra = 2*np.pi*(random(10000)-0.5)
dec = np.arcsin(1.-random(10000)*2)

fig_save = "/d/www/kianah/public_html/week4_map_projections.png"
plt.clf()

fig = plt.figure(figsize=[10,10])

#KFH Plot x-y
# Note: There are more points near the equator of the sphere when
# projected into the cartesian grid.
ax1 = fig.add_subplot(221)
ax1.scatter(ra, dec, marker='o', color='m', s=0.7, alpha=0.5)
ax1.grid(color='b', linestyle='solid', linewidth=0.5)
ax1.set_xlabel("RA (radians)")
ax1.set_ylabel("Dec (radians)")
ax1.set_title("Cartesian projection")

#KFH Plot aitoff projection
ax = fig.add_subplot(222, projection='aitoff')
ax.scatter(ra, dec, marker='o', color='m', s=0.7, alpha=0.5)
xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
ax.set_xticklabels(xlab, weight=800)
ax.grid(color='b', linestyle='solid', linewidth=0.5)
ax.set_title("Aitoff Projection")

#KFH Plot Lambert projection
ax2 = fig.add_subplot(223, projection='lambert')
ax2.scatter(ra, dec, marker='o', color='m', s=0.7, alpha=0.5)
xlab = ['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h']
ax2.set_xticklabels(xlab, weight=800)
ax2.grid(color='b', linestyle='solid', linewidth=0.5)
ax2.set_title("Lambert Projection")

fig.savefig(fig_save)

print("figure saved at " + fig_save) 
