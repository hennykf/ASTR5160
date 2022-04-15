from astropy.table import Table
import  numpy as np
import matplotlib.pyplot as plt
from scipy.stats import chi2
import pandas as pd


def y_predict(mb, xvals):
	"""
	Inputs a list of [b, m] values and x values and generates the
	y values based on the formula y=b+mx
	"""
	y_vals = mb[1]*xvals + mb[0]
	return(y_vals)

# KFH read in the data and find means and variance
data = Table.read('/d/scratch/ASTR5160/week13/line.data', format="ascii")
means = [np.mean(data[col]) for col in data.columns]
var = [np.var(data[col], ddof=1) for col in data.columns]
print("Means:",means)
print("Variance:", var)

# KFH generate test x vals to plot
test_x = np.linspace(0,10,11)

# KFH Plot data
plt.scatter(np.linspace(0.5, 9.5, 10), means, label='means')
plt.plot(test_x, 5+3*test_x, label='5+3x') 
plt.plot(test_x, 4.5+3*test_x, label='4.5+3x')
plt.plot(test_x, 4+3*test_x, label='4+3x')
plt.plot(test_x, 5.5+3*test_x, label='5.5+3x')
plt.plot(test_x, 4.5+3.1*test_x, label='4.5+3.1x')
plt.ylim([0,35])
plt.legend()
plt.show(block=False)

# KFH Create some test m and b values based on plotting
mb = np.asarray([[5,3], [4.5,3],[4, 3], [5.5, 3], [4.5, 3.1]])

# KFH find chi2 for each parameter pair
yvals_predicted = [y_predict(i, np.linspace(0.5, 9.5, 10)) for i in mb]
chi2 = [sum(((list(means) - i)**2) / var) for i in yvals_predicted]

# KFH Print out m and b parameters and chi2 for that fit
for i in np.arange(len(chi2)):
	print(mb[i], chi2[i])

# KFH As indicated by the low chi2 value, 5+3x is the best fit.
		
	
