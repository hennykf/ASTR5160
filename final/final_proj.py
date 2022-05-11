import numpy as np
import matplotlib.pyplot as plt
import emcee
import corner
from astropy.table import  Table
from IPython.display import display, Math
from matplotlib.gridspec import GridSpec

def linear_eq(m, b, x):
	"""
	Finds y values for the equation y=mx+b based on x values, and
	m and b parameters
	"""
	y = m*x + b
	return(y)
	
def quad_eq(a2, a1, a0, x):
	"""
	Finds y values for the equation y = a2*x**2 + a1*x + a0, based
	on given x values, and a2, a1, and a0 parameters
	"""
	y = a2*(x**2) + a1*x + a0
	return(y)

def log_prior_lin(m, b, m_min_max, b_min_max):
	"""
	Sets the log prior value for a linear fit of the form
	y=mx+b
	
	Inputs
	------
	m: float
	- The value m in the fit y=mx+b
	b: float
	- The value b in the fit y=mx+b
	m_min_max: list
	- list of the minimum and maximum permitted values for m; [mmin, mmax]
	b_min_max: list
	- list of minimum and maximum permitted values for b; [bmin, bmax]
	
	Returns
	-------
	Value of 0 or negative infinity depending on if value of m and b are within priors
	"""
	m_min, m_max = m_min_max
	b_min, b_max = b_min_max
	
	# KFH Return value based on fit info
	if m_min <= m <= m_max and b_min <= b <= b_max:
		return 0.0
	
	return -np.inf
	
def log_prior_quad(a2, a1, a0, a2_min_max, a1_min_max, a0_min_max):
	"""
	Sets the log prior value for a linear fit of the form
	y=mx+b
	
	Inputs
	------
	a2: float
	- The value a2 in the fit y=a2*x^2 + a1*x + a0
	a1: float
	- The value a1 in the fit y=a2*x^2 + a1*x + a0
	a0: float
	- The value a0 in the fit y=a2*x^2 + a1*x + a0
	a2_min_max: list
	- list of the minimum and maximum permitted values for a2; [a2min, a2max]
	a1_min_max: list
	- list of the minimum and maximum permitted values for a1; [a1min, a1max]
	a0_min_max: list
	- list of the minimum and maximum permitted values for a0; [a0min, a0max]
	
	Returns
	-------
	Value of 0 or negative infinity depending on whether a0, a1, and a2 are within
	priors
	"""
	a2_min, a2_max = a2_min_max
	a1_min, a1_max = a1_min_max
	a0_min, a0_max = a0_min_max
	
	# KFH Return value based on fit info
	if (a2_min <= a2 <= a2_max and a1_min <= a1 <= a1_max) and a0_min <= a0 <= a0_max:
		return 0.0
	
	return -np.inf
	
def log_posterior_lin(mb, xvals, yvals, var, m_min_max, b_min_max):
	"""
	Calculates the log of the posterior probability for a linear fit
	
	Inputs
	------
	mb: list
	- Initial m and b values in the equation y=mx+b, ex. [m, b]
	xvals: array or list
	- generated x values over  span of data
	yvals: array or list
	- Actual y data
	var: array or list
	- variance, or the square of yerr
	m_min_max: list
	- list of minimum and maximum values for m, ex. [0,5]
	b_min_max: list
	- list  of minimum and maximum values for b, ex. [0,5]
	"""
	m, b = mb
	
	# KFH Find log priors
	ln_prior = log_prior_lin(m, b, m_min_max, b_min_max)
	
	# KFH calculate y values based on model
	y_predict = linear_eq(m, b, xvals)
	
	# KFH Calculate ln of likelihood and posterior prob
	ln_likelihood = (-0.5)*sum(((yvals - y_predict)**2)/var + np.log(2*np.pi*var))
	ln_posterior = ln_prior + ln_likelihood
	
	return(ln_posterior)
	
def log_posterior_quad(a2a1a0, xvals, yvals, var, a2_min_max, a1_min_max, a0_min_max):
	"""
	Calculates the log of the posterior probability for a quadratic fit
	
	Inputs
	------
	a2a1a0: list
	- Initial a2, a1, a0 values in the equation y=a2*x^2 + a1*x + a0
	xvals: array or list
	- generated x values over  span of data
	yvals: array or list
	- Actual y data
	var: array or list
	- variance, or the square of yerr
	a2_min_max: list
	- list of minimum and maximum values for a2, ex. [0,5]
	a1_min_max: list
	- list  of minimum and maximum values for a1, ex. [0,5]
	a0_min_max: list
	- list  of minimum and maximum values for a0, ex. [0,5]
	"""
	a2, a1, a0 = a2a1a0
	
	# KFH Find log priors
	ln_prior = log_prior_quad(a2, a1, a0, a2_min_max, a1_min_max, a0_min_max)
	
	# KFH calculate y values based on model
	y_predict = quad_eq(a2, a1, a0, xvals)
	
	# KFH Calculate ln of likelihood and posterior prob
	ln_likelihood = (-0.5)*sum(((yvals - y_predict)**2)/var + np.log(2*np.pi*var))
	ln_posterior = ln_prior + ln_likelihood
	
	return(ln_posterior)
	
def make_plots(data, labels, samples, flat_samples, coeffs, lin_quad, nd):
	"""
	Generates two plot windows, one with the random walk plots for each variable and
	the scatter plot of data with the 50th percentile best fit overlaid, and a second
	plot which is a corner plot showing the distribution of likelihood of variable values
	
	Inputs
	------
	data: astropy Table
	- An astropy table containing the columns x, y, and yerr
	labels: list
	- List of the names of the coefficients, for example ['m', 'b']
	samples: ndarray
	- Array of samples from doing sampler.get_chain()
	flat_samples: ndarray
	- Flattened array of samples from sampler.get_chain()
	coeffs: list
	- 50th percentile values for coefficients obtained from mcmc
	lin_quad: str
	- "lin" shows a linear fit, using coefficients from doing an mcmc linear simulation
	- "quad" shows a quadratic fit, using coefficients from doing an mcmc quadratic simulation
	nd: int
	- The ndim used in the mcmc simulation
	
	Outputs
	-------
	Two plots, which can be saved or shown outside the function
	"""
	
	# KFH Initialize grid layout
	fig, axs = plt.subplots(2, 2)
	gs = axs[0,0].get_subplotspec().get_gridspec()
	
	for a in axs[:,:].flat:
		a.remove()

	# KFH Plot scatter plot with best fit line
	scatter_subfig = fig.add_subfigure(gs[1,:])
	axs_scat = scatter_subfig.subplots(1)
	axs_scat.errorbar(data['x'], data['y'], yerr = data['yerr'],
							 fmt='o', label='data')
							 
	# KFH Plot depending on whether linear or quadratic
	if lin_quad == 'lin':
		axs_scat.plot(x, linear_eq(coeffs[0], coeffs[1], x),
					 label='fit: y = {0:.3f}x + {1:.3f}'.format(coeffs[0], coeffs[1]))
	elif lin_quad == 'quad':
		axs_scat.plot(x, quad_eq(coeffs[0], coeffs[1], coeffs[2], x),
					 label='fit: y = {0:.3f}x^2 + {1:.3f}x + {2:.3f}'.format(coeffs[0],
					 											 coeffs[1], coeffs[2]))
	axs_scat.legend()
	
	# KFH get chain results
	subfig = fig.add_subfigure(gs[0,:])
	ax_sub = subfig.subplots(nd, sharex=True)
	
	# KFH Make a distribution plot for each parameter
	for i in range(nd):
		ax = ax_sub[i]
		ax.plot(samples[:,:,i], "k", alpha=0.3)
		ax.set_xlim(0, len(samples))
		ax.set_ylabel(labels[i])
		ax.yaxis.set_label_coords(-0.1, 0.5)
	
	ax_sub[-1].set_xlabel('step number');
	
	# KFH Make a corner plot
	r = corner.corner(flat_samples, labels=labels);
	#plt.show()


if __name__=='__main__':

	# KFH Read in data and find variance
	data = Table.read('/d/scratch/ASTR5160/final/dataxy.fits')
	var = np.asarray(data['yerr'])**2
	x = np.linspace(1, 19, 10)

	# KFH Run mcmc for linear fit, starting with eyeballed values from raw data
	m_start = -.875
	b_start = 3
	
	pos = np.array([m_start, b_start]) + 1e-4 * np.random.randn(32, 2)
	nwalkers, ndim = pos.shape
	
	sampler = emcee.EnsembleSampler(nwalkers, ndim, log_posterior_lin,
									 args = (x, np.asarray(data['y']),
									 		 var, [-2, 0], [-3, 8]))
	
	sampler.run_mcmc(pos, 5000, progress=True)
	samples = sampler.get_chain()
	labels = ['m', 'b']
	flat_samples = sampler.get_chain(discard=100, thin=15, flat=True)
	
	# Get final values and uncertainty interval and print
	final_vals = []
	for i in range(ndim):
		mcmc = np.percentile(flat_samples[:,i], [16,50,84])
		q = np.diff(mcmc)
		txt = "{0} parameter, 50th={1:.3f}, +{2:.3f}, -{3:.3f}".format(labels[i],
																 mcmc[1], q[0], q[1])
		final_vals.append(mcmc[1])
		print(txt)
	
	# Run mcmc for quadratic fit, starting with eyeballed values from raw data
	a2_start = .05
	a1_start  = -2
	a0_start = 0
	
	pos_quad = np.array([a2_start, a1_start, a0_start]) + 1e-4 * np.random.randn(32, 3)
	nwalkers_quad, ndim_quad = pos_quad.shape
	
	sampler_quad = emcee.EnsembleSampler(nwalkers_quad, ndim_quad, log_posterior_quad,
									 	args = (x, np.asarray(data['y']), var,
									 	 [-10, 10], [-10, 10], [-10,20]))
	sampler_quad.run_mcmc(pos_quad, 5000, progress=True)
	samples_quad = sampler_quad.get_chain()
	labels_quad = ['a2', 'a1', 'a0']
	flat_samples_quad = sampler_quad.get_chain(discard=100, thin=15, flat=True)
	
	# KFH Get final values and uncertainty interval and print
	final_vals_quad = []
	for i in range(ndim_quad):
		mcmc = np.percentile(flat_samples_quad[:,i], [16,50,84])
		q = np.diff(mcmc)
		txt = "{0} parameter, 50th={1:.3f}, +{2:.3f}, -{3:.3f}".format(labels_quad[i],
																 mcmc[1], q[0], q[1])
		final_vals_quad.append(mcmc[1])
		print(txt)
		
	print("If we look at the corner plot for the quadratic fit, we  can see that the a2\
 distribution is ofset from 0 by quite a bit. Only the lefthand tail of the\
 distribution overlaps with 0 at all. This indicates that most  fits for a2\
 are greater than zero, which shows us that a quadratic fit is most appropriate\
 for this dataset. If the  distribution of a2 values was centered at 0, we could\
 say with some confidence that a linear fit would be more appropriate, but that is\
 clearly not the case with this data.")

	# KFH Plot linear and quadratic  plots
	make_plots(data, labels, samples, flat_samples, final_vals, lin_quad='lin', nd=ndim)	
	make_plots(data, labels_quad, samples_quad, flat_samples_quad,
				 final_vals_quad, lin_quad='quad', nd=ndim_quad)
	plt.show()
	

