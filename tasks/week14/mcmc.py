import numpy as np
import  random

data = np.loadtxt('/d/scratch/ASTR5160/week13/line.data', unpack=True)

# KFH find variance and means of data
var = [np.var(col, ddof=1) for col in data]
means = [np.mean(col) for col in data]


def post_prob(m, b, means, var, m_range, b_range):
	"""
	m_range and b_range should be of the format [m_min, m_max] and [b_min, b_max]
	"""
	
	# KFH Calculate y vals from average x  vals
	x_vals = np.linspace(0.5, 9.5, 10)
	y_vals = (m*np.asarray(x_vals)) + b
	
	if (m > m_range[1] or m < m_range[0]) or (b > b_range[1] or b < b_range[0]):
		ln_postprob = -np.inf
		
	else:
		# KFH Compute the ln of the posterior probability
		ln_postprob = -0.5 * sum((((means - y_vals)**2) / var) \
		 + np.log(2*np.pi*np.asarray(var)))
	
	return(np.exp(ln_postprob), ln_postprob)

def proposal_funct(m, b, means, var, m_range, b_range):
	"""
	Start with m, b, and run post_prob on variations by delta m and b.
	If the new m and b results in a more likely posterior probability
	
	Inputs
	------
	m: float
	- A value corresponding to the slope of mx+b
	b: float
	- A value corresponding to the y-intercept of mx+b
	means: list
	- Mean y values of bins of data
	var: list
	- Variance in y values of bins of data
	m_range: list
	- A list of floats of the form [m_min, m_max], that limit the m priors
	b_range: list
	- A list of floats of the form [b_min, b_max], that limit the b priors
	
	Returns
	-------
	m: float
	- The value of m selected by the metropolis-hastings method
	b: float
	- The value of b selected by the metropolis-hastings method
	ln_init: float
	- The ln of the posterior probability with the initial m and b
	ln_new: float
	- The ln of the posterior probability with the returned values of m and b 
	"""
	# KFH Find initial posterior probability
	init_prob, ln_init = post_prob(m, b, means, var, m_range, b_range)
	
	new_m = m + np.random.normal(0, 0.27)
	new_b = b + np.random.normal(0,0.27)
	
	# KFH Find new posterior probability
	new_prob, ln_new = post_prob(new_m, new_b, means, var, m_range, b_range)
	
	# KFH Determine whether to keep the calculated new_m, new_b, or original vals
	r = new_prob / init_prob
	if r > 1:
		m, b = new_m, new_b
	
	else:
		if random.random() < r:
			m, b = new_m, new_b
		
	# KFH return chain
	return(m, b, ln_init, ln_new)

if __name__=='__main__':

	# KFH Find the inital likelihood and log likelihood with m=3, b=5
	postprob = post_prob(3, 5, means, var, [0,5], [0,8])
	print("likelihood and ln(likelihood)", postprob)

	# KFH Run proposal function 1000 times
	num_trials = 1000
	new = np.asarray([proposal_funct(3, 5, means, var, [0,5], [0,8]) \
			for i in np.arange(num_trials)])
			
	# KFH Calculate how many of these runs resulted in accepted m, b values
	print("Percent runs with new parameters:",
	 1-np.count_nonzero(new.T[0] == 3)/num_trials)
	


	
