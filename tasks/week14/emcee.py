import numpy as np
import random
import emcee
import corner
import matplotlib.pyplot as plt
from IPython.display import display, Math

data = np.loadtxt('/d/scratch/ASTR5160/week13/line.data', unpack=True)

# KFH find variance and means of data
var = [np.var(col, ddof=1) for col in data]
means = [np.mean(col) for col in data]

m, b = 3, 5
n_steps = 10

# KFH Pulling these functions from the tutuorial
def log_prior(theta):
    m, b, log_f = theta
    if 0 < m < 0.5 and 0.0 < b < 10.0 and -10.0 < log_f < 1.0:
        return 0.0
    return -np.inf

def log_probability(theta, x, y, yerr):
    lp = log_prior(theta)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(theta, x, y, yerr)

m_new = [m + np.random.normal() for i in np.arange(n_steps)]
b_new = [b + np.random.normal() for i in np.arange(n_steps)]
full_mb = np.asarray([m_new, b_new]).T
print(full_mb.shape)

# KFH For some reason I'm getting the message that:
# module 'emcee' has no attribute 'EnsembleSampler'...:(
sampler = emcee.EnsembleSampler(n_steps, 2,log_probability, args=(x, y, yerr))

# KFH Run 1000 steps of mcmc
sampler.run_mcmc(full_mb, 1000, progress=True)
samples = sampler.get_chain()

# KFH Find how many steps are needed for the chain to forget starting point
tau = sampler.get_autocorr_time()
print('tau', tau)

# KFH Discard first 100 steps and thin the autocorrelation
flat_samples = sampler.get_chain(discard=100, thin=15, flat=True)
print(flat_samples.shape)

# KFH Create corner plot of resutls
fig = corner.corner(flat_samples, labels=labels, truths=[m_true, b_true])
fig.show()

# KFH Determine uncertainty values at 16th, 50th, and 84th percentile
# KFH of samples in the marginalized distributions
for i in range(2):
	mcmc = np.percentile(flat_samples[:, i], [16, 50, 84])
	q = np.diff(mcmc)
	txt = "\mathrm{{{3}}} = {0:.3f}_{{-{1:.3f}}}^{{{2:.3f}}}"
	txt = txt.format(mcmc[1], q[0], q[1], labels[i])
	print(txt)



