import numpy as np
import matplotlib.pyplot as plt

def randomization_test(dat1, dat2, N, stat = 'mean', display = False):

	cdat = np.concatenate((dat1, dat2))
	npoints = cdat.size
	half_ix = np.ceil(npoints / 2.)
	stat_obt = np.abs(dat1.mean() - dat2.mean())
	print 'Observed %s separation : %5.5f' % (stat, stat_obt)
	i = 0
	
	if stat == 'mean':
		func = np.mean
	if stat == 'median':
		func = np.median
	
	stat_i = np.empty(N)
	
	# run N randomization trials
	for n in range(N):
		samp_ix = np.random.permutation(np.arange(npoints))
		samp1 = cdat[samp_ix[:half_ix]]
		samp2 = cdat[samp_ix[half_ix:]]
		stat_i[n] = func(samp1) - func(samp2)
		# if the difference is means is larger than the obtained difference, increment
		if np.abs(stat_i[n]) > stat_obt:
			i += 1
		
	p = np.float32(i) / N
	

	if display:
		fig = plt.figure();
		ax = fig.add_subplot(111)
		ax.hist(stat_i)
		ax.set_title('Difference in %s for N = %i iterations' % (stat, N))
		ax.set_xlabel('|%s_1 - %s_2|' % (stat, stat))
		ax.set_ylabel('Count')
		ax.axvline(stat_obt, color = 'r', ls = '--')

	
	return p