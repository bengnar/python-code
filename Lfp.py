import itertools
import RF
import numpy as np
import matplotlib.pylab as plt
from matplotlib.mlab import specgram
from scipy.signal import butter, filtfilt
#from nitime.algorithms import multi_taper_psd




def filter_lfp(lfp, Fs = 384., cutoff = 80., N = 10, btype = 'low'):

	nyq = Fs/2.
	Wn = cutoff/nyq
	b, a, = butter(N, Wn, btype = btype)

	ntrials, npts = lfp.shape
	lfp_filt = np.empty_like(lfp)
	for i in range(ntrials):
		lfp_filt[i, :] = filtfilt(b, a, lfp[i, :])

	return lfp_filt

def dpss(s, Fs):
	F, S, Serr = multi_taper_psd(s, Fs = Fs)
	return F, S, Serr

def mean_dpss(ss, Fs):
	
	S = []
	ntrials, nbins = ss.shape
	for i in range(ntrials):
		F, S_, Serr = dpss(ss[i, :], Fs)
		S.append(S_)

	# S = np.asarray(S)
	return S

def calc_power_spectrum(lfp):
	
	npts, ntrials = lfp.shape
	Fs = 384.
	X = np.zeros()
	x = np.zeros((129, 7, ntrials))
	for i in range(ntrials):
		x[:, :, i] = specgram(lfp[:127, i], Fs = Fs)[0]

def calc_lfp_by_stim(rast, stimparams):

	nstimparams = stimparams.shape[1]
	
	usp = []
	for i in range(nstimparams):	
		usp.append(list(np.unique(stimparams[:, i])))

	nparamlevels = np.empty(nstimparams, dtype = np.int32)
	for i in range(nstimparams):
		nparamlevels[i] = len(usp[i])

	ntrials_per_stim = np.zeros(nparamlevels)

	'''
	compute
	nbins	:	the number of bins
	'''
	dur_ms = rast.shape[1] # number of milliseconds
	t_ms = np.arange(dur_ms) # time indices in ms
	nbins = rast.shape[-1]
	assert np.unique(ntrials_per_stim).size == 1

	ntrials = np.int32(ntrials_per_stim[0])
	
	psth_shape = np.hstack((nparamlevels, nbins))
	psth = np.zeros(psth_shape)
	
	combinations = []
	combinations_ix = []
	for i in itertools.product(*usp):
		combinations.append(i)
		combinations_ix_ = []
		for j, i_ in enumerate(i):
			q = (np.array(usp[j])==i_).nonzero()[0][0]
			combinations_ix_.append(q)
		combinations_ix.append(combinations_ix_)
		
	for m, n in zip(combinations, combinations_ix):
		ix = RF.get_trials(stimparams, m)
		ntrials = ix.size
		lfp_ = rast[ix, :]
		psth[tuple(n)] = lfp_.sum(0)
				
	return psth, usp