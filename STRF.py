import h5py, os, glob, re
from matplotlib.mlab import specgram
import numpy as np
import Spikes
from STRF_util import make_delayed
from ridge import ridge, bootstrap_ridge, find_best_alpha, set_weights, avg_corr_mats

studydir = '/Volumes/BOB_SAGET/Fmr1_voc/'

fs_stim = 250000 # stimulus sample rate
fs_strf = 200. # corrresponds to 5-ms time bins
per_strf = 1./fs_strf
stim_dur = 15.85 # duration (seconds) of the sound stimulus
resp_dur = 16. # duration (seconds) of the response
bins_strf = np.arange(0, resp_dur+per_strf, per_strf)
nbins_strf = bins_strf.size-1


# fig = plt.figure()
# ax = []
# for i in range(3):
# 	ax.append(fig.add_subplot(3, 1, i+1))
# 	voc.plot_spec(P[i]**0.3, T = T[i], F = F[i], ax = ax[-1])

def calc_strfs(experiment, pennum, prefix = 'VOC', nchunks = 5):

	resps = load_response(experiment, pennum, prefix = prefix)
	dStims = load_stimulus()
	
	resplen = resps.shape[0]
	assert (resplen%nchunks) == 0
	chunklen = resplen / nchunks # make sure nothing is getting left out

	# make chunk indices
	chunks = np.reshape(np.arange(resplen), (nchunks, resplen/nchunks))

	alphas = np.logspace(0, 5, 20) # 10^0 to 10^5, 45 times
	nboots = 25

	allinds = np.zeros((resps.shape[0],)).astype(bool) ## This is an array of boolean Falses
	allwts, allcorrs, allstrfs = list(), list(), list()

	for c in range(nchunks):
	    print "Using chunk %d ..."%c
	    Pinds = allinds.copy() ## Create a copy (so changing it doesn't change allinds, if you want to do this many times)
	    Pinds[chunks[c, :]] = True ## Indices of the test data are True
	    Rinds = ~Pinds ## Only the indices of the training data are True
	    Rresp = resps[Rinds, :][..., np.newaxis] ## Selects rows of resps where Rinds is True (training)
	    Presps = resps[Pinds, :][..., np.newaxis] # test set
    
	    ridge_corrs = dict()
    
	    Rstim = dStims[Rinds, :] # training stimulus
	    Pstim = dStims[Pinds, :] # validation stimulus
    
	    ridge_corrs = bootstrap_ridge(Rstim, Rresp, Pstim, Presps, alphas, nboots, 100, nchunks, corrmin = 0.2, joined=None)
	    print "Finding average correlations for VOC and VBK..."
    
	    valphas = find_best_alpha(Rresp, ridge_corrs, alphas, nboots, joined=None)
	    print valphas
    
	    wt, corrs, strfs = dict(), dict(), dict()
    
	    Rstim = dStims[Rinds] # training stimulus
	    Pstim = dStims[Pinds] # validation stimulus
    
	    wt, corrs = set_weights(Rstim, Rresp, Pstim, Presps, alphas, valphas)
	    strfs = [w.reshape((len(delays), stims.shape[1])).T for w in wt.T]
    
	    allwts.append(wt)
	    allcorrs.append(corrs)
	    allstrfs.append(strfs)

	print "Finding the mean STRF and mean STRF performance..."
	meanstrfs, meancorrs = dict(), dict()

	zscore = lambda m: (m-m.mean(0))/m.std(0) # Zscores rows

	meanstrfs = np.mean([a for a in allstrfs], 0)
	meancorrs = np.mean([a for a in allcorrs], 0)

	# allvalidstrfs = np.zeros((len(pens),stims.shape[1],ndelays,len(chunks_sent)))
	# 
	# allvalidcorrs = np.zeros((len(pens),len(chunks_sent)))

	return allwts, allcorrs, allstrfs


def load_response(experiment, pennum, prefix = 'VOC'):
	
	'''RESPONSE'''
	f = h5py.File(os.path.join(studydir, 'voc_ko_nai_20130116', 'fileconversion', '%s%3.3u.h5' % (prefix, pennum)), 'r')
	rast = f['rast'].value
	stimparams = f['stimID'].value
	stimparams = stimparams[:, 0]
	stimparams = stimparams[:, np.newaxis]
	f.close()
	
	resps_, _ = Spikes.calc_psth_by_stim(rast, stimparams, bins = bins_strf)
	# resps = resps[..., np.newaxis]
	# resps = resps.T

	# make concatenated vector of responses
	# resps = resps_.reshape(np.prod(resps_.shape), order = 'C')
	# assert (resps[:resps_.shape[1]]==resps_[0, :]).all()
	
	return resps
	
def load_stimulus(ndelays = 40, delaytime = 5):
	'''
	P should have been saved as a (no. time points x no. frequencies x no.stimuli)
	Output:
		dStims
	'''
	
	'''load already-processed stimulus spectrogram'''
	tmp = np.load(os.path.join(studydir, 'voc.npz'))
	P = tmp['P']
	F = tmp['F']
	T = tmp['T']

	if len(P.shape)==2:
		P = P[..., np.newaxis]
	
	nbins, nfreqs, nstims = P.shape

	'''create delayed stimulus matrix'''
	dStims = np.empty((nbins, nfreqs*ndelays, nstims))
	delays = range(ndelays) # 0 - 200 msec
	for i in range(nstims):
		dStims[:, :, i] = make_delayed(P[:, :, i], delays)
	return dStims
	
def process_stim(fs_stim = 250000, fs_strf = 200):
	'''STIMULUS'''
	pos = np.loadtxt(os.path.join(studydir, 'stims', 'pos.txt'))
	npts = np.loadtxt(os.path.join(studydir, 'stims', 'npts.txt'))
	nstims = pos.size

	# load stimulus .wav file
	stim_fname = os.path.join(studydir, 'stims', 'voc.txt')
	s = np.loadtxt(stim_fname)

	increment = fs_stim / fs_strf
	NFFT = 128
	noverlap = NFFT-increment
	stims = np.empty((nbins_strf, (NFFT/2)+1, nstims))
	for i, (pos_, npts_) in enumerate(zip(pos, npts)):
		# print pos_-1, pos_-1 + npts_
		s_ = s[(pos_-1):((pos_-1) + npts_)]
		nsamp_t = s_.size
		(P_, F, T) = specgram(s_, Fs = fs_stim, NFFT = NFFT, noverlap = noverlap)
		P_ = np.hstack((P_, np.zeros((P_.shape[0], nbins_strf-P_.shape[1]))))
		stims[:, :, i] = P_.T
	
	f_hp = 20000
	stims = stims[:, F>f_hp, :]
	F = F[F>f_hp]
	
	np.savez(os.path.join(studydir, 'voc.npz'), P = stims, F = F, T = T)
	del(s)	
	
	