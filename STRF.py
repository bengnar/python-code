import h5py, os, glob, re
from matplotlib.mlab import specgram
import Spikes
from STRF_util import make_delayed
from ridge import ridge, bootstrap_ridge, find_best_alpha, set_weights, avg_corr_mats

basedir = '/Volumes/BOB_SAGET/Fmr1_voc/'

fs_stim = 250000 # stimulus sample rate
fs_strf = 200. # corrresponds to 5-ms time bins
per_strf = 1./fs_strf
stim_dur = 15.85 # duration (seconds) of the sound stimulus
resp_dur = 20. # duration (seconds) of the response
bins_strf = np.arange(0, resp_dur+per_strf, per_strf)
nbins_strf = bins_strf.size


def process_stim():
	'''STIMULUS'''
	pos = np.loadtxt(os.path.join(basedir, 'stims', 'pos.txt'))
	npts = np.loadtxt(os.path.join(basedir, 'stims', 'npts.txt'))

	# load stimulus .wav file
	stim_fname = os.path.join(basedir, 'stims', 'voc.txt')
	s = np.loadtxt(stim_fname)

	increment = fs_stim / fs_strf
	NFFT = 128
	noverlap = NFFT-increment
	bins = np.arange(0, 20, 0.005)
	nsamp_f = bins.size
	stims = np.empty((0, (NFFT/2)+1))
	for pos_, npts_ in zip(pos, npts):
		# print pos_-1, pos_-1 + npts_
		s_ = s[(pos_-1):((pos_-1) + npts_)]
		nsamp_t = s_.size
		(P_, F, T) = specgram(s_, Fs = fs_stim, NFFT = NFFT, noverlap = noverlap)
		P_ = np.hstack((P_, np.zeros((P_.shape[0], nbins_strf-P_.shape[1]))))
		stims = np.vstack((stims, P_.T))
	
	f_hp = 20000
	stims = stims[:, F>f_hp]
	F = F[F>f_hp]
	del(s)


'''load already-processed stimulus spectrogram'''
tmp = np.load(os.path.join(basedir, 'voc.npz'))
P = tmp['P']
F = tmp['F']
T = tmp['T']

'''create delayed stimulus matrix'''
delaylen = 40 #40 # number of delays
delaytime = 5 # in msec
delays = range(delaylen) # 0 - 200 msec

dStims = make_delayed(stims, delays)

# fig = plt.figure()
# ax = []
# for i in range(3):
# 	ax.append(fig.add_subplot(3, 1, i+1))
# 	voc.plot_spec(P[i]**0.3, T = T[i], F = F[i], ax = ax[-1])
	
'''RESPONSE'''
f = h5py.File(os.path.join(basedir, 'voc_ko_nai_20130116', 'fileconversion', 'VOC002.h5'))
rast = f['rast'].value
stimparams = f['stimID'].value
stimparams = stimparams[:, 0]
stimparams = stimparams[:, np.newaxis]
f.close()
(resps, _) = Spikes.calc_psth_by_stim(rast, stimparams, bins = bins_strf)
resps = resps[..., np.newaxis]
# resps = resps.T

# make concatenated vector of responses
resps = resps.reshape(np.prod(resps.shape), order = 'C')

nchunks = 5
resplen = resps.shape[0]
assert (resplen%nchunks) == 0
chunklen = resplen / nchunks
 # make sure nothing is getting left out
chunks = np.reshape(np.arange(resplen), (nchunks, resplen/nchunks))



alphas = np.logspace(0, 5, 20) # 10^0 to 10^5, 45 times
nBoots = 25

allinds = np.zeros((resps.shape[0],)).astype(bool) ## This is an array of boolean Falses
allwts, allcorrs, allstrfs = list(), list(), list()

for c in range(nchunks):
    print "Using chunk %d ..."%c
    Pinds = allinds.copy() ## Create a copy (so changing it doesn't change allinds, if you want to do this many times)
    Pinds[chunks[c, :]] = True ## Only the indices of the test data are True
    Rinds = ~Pinds ## Only the indices of the training data are True
    
    Rresps = resps[Rinds, :][..., np.newaxis] ## Selects rows of resps where Rinds is True (training)
    Presps = resps[Pinds, :][..., np.newaxis] # test set
    
    ridge_corrs = dict()
    
    Rstim = dStims[Rinds, :] # training stimulus
    Pstim = dStims[Pinds, :] # validation stimulus
    
    ridge_corrs = bootstrap_ridge(Rstim, Rresps, Pstim, Presps, alphas, nBoots, 100, nchunks, corrmin=0.2, joined=None)
    print "Finding average correlations for VOC and VBK..."
    
    valphas = find_best_alpha(Rresps, ridge_corrs, alphas, nBoots, joined=None)
    print valphas
    
    wt, corrs, strfs = dict(), dict(), dict()
    
    Rstim = dStims[Rinds] # training stimulus
    Pstim = dStims[Pinds] # validation stimulus
    
    wt, corrs = set_weights(Rstim, Rresps, Pstim, Presps, alphas, valphas)
    strfs = [w.reshape((len(delays), stims.shape[1])).T for w in wt.T]
    
    allwts.append(wt)
    allcorrs.append(corrs)
    allstrfs.append(strfs)

print "Finding the mean STRF and mean STRF performance..."
meanstrfs, meancorrs = dict(), dict()

zscore = lambda m: (m-m.mean(0))/m.std(0) # Zscores rows

meanstrfs = np.mean([a for a in allstrfs], 0)
meancorrs = np.mean([a for a in allcorrs], 0)

allvalidstrfs = np.zeros((len(pens),stims.shape[1],delaylen,len(chunks_sent)))
allvalidcorrs = np.zeros((len(pens),len(chunks_sent)))





