
'''Sparse STRF code'''
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import cm

import scipy.io as sio
from scipy import sparse
import numpy as np
import os, shutil, glob, re
import h5py
import Spikes
import misc

# from text.movie.util.util import save_table_file
from STRF_utils import make_delayed, sp_make_delayed, counter
from ridge import eigridge

from STRF_model_significance import correlation_pvalue, exact_correlation_pvalue, rolled_correlation_pvalue, shuffled_correlation_pvalue

studydir = '/Volumes/BOB_SAGET/Fmr1_voc'
nalphas = 60

dtype = np.dtype([('gen', 'S2'), ('exp', 'S3'), ('sess', 'S20'), ('unit', 'i4'), ('cf', 'f4'), ('alphas', '%if4' % nalphas), ('best_alpha', 'f4'), ('best_alpha_ix', 'i4'), ('wts', '190f4'), ('test_pred', '960f4'), ('train_corr', '%if4' % nalphas), ('test_corr', '%if4' % nalphas), ('test_corr_pval', '%if4' % nalphas)])

p = re.compile('(\d+)')

def fit_strfs_all(experiments):
	
	db = np.empty(0, dtype = dtype)
	
	for experiment in experiments:

		_, gen, exp, sess = experiment.split('_')
		stims = load_stimulus()
		if len(stims.shape)==2:
			stims = stims[..., np.newaxis]
		nstims = stims.shape[2]
	
		fpaths = glob.glob(os.path.join(studydir, 'Sessions', experiment, 'fileconversion', 'VOC*.h5'))
		for i, fpath in enumerate(fpaths):
			absol, relat = os.path.split(fpath)
			fname, _ = os.path.splitext(relat)
			unitnum = np.int32(p.findall(fname)[0])
			figtitle = os.path.splitext(fpath)[0]
			resps = load_responses(fpath)
			resps = (resps - resps.mean()) / resps.std() # zscore responses
			db_ = fit_strfs(stims, resps, figtitle = figtitle, nalphas = nalphas, ndelays = 10)
			db_['gen'] = gen
			db_['exp'] = exp
			db_['sess'] = sess
			db_['unit'] = unitnum
			db.resize(db.size+1)
			db[-1] = db_
	
	return db

def fit_strfs(stims, resps, figtitle = 'tmp', nalphas = 500, ndelays = 10, ax = None):
	
	## Run the STRF fitting
	alphas = np.hstack((0, np.logspace(-40, 10, nalphas-1))) # 10^-2 to 10^5, 45 times
	
	_, nfreqs, nstims = stims.shape

	# load the PSTHs
	[train_stim, train_resp], [test_stim, test_resp], [valid_stim, valid_resp] = split_stim_resp(stims, resps, ndelays = ndelays, test_chunk = 4)

	# ridge regression on training and validation stim/resp
	wt, train_corr = eigridge(train_stim, valid_stim, train_resp, valid_resp, alphas, saveallwts = True, verbose = False)
	train_corr
	wt = np.asarray(wt).squeeze()
	train_corr = np.asarray(train_corr).squeeze()
	

	# calculate prediction
	test_pred = test_stim.dot(np.array(wt).squeeze().T)
	# correlation between prediction and actual response
	test_corr = [np.corrcoef(tp, test_resp.squeeze())[0,1] for tp in test_pred.T]
	# significance of pred/actual correlation
	test_corr_pval = [shuffled_correlation_pvalue(tp, test_resp.squeeze(), nboots = 100000) for tp in test_pred.T] #usually 100,000 boots
	test_corr_pval = -np.ones(nalphas)
	
	test_corr = np.asarray(test_corr)

	best_alpha_ix = test_corr.argmax()
	best_wt = wt[best_alpha_ix, :]
	best_alpha = alphas[best_alpha_ix]
	best_test_pred = test_pred[..., best_alpha_ix]
	
	db_ = np.array(('None', 'None', 'None', -1, -1, alphas, best_alpha, best_alpha_ix, best_wt, best_test_pred, train_corr, test_corr, test_corr_pval), dtype = dtype)
					# ax[-1].plot(test_corr)
	
	return db_
	


def load_responses(fpath):
	
	'''RESPONSE'''
	f = h5py.File(fpath, 'r')
	rast = f['rast'].value
	stimparams = f['stimID'].value
	stimparams = stimparams[:, 0]
	stimparams = stimparams[:, np.newaxis]
	f.close()
	
	bins_strf = np.arange(0, 16.005, 0.005)
	
	(resps, _) = Spikes.calc_psth_by_stim(rast, stimparams, bins = bins_strf)
	
	return resps

	
def load_stimulus():
	
	'''load already-processed stimulus spectrogram'''
	tmp = np.load(os.path.join(studydir, 'stims', 'voc_43to80kHz.npz'))
	P = tmp['P']
	F = tmp['F']
	T = tmp['T']
	
	return P

def split_stim_resp(stims, resps, nchunks = 10, test_chunk = 0, ndelays = 10):
	# resp = resp[..., np.newaxis]
	# resp = resp.T

	# make concatenated vector of responses
	# resp = resp_.reshape(np.prod(resp_.shape), order = 'C')
	# assert (resp[:resp_.shape[1]]==resp_[0, :]).all()

	resplen, nfreqs, nstims = stims.shape
	
	chunks = np.reshape(np.arange(resplen), (nchunks, resplen/nchunks))
	chunk_ix = np.arange(nchunks)
	
	train_chunk = (test_chunk+1)%nchunks
	chunks_test = chunk_ix==test_chunk
	chunks_valid = chunk_ix==train_chunk
	chunks_train = np.vstack((chunk_ix!=test_chunk, chunk_ix!=train_chunk)).all(0)
	
	ix_test = chunks[chunks_test, :].ravel()
	ix_train = chunks[chunks_train, :].ravel()
	ix_valid = chunks[chunks_valid, :].ravel()

	testlen = ix_test.size;	trainlen = ix_train.size; validlen = ix_valid.size

	# initialize empty stimulus groups
	train_stim = np.empty((trainlen*nstims, nfreqs))
	test_stim = np.empty((testlen*nstims, nfreqs))
	valid_stim = np.empty((validlen*nstims, nfreqs))
	
	# loop through stimuli and concatenate them together in their proper groups
	for i in xrange(nstims):
		stim = stims[..., i]
		train_stim[(i*trainlen):((i+1)*trainlen)] = stim[ix_train, :]
		test_stim[(i*testlen):((i+1)*testlen)] = stim[ix_test, :]
		valid_stim[(i*validlen):((i+1)*validlen)] = stim[ix_valid, :]
		
	delays = range(1, 1 + ndelays)
	
	train_stim = sparse.csc_matrix(train_stim)
	test_stim = sparse.csc_matrix(test_stim)
	valid_stim = sparse.csc_matrix(valid_stim)
	
	train_stim = sp_make_delayed(train_stim, delays)
	test_stim = sp_make_delayed(test_stim, delays)
	valid_stim = sp_make_delayed(valid_stim, delays)

	train_resp = resps[:, ix_train].ravel(order = 'C')[..., np.newaxis]
	test_resp = resps[:, ix_test].ravel(order = 'C')[..., np.newaxis]
	valid_resp = resps[:, ix_valid].ravel(order = 'C')[..., np.newaxis]

	
	return [train_stim, train_resp], [test_stim, test_resp], [valid_stim, valid_resp]

def plot_strfs(strfs):
	
	nstrfs = len(strfs)
	ndelays, nfreqs = strfs[0].shape
	ncols = np.ceil(np.sqrt(nstrfs))
	fig = plt.figure()
	for i in xrange(nstrfs):
		ax = fig.add_subplot(ncols, ncols, i+1)
		ax.imshow(strfs[i], aspect = 'auto')
	
