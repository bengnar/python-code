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
nstims = 3

dtype = np.dtype([('gen', 'S2'), ('exp', 'S3'), ('sess', 'S20'), ('unit', 'i4'), ('cf', 'f4'), ('stim_psth', '(%i,3200)f4' % nstims), ('stim', '%ii4' % nstims), ('alphas', '%if4' % nalphas), ('best_alpha', 'f4'), ('wts', '(540,%i)f4' % nstims), ('test_pred', '(320,%i)f4' % nstims), ('train_corr', '(%i,%i)f4' % (nstims, nalphas)), ('test_corr', '(%i,%i)f4' % (nstims, nalphas)), ('test_corr_pval', '(%i,%i)f4' % (nstims, nalphas))])

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

def fit_strfs(stims, resps, figtitle = 'tmp', nalphas = 60, ndelays = 10, ax = None):
	
	# ax, fig = misc.axis_check(ax)
	## Run the STRF fitting
	alphas = np.hstack((0, np.logspace(-40, 10, nalphas-1))) # 10^-2 to 10^5, 45 times
	
	_, nfreqs, nstims = stims.shape

	Wts = np.empty((540, nstims, nalphas))
	Test_pred = np.empty((320, nstims, nalphas))
	Train_corr = np.empty((nstims, nalphas))
	Test_corr = np.empty((nstims, nalphas))
	Test_corr_pval = np.empty((nstims, nalphas))
	
	for i in range(nstims):
		stimnum = i+1
		stim = sparse.csc_matrix(stims[..., i])


		# load the PSTHs
		resp = resps[i, :]
		[train_stim, train_resp], [test_stim, test_resp], [valid_stim, valid_resp] = split_stim_resp(stim, resp, ndelays = ndelays, test_chunk = 4)
	
		# print train_stim.mean(), train_resp.mean()

		# ridge regression on training and validation stim/resp
		wt, train_corr = eigridge(train_stim, valid_stim, train_resp, valid_resp, alphas, saveallwts = True, verbose = False)
		# calculate prediction
		test_pred = test_stim.dot(np.array(wt).squeeze().T)
		# correlation between prediction and actual response
		test_corr = [np.corrcoef(tp, test_resp.squeeze())[0,1] for tp in test_pred.T]
		# significance of pred/actual correlation
		test_corr_pval = [shuffled_correlation_pvalue(tp, test_resp.squeeze(), nboots = 1000) for tp in test_pred.T] #usually 100,000 boots

		Wts[:, i, :] = np.array(wt).squeeze().T
		Test_pred[:, i, :] = test_pred
		Train_corr[i, :] = train_corr
		Test_corr[i, :] = test_corr
		Test_corr_pval[i, :] = test_corr_pval

		
	best_ix = Test_corr.mean(0).argmax()
	best_alpha = alphas[best_ix]
	best_test_pred = Test_pred[..., best_ix]
	best_wts = Wts[..., best_ix]
	
	db_ = np.array(('None', 'None', 'None', -1, -1, resps, np.arange(nstims)+1, alphas, best_alpha, best_wts, best_test_pred, Train_corr, Test_corr, Test_corr_pval), dtype = dtype)
					# ax[-1].plot(test_corr)
	
	return db_
	
	# misc.sameyaxis(ax)
	# fig.suptitle(figtitle)
	# fig.savefig('%s_strf_test_corr.png' % figtitle)
	# plt.close(fig)
		

	# # store the weights (STRFs), correlations
	# wts.append(np.array(wt).squeeze())
	# corrs.append(np.array(corr).ravel())
	# tcorrs.append(tcorr)
	# tcorrpvals.append(np.array(tcorrpval))

	# save_table_file(outfile, dict(wts=np.array(wts),
	# 				  corrs=np.array(corrs),
	# 				  tcorrs=np.array(tcorrs),
	# 				  tcorrpvals=np.array(tcorrpvals),
	# 				  alphas=alphas,
	# 				  delays=delays))
	# print "Saved file %s"%outfile


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
	tmp = np.load(os.path.join(studydir, 'voc.npz'))
	P = tmp['P']
	F = tmp['F']
	T = tmp['T']
	
	return P

def split_stim_resp(stim, resp, nchunks = 10, test_chunk = 0, ndelays = 10):
	# resp = resp[..., np.newaxis]
	# resp = resp.T

	# make concatenated vector of responses
	# resp = resp_.reshape(np.prod(resp_.shape), order = 'C')
	# assert (resp[:resp_.shape[1]]==resp_[0, :]).all()

	resplen = resp.shape[0]
	
	chunks = np.reshape(np.arange(resplen), (nchunks, resplen/nchunks))
	chunk_ix = np.arange(nchunks)
	
	train_chunk = (test_chunk+1)%nchunks
	ix_test = chunk_ix==test_chunk
	ix_valid = chunk_ix==train_chunk
	ix_train = np.vstack((chunk_ix!=test_chunk, chunk_ix!=train_chunk)).all(0)
	
	delays = range(1, 1 + ndelays)
	
	train_stim = sp_make_delayed(stim[chunks[ix_train, :].ravel(), :], delays)
	test_stim = sp_make_delayed(stim[chunks[ix_test, :].ravel(), :], delays)
	valid_stim = sp_make_delayed(stim[chunks[ix_valid, :].ravel(), :], delays)

	
	# train_stim = stim[ix_train, :]
	# test_stim = stim[ix_test, :]
	# valid_stim = stim[ix_valid, :]
	
	train_resp = resp[chunks[ix_train, :].ravel()][..., np.newaxis]
	test_resp = resp[chunks[ix_test, :].ravel()][..., np.newaxis]
	valid_resp = resp[chunks[ix_valid, :].ravel()][..., np.newaxis]

	
	return [train_stim, train_resp], [test_stim, test_resp], [valid_stim, valid_resp]

def plot_strfs(strfs):
	
	nstrfs = len(strfs)
	ndelays, nfreqs = strfs[0].shape
	ncols = np.ceil(np.sqrt(nstrfs))
	fig = plt.figure()
	for i in xrange(nstrfs):
		ax = fig.add_subplot(ncols, ncols, i+1)
		ax.imshow(strfs[i], aspect = 'auto')
	
