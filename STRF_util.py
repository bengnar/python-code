import numpy as np
import tables
from matplotlib.pyplot import figure, show
import scipy.linalg

def make_delayed(stim, delays, circpad=False):
	"""Creates non-interpolated concatenated delayed versions of [stim] with the given [delays] 
	(in samples).
	
	If [circpad], instead of being padded with zeros, [stim] will be circularly shifted.
	"""
	nt, ndim = stim.shape
	dstims = []
	ndelays = len(delays)
	dstims = np.empty((nt, ndim*ndelays))
	for di, d in enumerate(delays):
		dstim = np.zeros((nt, ndim))
		if d<0: ## negative delay
			dstim[:d,:] = stim[-d:,:]
			if circpad:
				dstim[d:,:] = stim[:-d,:]
		elif d>0:
			dstim[d:,:] = stim[:-d,:]
			if circpad:
				dstim[:d,:] = stim[-d:,:]
		else: ## d==0
			dstim = stim.copy()
		dstims[:, (di*ndim):((di+1)*ndim)] = dstim

	return dstims

def best_corr_vec(wvec, vocab, SU, n=10):
	"""Returns the [n] words from [vocab] most similar to the given [wvec], where each word is represented
	as a row in [SU].  Similarity is computed using correlation."""
	wvec = wvec - np.mean(wvec)
	nwords = len(vocab)
	corrs = np.nan_to_num([np.corrcoef(wvec, SU[wi,:]-np.mean(SU[wi,:]))[1,0] for wi in range(nwords-1)])
	scorrs = np.argsort(corrs)
	words = list(reversed([(corrs[i],vocab[i]) for i in scorrs[-n:]]))
	return words

def get_word_prob():
	"""Returns the probabilities of all the words in the mechanical turk video labels.
	"""
	import constants as c
	import cPickle
	data = cPickle.load(open(c.datafile)) # Read in the words from the labels
	wordcount = dict()
	totalcount = 0
	for label in data:
		for word in label:
			totalcount += 1
			if word in wordcount:
				wordcount[word] += 1
			else:
				wordcount[word] = 1
	
	wordprob = dict([(word, float(wc)/totalcount) for word, wc in wordcount.items()])
	return wordprob

def best_prob_vec(wvec, vocab, space, wordprobs):
	"""Orders the words by correlation with the given [wvec], but also weights the correlations by the prior
	probability of the word appearing in the mechanical turk video labels.
	"""
	words = best_corr_vec(wvec, vocab, space, n=len(vocab)) ## get correlations for all words
	## weight correlations by the prior probability of the word in the labels
	weightwords = []
	for wcorr,word in words:
		if word in wordprobs:
			weightwords.append((wordprobs[word]*wcorr, word))
	
	return sorted(weightwords, key=lambda ww: ww[0])

def find_best_words(vectors, vocab, wordspace, actual, display=True, num=15):
	cwords = []
	for si in range(len(vectors)):
		cw = best_corr_vec(vectors[si], vocab, wordspace, n=num)
		cwords.append(cw)
		if display:
			print "Closest words to scene %d:" % si
			print [b[1] for b in cw]
			print "Actual words:"
			print actual[si]
			print ""
	return cwords

def find_best_stims_for_word(wordvector, decstims, n):
	"""Returns a list of the indexes of the [n] stimuli in [decstims] (should be decoded stimuli)
	that lie closest to the vector [wordvector], which should be taken from the same space as the
	stimuli.
	"""
	scorrs = np.array([np.corrcoef(wordvector, ds)[0,1] for ds in decstims])
	scorrs[np.isnan(scorrs)] = -1
	return np.argsort(scorrs)[-n:][::-1]

def princomp(x, use_dgesvd=False):
	"""Does principal components analysis on [x].
	Returns coefficients, scores and latent variable values.
	Translated from MATLAB princomp function.  Unlike the matlab princomp function, however, the
	rows of the returned value 'coeff' are the principal components, not the columns.
	"""
	
	n,p = x.shape
	#cx = x-np.tile(x.mean(0), (n,1)) ## column-centered x
	cx = x-x.mean(0)
	r = np.min([n-1,p]) ## maximum possible rank of cx

	if use_dgesvd:
		from svd_dgesvd import svd_dgesvd
		U,sigma,coeff = svd_dgesvd(cx, full_matrices=False)
	else:
		U,sigma,coeff = np.linalg.svd(cx, full_matrices=False)
	
	sigma = np.diag(sigma)
	score = np.dot(cx, coeff.T)
	sigma = sigma/np.sqrt(n-1)
	
	latent = sigma**2

	return coeff, score, latent

def eigprincomp(x, npcs=None):
	"""Does principal components analysis on [x].
	Returns coefficients (eigenvectors) and eigenvalues.
	If given, only the [npcs] greatest eigenvectors/values will be returned.
	"""
	n,p = x.shape
	#cx = x-np.tile(x.mean(0), (n,1)) ## column-centered x
	cx = x-x.mean(0)
	r = np.min([n-1,p]) ## maximum possible rank of cx

	xcov = np.cov(cx.T)
	if npcs is not None:
		latent,coeff = scipy.linalg.eigh(xcov, eigvals=(p-npcs,p-1))
	else:
		latent,coeff = np.linalg.eigh(xcov)

	## Transpose coeff, reverse its rows
	return coeff.T[::-1], latent[::-1]

def fixPCs(orig, new):
	"""Finds and fixes sign-flips in PCs by finding the coefficient with the greatest
	magnitude in the [orig] PCs, then negating the [new] PCs if that coefficient has
	a different sign.
	"""
	flipped = []
	for o,n in zip(orig, new):
		maxind = np.abs(o).argmax()
		if o[maxind]*n[maxind]>0:
			## Same sign, no need to flip
			flipped.append(n)
		else:
			## Different sign, flip
			flipped.append(-n)
	
	return np.vstack(flipped)


def plot_model_comparison(corrs1, corrs2, name1, name2, thresh=0.35):
	fig = figure(figsize=(8,8))
	ax = fig.add_subplot(1,1,1)
	
	good1 = corrs1>thresh
	good2 = corrs2>thresh
	better1 = corrs1>corrs2
	#both = np.logical_and(good1, good2)
	neither = np.logical_not(np.logical_or(good1, good2))
	only1 = np.logical_and(good1, better1)
	only2 = np.logical_and(good2, np.logical_not(better1))
	
	ptalpha = 0.3
	ax.plot(corrs1[neither], corrs2[neither], 'ko', alpha=ptalpha)
	#ax.plot(corrs1[both], corrs2[both], 'go', alpha=ptalpha)
	ax.plot(corrs1[only1], corrs2[only1], 'ro', alpha=ptalpha)
	ax.plot(corrs1[only2], corrs2[only2], 'bo', alpha=ptalpha)
	
	lims = [-0.5, 1.0]
	
	ax.plot([thresh, thresh], [lims[0], thresh], 'r-')
	ax.plot([lims[0], thresh], [thresh,thresh], 'b-')
	
	ax.text(lims[0]+0.05, thresh, "$n=%d$"%np.sum(good2), horizontalalignment="left", verticalalignment="bottom")
	ax.text(thresh, lims[0]+0.05, "$n=%d$"%np.sum(good1), horizontalalignment="left", verticalalignment="bottom")
	
	ax.plot(lims, lims, '-', color="gray")
	ax.set_xlim(lims)
	ax.set_ylim(lims)
	ax.set_xlabel(name1)
	ax.set_ylabel(name2)
	
	show()
	return fig

def plot_model_comparison_rois(corrs1, corrs2, name1, name2, roivoxels, roinames, thresh=0.35):
	"""Plots model correlation comparisons per ROI.
	"""
	fig = figure()
	ptalpha = 0.3
	
	for ri in range(len(roinames)):
		ax = fig.add_subplot(4, 4, ri+1)
		ax.plot(corrs1[roivoxels[ri]], corrs2[roivoxels[ri]], 'bo', alpha=ptalpha)
		lims = [-0.3, 1.0]
		ax.plot(lims, lims, '-', color="gray")
		ax.set_xlim(lims)
		ax.set_ylim(lims)
		ax.set_title(roinames[ri])
	
	show()
	return fig

def save_table_file(filename, filedict):
	"""Saves the variables in [filedict] in a hdf5 table file at [filename].
	"""
	hf = tables.openFile(filename, mode="w", title="save_file")
	for vname, var in filedict.items():
		hf.createArray("/", vname, var)
	hf.close()

