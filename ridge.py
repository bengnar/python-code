#import scipy
import numpy as np
M = np.matrix

def ridge(Rstim, Pstim, Rresp, Presp, alphas, rPresp=None, rPstim=None, saveallwts=True,
		  stop_early=False, dtype=np.single, corrmin=0.2):
	"""Ridge regresses [Rstim] onto [Rresp] for each ridge parameter in [alpha].  Returns the fit
	linear weights for each alpha, as well as the distributions of correlations on a held-out test
	set ([Pstim] and [Presp]).  Note that these should NOT be the "real" held-out test set, only a
	small test set used to find the optimal ridge parameter.
	
	If an [rPresp] and [rPstim], or 'real' Presp and Pstim, are given, correlations on that dataset
	will be computed and displayed for each alpha.
	
	If [savallewts] is True, all weights will be returned.  Otherwise only the best weights will be
	returned.
	
	If [stop_early] is True, the weights and correlations will be returned as soon as the mean
	correlation begins to drop.  Does NOT imply early-stopping in the regularized regression sense.
	
	The given [dtype] will be applied to the regression weights as they are computed.
	"""
	## Precalculate SVD to do ridge regression
	print "Doing SVD..."
	U,S,Vh = np.linalg.svd(Rstim, full_matrices=False)
	#mU = M(U.T) ## Convert U and Vh to numpy matrices
	#mVh = M(Vh.T)
	#del U ## Delete array representations for memory's sake
	#del Vh
	Prespnorms = np.apply_along_axis(np.linalg.norm, 0, Presp) ## Precompute test response norms
	Rcorrs = []
	Pcorrs = []
	wts = []
	bestcorr = -1.0
	UR = np.dot(U.T, Rresp)
	for a in alphas:
		D = np.diag(S/(S**2+a**2)) ## Reweight singular vectors by the ridge parameter 
		#wt = (mVh*M(D)*mU*M(Rresp)).astype(dtype) ## Expand pseudoinverse with new singular values,
		#		 d								 ## multiply by response
		#wt = reduce(np.dot, [Vh.T, D, U.T, Rresp]).astype(dtype)
		wt = reduce(np.dot, [Vh.T, D, UR]).astype(dtype)
		pred = np.dot(Pstim, wt) ## Predict test responses
		prednorms = np.apply_along_axis(np.linalg.norm, 0, pred) ## Compute predicted test response norms
		#Rcorr = np.array([np.corrcoef(Presp[:,ii], pred[:,ii].ravel())[0,1] for ii in range(Presp.shape[1])]) ## Slowly compute correlations
		Rcorr = np.array(np.sum(np.multiply(Presp, pred), 0)).squeeze()/(prednorms*Prespnorms) ## Efficiently compute correlations
		Rcorr[np.isnan(Rcorr)] = 0
		Rcorrs.append(Rcorr)
		
		if saveallwts:
			wts.append(wt)
		elif Rcorr.mean()>bestcorr:
			bestcorr = Rcorr.mean()
			del wt
			#wts = wt
		
		print "Training: alpha=%0.3f, mean corr=%0.3f, max corr=%0.3f, over-under(%0.2f)=%d" % (a, np.mean(Rcorr), np.max(Rcorr), corrmin, (Rcorr>corrmin).sum()-(-Rcorr>corrmin).sum())
		
		## Test alpha on real test set if given
		if rPresp is not None and rPstim is not None:
			rpred = np.dot(rPstim, wt)
			Pcorr = np.array([np.corrcoef(rPresp[:,ii], rpred[:,ii].ravel())[0,1] for ii in range(rPresp.shape[1])])
			Pcorr[np.isnan(Pcorr)] = 0.0
			print "Testing: alpha=%0.3f, mean corr=%0.3f, max corr=%0.3f" % (a, np.mean(Pcorr), np.max(Pcorr))
			Pcorrs.append(Pcorr)
			if sum(np.isnan(Pcorr)):
				raise Exception("nan correlations")

		## Quit if mean correlation decreases
		if stop_early and Rcorr.mean()<bestcorr:
			break

	if rPresp is not None and rPstim is not None:
		return wts, Rcorrs, Pcorrs
	else:
		return wts, Rcorrs

def eigridge(Rstim, Pstim, Rresp, Presp, alphas, rPresp=None, rPstim=None, saveallwts=True,
			 stop_early=False, dtype=np.single, corrmin=0.2):
	"""Ridge regresses [Rstim] onto [Rresp] for each ridge parameter in [alpha].  Returns the fit
	linear weights for each alpha, as well as the distributions of correlations on a held-out test
	set ([Pstim] and [Presp]).  Note that these should NOT be the "real" held-out test set, only a
	small test set used to find the optimal ridge parameter.
	
	If an [rPresp] and [rPstim], or 'real' Presp and Pstim, are given, correlations on that dataset
	will be computed and displayed for each alpha.
	
	If [savallewts] is True, all weights will be returned.  Otherwise only the best weights will be
	returned.
	
	If [stop_early] is True, the weights and correlations will be returned as soon as the mean
	correlation begins to drop.  Does NOT imply early-stopping in the regularized regression sense.
	
	The given [dtype] will be applied to the regression weights as they are computed.
	"""
	## Precalculate SVD to do ridge regression
	print "Doing Eigenvalue decomposition..."
	cmode = Rstim.shape[0]<Rstim.shape[1]
	Prespnorms = np.apply_along_axis(np.linalg.norm, 0, Presp) ## Precompute test response norms
	
	if cmode:
		[S,U] = np.linalg.eigh(np.dot(Rstim.astype(dtype), Rstim.T.astype(dtype)))
		U1 = np.dot(Rstim.T, U)
		U2 = np.dot(U.T, Rresp)
	else:
		[S,U] = np.linalg.eigh(np.dot(Rstim.T.astype(dtype), Rstim.astype(dtype)))
		Usr = reduce(np.dot, [U.T, Rstim.T, Rresp])

	print "Cmode=",cmode
	print "Eigenvector shape=",U.shape
	
	Rcorrs = []
	Pcorrs = []
	wts = []
	bestcorr = -1.0
	for a in alphas:
		print "Running alpha %0.3f"%a
		D = np.diag(1/(S+a)).astype(dtype)
		
		if cmode:
			#wt = reduce(np.dot, [Rstim.T, reduce(np.dot, [U, D, U.T]).astype(dtype), Rresp]).astype(dtype)
			wt = reduce(np.dot, [U1, D, U2]).astype(dtype)
		else:
			#wt = reduce(np.dot, [reduce(np.dot, [U, D, U.T]).astype(dtype), Rstim.T, Rresp]).astype(dtype)
			wt = reduce(np.dot, [U, D, Usr]).astype(dtype)
		
		pred = np.dot(Pstim, wt) ## Predict test responses
		prednorms = np.apply_along_axis(np.linalg.norm, 0, pred) ## Compute predicted test response norms
		#Rcorr = np.array([np.corrcoef(Presp[:,ii], pred[:,ii].ravel())[0,1] for ii in range(Presp.shape[1])]) ## Slowly compute correlations
		Rcorr = np.array(np.sum(np.multiply(Presp, pred), 0)).squeeze()/(prednorms*Prespnorms) ## Efficiently compute correlations
		Rcorr[np.isnan(Rcorr)] = 0
		Rcorrs.append(Rcorr)
		
		if saveallwts:
			wts.append(wt)
		elif Rcorr.mean()>bestcorr:
			bestcorr = Rcorr.mean()
			wts = wt
		
		print "Training: alpha=%0.3f, mean corr=%0.3f, max corr=%0.3f, over-under(%0.2f)=%d" % (a, np.mean(Rcorr), np.max(Rcorr), corrmin, (Rcorr>corrmin).sum()-(-Rcorr>corrmin).sum())
		
		## Test alpha on real test set if given
		if rPresp is not None and rPstim is not None:
			rpred = np.dot(rPstim, wt)
			Pcorr = np.array([np.corrcoef(rPresp[:,ii], rpred[:,ii].ravel())[0,1] for ii in range(rPresp.shape[1])])
			Pcorr[np.isnan(Pcorr)] = 0.0
			print "Testing: alpha=%0.3f, mean corr=%0.3f, max corr=%0.3f" % (a, np.mean(Pcorr), np.max(Pcorr))
			Pcorrs.append(Pcorr)
			if sum(np.isnan(Pcorr)):
				raise Exception("nan correlations")

		## Quit if mean correlation decreases
		if stop_early and Rcorr.mean()<bestcorr:
			break

	if rPresp is not None and rPstim is not None:
		return wts, Rcorrs, Pcorrs
	else:
		return wts, Rcorrs


import random
import itertools as itools

def bootstrap_ridge(Rstim, Rresp, Pstim, Presp, alphas, nboots, chunklen, nchunks, dtype=np.single, corrmin=0.2, joined=None):
	"""Uses ridge regression with a bootstrapped held-out set to get optimal alpha values for each voxel.
	[nchunks] random chunks of length [chunklen] will be taken from [Rstim] and [Rresp] for each regression
	run.  [nboots] total regression runs will be performed.  The best alpha value for each voxel will be
	averaged across the bootstraps to estimate the best alpha for that voxel.
	
	If [joined] is given, it should be a list of lists where the STRFs for all the voxels in each sublist 
	will be given the same regularization parameter (the one that is the best on average).
	"""
	nresp, nvox = Rresp.shape

	Rcmats = []
	for bi in range(nboots):
		print "Selecting held-out test set.."
		allinds = range(nresp)
		indchunks = zip(*[iter(allinds)]*chunklen)
		random.shuffle(indchunks)
		heldinds = list(itools.chain(*indchunks[:nchunks]))
		notheldinds = list(set(allinds)-set(heldinds))
		
		RRstim = Rstim[notheldinds,:]
		PRstim = Rstim[heldinds,:]
		RRresp = Rresp[notheldinds,:]
		PRresp = Rresp[heldinds,:]
		
		## Run ridge regression using this test set
		Rwts, Rcorrs = ridge(RRstim, PRstim, RRresp, PRresp, alphas, saveallwts=False, dtype=dtype, corrmin=corrmin)
		
		Rcmat = np.vstack(Rcorrs)
		Rcmats.append(Rcmat)
		#bestainds = np.array(map(np.argmax, Rcmat.T))
		#bestalphas[bi,:] = alphas[bestainds]
	
	return np.dstack(Rcmats) # stack 2D matrices into 3D matrices

def find_best_alpha(Rresp, corrs, alphas, nboots, joined):
	nresp, nvox = Rresp.shape
	bestalphas = np.zeros((nboots, nvox))  ## Will hold the best alphas for each voxel
	
	print "Finding best alpha for each voxel.."
	if joined is None:
		## Find best alpha for each voxel
		meanbootcorrs = corrs.mean(2) # CHANGED 8/2/11 to work with VOC+VBK, CHECK THIS!!! LH
		#meanbootcorrs = corrs.mean((corrs.ndim-1))
		#print meanbootcorrs
		bestalphainds = np.argmax(meanbootcorrs, 0)
		valphas = alphas[bestalphainds]
	else:
		## Find best alpha for each group of voxels
		valphas = np.zeros((nvox,))
		for jl in joined:
			jcorrs = corrs[:,jl,:].mean(1).mean(1) ## Mean across voxels in the set, then mean across bootstraps
			bestalpha = np.argmax(jcorrs)
			valphas[jl] = alphas[bestalpha]
	return valphas
		
def set_weights(Rstim, Rresp, Pstim, Presp, alphas, valphas):
	## Find weights for each voxel
	U,S,Vh = np.linalg.svd(Rstim, full_matrices=False)
	UR = np.dot(U.T, np.nan_to_num(Rresp))
	pred = np.zeros(Presp.shape)
	wt = np.zeros((Rstim.shape[1], Rresp.shape[1]))
	for ai,alpha in enumerate(alphas):
		selvox = np.nonzero(valphas==alpha)[0]
		awt = reduce(np.dot, [Vh.T, np.diag(S/(S**2+alpha**2)), UR[:,selvox]])
		pred[:,selvox] = np.dot(Pstim, awt)
		wt[:,selvox] = awt

	## Find test correlations
	nnpred = np.nan_to_num(pred)
	corrs = np.nan_to_num(np.array([np.corrcoef(Presp[:,ii], nnpred[:,ii].ravel())[0,1] for ii in range(Presp.shape[1])]))
	return wt, corrs
	
def avg_corr_mats(corrs):
	avg_corr = sum(corrs)/len(corrs)
	return avg_corr
	
def plot_corr_quantiles(corrs, alphas, quants=(0.05, 0.25, 0.5, 0.75, 0.95)):
	"""Plots quantiles of the correlations for each alpha."""
	from matplotlib import cm
	nq = len(quants)
	colors = cm.jet(np.linspace(0, 1, nq))

	from matplotlib.pyplot import figure, show
	fig = figure()
	ax = fig.add_subplot(1,1,1)
	for qi,color in zip(quants, colors):
		ax.plot([quantile(cor, qi) for cor in corrs], 'o-', color=color)
	return ax

def quantile(values, prob):
	"""Returns the [prob] quantile of [values]."""
	sorted_values = np.sort(values)
	return sorted_values[np.max([int(prob*values.size-1), 0])]

def plot_ridge_report(wts, corrs, alphas):
	"""Plots a nice report of a ridge regression run."""
	from matplotlib.pyplot import figure, show
	fig = figure()
	ax1 = fig.add_subplot(1,1,1)
	
	#ax1.boxplot(np.vstack([c[c!=0] for c in corrs]).T, notch=1, whis=2.0)
	violin_plot(ax1, np.vstack([c[c!=0] for c in corrs]), np.arange(len(corrs)))

from matplotlib.pyplot import figure, show
from scipy.stats import gaussian_kde
from numpy.random import normal
from numpy import arange

def violin_plot(ax, data, pos, bp=False):
	"""
	create violin plots on an axis
	"""
	dist = max(pos)-min(pos)
	w = min(0.15*max(dist,1.0),0.5)
	for d,p in zip(data,pos):
		k = gaussian_kde(d) #calculates the kernel density
		m = k.dataset.min() #lower bound of violin
		M = k.dataset.max() #upper bound of violin
		x = arange(m,M,(M-m)/100.) # support for violin
		v = k.evaluate(x) #violin profile (density curve)
		v = v/v.max()*w #scaling the violin to the available space
		ax.fill_between(x,p,v+p,facecolor='y',alpha=0.3)
		ax.fill_between(x,p,-v+p,facecolor='y',alpha=0.3)
	if bp:
		ax.boxplot(data,notch=1,positions=pos,vert=1)

