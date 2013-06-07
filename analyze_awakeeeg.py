import os, glob, h5py
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
import Lfp, Spikes, voc

import misc

studydir = '/Volumes/BOB_SAGET/TNFalpha/tinnitus/awakeeeg'
stim_duration = 3.
trial_duration = 4.5
subjIDs = ['Red102', 'Green102', 'Blue102', 'Red101']
epochs = ['Pre', 'Post']
hemis = 'rl'

def post_lesion_silent_comparison():

	for subjID in subjIDs:
		for epoch in epochs:
			d = pd.HDFStore(os.path.join(studydir, 'Sessions', epoch, 'fileconversion', '%s_silent.h5' % subjID))
			df.append(d['df'])
			d.close()

	
def make_silent_sheets(epoch = 'Pre'):

	sesspaths = glob.glob(os.path.join(studydir, 'Sessions', epoch, '*'))
	for sesspath in sesspaths:
		fpaths = glob.glob(os.path.join(sesspath, 'fileconversion', '*.h5'))

		for fpath in fpaths:

			absol, relat = os.path.split(fpath)
			fname, _ = os.path.splitext(relat)
			
			f = pd.HDFStore(fpath, 'r')
			df = f['df']
			f.close()

			df = df[df.rr==-1]
			df = df[df.good]
			gp = df.groupby('hemi')

			nsamp = df.lfp.values[0].size
			Fs = nsamp / trial_duration

			lfp_mean = gp.lfp.apply(np.mean)
			lfp_err = gp.lfp.apply(np.std)

			fft_mean = gp.fft.apply(np.mean)
			fft_err = gp.fft.apply(np.std)

			t = np.linspace(0, trial_duration, nsamp)
			f = np.linspace(0, Fs/2., df.fft.values[0].size)

			fig = plt.figure()
			ax_lfp = []; ax_fft = []
			for i, ((k, lfp_mean_), (k, lfp_err_), (k, fft_mean_), (k, fft_err_)) in enumerate(zip(lfp_mean.iterkv(), lfp_err.iterkv(),\
				fft_mean.iterkv(), fft_err.iterkv())):
				# lfp_err_ = lfp_err[k]
				ax_lfp.append(fig.add_subplot(2, 2, i+1));
				misc.errorfill(t, lfp_mean_, lfp_err_, ax = ax_lfp[-1])
				ax_lfp[-1].set_title(k)

				ax_fft.append(fig.add_subplot(2, 2, i+3));
				misc.errorfill(f, np.log(fft_mean_), np.log(fft_err_), ax = ax_fft[-1])

			misc.sameyaxis(ax_lfp); misc.sameyaxis(ax_fft)
			misc.samexaxis(ax_lfp, [0, stim_duration]); misc.samexaxis(ax_fft, [f.min(), f.max()])

			figpath = os.path.join(studydir, 'Sheets', '%s_silent.png' % fname)
			fig.savefig(figpath)
			plt.close(fig)

def make_noise_rr_sheets():

	fpaths = glob.glob(os.path.join(studydir, 'Sessions', 'Pre', 'both',  '*.h5'))

	for fpath in fpaths:

		absol, relat = os.path.split(fpath)
		fname, _ = os.path.splitext(relat)

		f = pd.HDFStore(fpath, 'r')
		df = f['df']
		f.close()
		df = df[df.rr>0] # subset only non-silent trials
		df = df[df.good] # use only good trials
		urrs = np.unique(df.rr)
		nrrs = urrs.size

		gp = df.groupby(('rr', 'hemi'))

		nsamp = df.lfp.values[0].size
		Fs = nsamp / trial_duration

		lfp_mean = gp.lfp.apply(np.mean) 
		lfp_err = gp.lfp.apply(np.std)

		fft_mean = gp.fft.apply(np.mean)
		fft_err = gp.fft.apply(np.std)

		t = np.linspace(0, trial_duration, nsamp)
		f = np.linspace(0, Fs/2., df.fft.values[0].size)

		fig = plt.figure(figsize = (14, 8.8))
		ax_lfp = []; ax_fft = []
		for i, ((k, lfp_mean_), (k, lfp_err_), (k, fft_mean_), (k, fft_err_)) in enumerate(zip(lfp_mean.iterkv(), lfp_err.iterkv(),\
			fft_mean.iterkv(), fft_err.iterkv())):

			rr, hemi = k
			if hemi=='l':
				j=1
			else: j=2

			ax_lfp.append(fig.add_subplot(nrrs, 4, ((i/2)*4)+j))
			misc.errorfill(t, lfp_mean_, lfp_err_, ax = ax_lfp[-1])
			ax_lfp[-1].set_title(k)

			ax_fft.append(fig.add_subplot(nrrs, 4, ((i/2)*4)+2+j))
			misc.errorfill(f, np.abs(fft_mean_), np.abs(fft_err_), ax = ax_fft[-1])


		misc.sameyaxis(ax_lfp); misc.sameyaxis(ax_fft)
		misc.samexaxis(ax_lfp, [0, stim_duration]); misc.samexaxis(ax_fft, [0, 100.])
		[a.set_yticklabels('') for a in ax_lfp[1::2]]
		fig.savefig(os.path.join(studydir, 'Sheets', '%s_noise.png'%fname))
		plt.close(fig)

def make_noise_singleburst_sheets():

	fpaths = glob.glob(os.path.join(studydir, 'Sessions', 'Pre', 'both',  '*.h5'))

	fig = plt.figure(figsize = (14, 8.8))
	ax_lfp = []

	for j, fpath in enumerate(fpaths):

		absol, relat = os.path.split(fpath)
		fname, _ = os.path.splitext(relat)

		f = pd.HDFStore(fpath, 'r')
		df = f['df']
		f.close()
		df = df[np.logical_and(df.rr>0, df.rr<16, df.good)] # subset only non-silent trials

		gp = df.groupby('hemi')

		nsamp = df.lfp.values[0].size
		Fs = nsamp / trial_duration

		lfp_mean = gp.lfp.apply(np.mean) 
		lfp_err = gp.lfp.apply(np.std)

		fft_mean = gp.fft.apply(np.mean)
		fft_err = gp.fft.apply(np.std)

		t = np.linspace(0, trial_duration, nsamp)
		f = np.linspace(0, Fs/2., df.fft.values[0].size)

		for i, ((k, lfp_mean_), (k, lfp_err_), (k, fft_mean_), (k, fft_err_)) in enumerate(zip(lfp_mean.iterkv(), lfp_err.iterkv(),\
			fft_mean.iterkv(), fft_err.iterkv())):

			ax_lfp.append(fig.add_subplot(len(fpaths), 2, (j*2)+i+1))
			misc.errorfill(t, lfp_mean_, lfp_err_, ax = ax_lfp[-1])
			if i==0:
				ax_lfp[-1].set_title('%s / %s' % (fname, k))
			else:
				ax_lfp[-1].set_title(k)

	misc.sameyaxis(ax_lfp, [-0.0004, 0.0004])
	misc.samexaxis(ax_lfp, [0.05, 0.175])
	[a.set_yticklabels('') for a in ax_lfp[1::2]]
	fig.savefig(os.path.join(studydir, 'Sheets', 'noise_singleburst.png'))
	plt.close(fig)

def add_ffts():

	fpaths = glob.glob(os.path.join(studydir, 'Sessions', 'Pre', 'both', '*.h5'))
	for fpath in fpaths:
		f = pd.HDFStore(fpath)

		df = f['df']
		df = add_fft(df)
		f['df'] = df

		f.close()

def add_fft(df):

	ntrials = len(df)

	lfp = misc.objectarray2floatarray(df.lfp)
	
	lfp_fft = np.abs(np.fft.rfft(lfp))
	FFT = []
	[FFT.append(lfp_fft[i, :]) for i in xrange(ntrials)]
	df['fft'] = FFT

	return df

def plot_signal_quality():
	'''
	Plot "signal quality" measure for each animal.
	This is simply the histogram of maximum trial amplitudes.
	Movement artifacts are usually large amplitude (>0.002) spikes in the signal.
	'''
	fpaths = glob.glob(os.path.join(studydir, 'Sessions', 'Pre', 'both', '*.h5'))

	axs = []
	fig = plt.figure()
	npaths = len(fpaths)
	nrows = np.ceil(np.sqrt(npaths))
	for i, fpath in enumerate(fpaths):
		d = pd.HDFStore(fpath, 'r')
		df = d['df']
		d.close()
		axs.append(fig.add_subplot(nrows, nrows, i+1))
		axs[-1].hist(df.lfp.apply(np.max), bins = np.arange(0, 0.01, 0.0005))
		axs[-1].set_title(os.path.split(fpath)[-1])

	misc.samexaxis(axs)
	plt.show()
	fig.savefig(os.path.join(studydir, 'Analysis', 'signal_quality.png'))


def mark_good_trials(epoch = 'Pre'):

	sesspaths = glob.glob(os.path.join(studydir, 'Sessions', epoch, '*'))
	for sesspath in sesspaths:
		fpaths =  glob.glob(os.path.join(sesspath, 'fileconversion', '*.h5'))

		for fpath in fpaths:
		
			f = pd.HDFStore(fpath, 'r')

			df = f['df']
			df = mark_good(df)
			f['df'] = df

			f.close()

def mark_good(df):

	df['lfp_max'] = df.lfp.apply(np.max)
	df['good'] = df.lfp_max<0.0015

	return df

def plot_trials(df):

	plt.close('all');
	fig1 = plt.figure();
	fig2 = plt.figure();
	ax1 = []; ax2 = []
	hemis = ['r', 'l']
	for j, hemi in enumerate(hemis[0]):
		df_ = df[df.hemi==hemi]
		ntrials = df_.lfp.size
		nsamp = df_.lfp[df_.index[0]].size
		for k, i in enumerate(range(ntrials)[::40]):
			ax1.append(fig1.add_subplot(5, 5, k+1))
			ix = df_.index[i]
			ax1[-1].hist(df_.lfp[ix], bins = 100)
			ax1[-1].set_xticklabels('')
			ax1[-1].set_yticklabels('')
			ax1[-1].set_title(ix)

			ax2.append(fig2.add_subplot(5, 5, k+1))
			ix = df_.index[i]
			ax2[-1].plot(df_.lfp[ix])
			ax2[-1].set_xticklabels('')
			ax2[-1].set_yticklabels('')

	misc.samexaxis(ax1)
	misc.sameyaxis(ax2)

def convert_to_pd(epoch = 'Pre'):

	sesspaths = glob.glob(os.path.join(studydir, 'Sessions', epoch, '*'))

	for sesspath in sesspaths:

		absol, sessID = os.path.split(sesspath)
		fpaths = glob.glob(os.path.join(sesspath, 'fileconversion', '*.h5'))

		for fpath in fpaths:

			absol, relat = os.path.split(fpath)
			sessID, _ = os.path.splitext(relat)
			subjID, blocknum, stimtype = sessID.split('_')

			LFP = []; STIMID = []; HEMI = []; SUBJID = []; EPOCH = []; SESSID = [];
			if stimtype == 'noise':
				ATTEN = [];

			fin = h5py.File(fpath, 'r')
			lfp = fin['lfp'].value
			stimparams = fin['stimID'].value
			fin.close()

			nchan, ntrials, nsamp = lfp.shape

			for i in xrange(nchan):

				[LFP.append(lfp[i, j, :]) for j in xrange(ntrials)]
				if stimtype == 'noise':
					ATTEN.extend(stimparams[:, 1].tolist())
				HEMI.extend([hemis[i]]*ntrials)
				SUBJID.extend([subjID]*ntrials)
				EPOCH.extend([epoch]*ntrials)
				SESSID.extend([sessID]*ntrials)

			if stimtype == 'noise':
				d = dict(lfp = LFP, atten = ATTEN, hemi = HEMI, subj = SUBJID, epoch = EPOCH, sess = SESSID)
			elif stimtype == 'silent':
				d = dict(lfp = LFP, hemi = HEMI, subj = SUBJID, epoch = EPOCH, sess = SESSID)
			df = pd.DataFrame(d)

			outpath = fpath.replace('fileconversion', 'both')
			fout = pd.HDFStore(outpath)
			fout['df'] = df
			fout.close()

