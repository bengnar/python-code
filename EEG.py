import numpy as np
import pandas as pd
import os, glob, h5py
from matplotlib import pyplot as plt
from scipy.signal import butter, filtfilt
import Spikes, Lfp, RF
import misc
from scipy.io import loadmat
from matplotlib.mlab import specgram

studydir = '/Volumes/BOB_SAGET/TNFalpha/tinnitus/'


def fileconvert_all(experiment = 'unilateralstim', epoch = 'Pre'):

	sesspaths = glob.glob(os.path.join(studydir, experiment, 'Sessions', epoch, '*'))
	for sesspath in sesspaths:
		matpaths = glob.glob(os.path.join(sesspath, 'data', '[A-Za-z]*.mat'))	
		for matpath in matpaths:
		    print matpath
	            fileconvert(matpath)

def fileconvert(matpath):

	outpath = matpath.replace('data', 'fileconversion').replace('.mat', '.h5')
	if not os.path.exists(outpath):
		tmpfile = loadmat(matpath)
		Data0 = tmpfile['Data']
		u_ = h5py.File(outpath, 'w')
		
		lfp, stimID = export_unit(Data0)
		# save out to file
		u_.create_dataset('lfp', data = lfp, compression = 'gzip')
		u_.create_dataset('stimID', data = stimID)
		u_.close()


def export_unit(Data0):
	
	# number of time bins to include in the LFP array
	nlfpsamp = 0
	for tt, trial in enumerate(Data0['trial'][0][0][0]):

		thislfpsamp = trial['LFP'].shape[1]
		if thislfpsamp>nlfpsamp:
			nlfpsamp = thislfpsamp
	
	ntrials = Data0['trial'][0][0][0].size # find number of trials
	nstimID = Data0['trial'][0][0][0][0]['Epoch_Value'][0].size
	
	# initialze LFP
	lfp = np.ndarray((2, ntrials, nlfpsamp), dtype = 'float32')
	# initialize frequency and attenuation IDs
	stimID = np.ndarray((ntrials, nstimID), dtype = 'float32')

	for tt in range(ntrials):
		trial = Data0['trial'][0][0][0][tt]

		thisstimID = np.float32(trial['Epoch_Value'][0])
		# get the LFP for this trial and pad it with nans so it can fit in a matrix (since some of the trials have +/-1 data point for LFP)
		for cc in range(2):
			lfpchannel = trial['LFP'][cc]
			lfp[cc, tt, :len(lfpchannel)] = lfpchannel
			
		# add to Epoch_Value
		stimID[tt, :] = thisstimID

	remID = np.array([0., 0.])
	trial_mask = RF.make_trial_mask(stimID, remID)
	lfp = lfp[:, ~trial_mask, :]
	stimID = stimID[~trial_mask, :]

	return lfp, stimID

def add_psd(experiment = 'awakeeeg', epoch = 'Pre'):
	duration = 6.
	fpaths = glob.glob(os.path.join(studydir, experiment, 'Sessions', epoch,  'fileconversion', '*.h5'))
	for fpath in fpaths:
		print fpath
		f = h5py.File(fpath)
		S = []
		nchan, ntrials, nsamp = f['lfp'].shape
		Fs = nsamp / duration
		for j in range(nchan):
			S.append([])
			for k in range(ntrials):
				F, S_, Serr = multi_taper_psd(f['lfp'][j, k, :], Fs = Fs)
				S[-1].append(S_)
		S = np.asarray(S)
		f.create_dataset('psd', data = S, compression = 'gzip')
		f.create_dataset('F', data = F)
		f.close()

def add_spectrogram(experiment = 'awakeeeg', epoch = 'Pre'):
	
	duration = 6.
	fpaths = glob.glob(os.path.join(studydir, experiment, 'Sessions', epoch,  'fileconversion', '*.h5'))
	for fpath in fpaths:
		print fpath
		f = h5py.File(fpath)
		S = []
		nchan, ntrials, nsamp = f['lfp'].shape
		Fs = nsamp / duration
		for j in range(nchan):
			S.append([])
			for k in range(ntrials):
				S_, F, T = specgram(f['lfp'][j, k, :], Fs = Fs, NFFT = 128, noverlap = 64)
				S[-1].append(S_)
		S = np.asarray(S)
		g = f.create_group('spectrogram')
		g.create_dataset('S', data = S, compression = 'gzip')
		g.create_dataset('F', data = F)
		g.create_dataset('T', data = T)
		f.close()

def remove_60hz(lfp, fs = 384.384384384):

	fs_nyq = fs / 2.

	b, a = butter(5, [59./fs_nyq, 61./fs_nyq], btype = 'bandstop', output = 'ba')
	lfp_filt = filtfilt(b, a, lfp)

	return lfp_filt

# lfp = df.lfp[0]
# lfpfilt = df.lfpfilt[0]

# fs = 
# LFP = np.fft.fftshift(np.fft.rfft(lfp))
# LFPFILT = np.fft.fftshift(np.fft.rfft(lfpfilt))

# freqs = np.fft.fftfreq(d=1./fs)
