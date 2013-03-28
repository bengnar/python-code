import os, glob, h5py, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import Spikes, RF, RR, misc
import fileconversion

# studydir = '/Volumes/BOB_SAGET/Fmr1_voc'
studydir = '/Users/robert/Desktop/Fmr1_voc'
sessions = glob.glob(os.path.join(studydir, 'Sessions', '[A-Za-z]*'))

ix2freq = 1000 * 2**(np.arange(0, 64)/10.)
nrrs = 7
t_rrtf = np.arange(0, 6, 0.001)

def rr_make_contactsheets():
	'''
	loop through all the sessions and plot the rrtfs
	'''
	fig = plt.figure();
	ax_rf = fig.add_subplot(4, 3, 6);
	ax_rfrast = fig.add_subplot(4, 3, 9);
	ax_rfpsth = fig.add_subplot(4, 3, 12);
	ax_cfrr = [fig.add_subplot(nrrs, 3, i) for i in np.arange(2, 3*nrrs, 3)]
	ax_noiserr = [fig.add_subplot(nrrs, 3, i) for i in np.arange(1, 3*nrrs, 3)]

	for session in sessions[:1]:

		unitinfos = fileconversion.get_session_unitinfo(session, onlycomplete = ('RF', 'RR', 'VOC'))
			
		for unitkey in unitinfos.keys()[:1]:
			
			unitinfo = unitinfos[unitkey]

			rf_ix = unitinfo['stimtype'].index('RF')
			
			f_rf = h5py.File(unitinfo['fpath'][rf_ix], 'r')
			rf_rast = f_rf['rast'].value
			rf_stimparams = f_rf['stimID'].value
			cf_ix = f_rf['cf'].value
			f_rf.close()
			
			cf = ix2freq[20:][int(cf_ix)]

			''' calculate and plot RF, psth, and sorted raster'''
			rf = RF.calc_rf(rf_rast, rf_stimparams)
			rf_psth = Spikes.calc_psth(rf_rast)
			RF.plot_rf(rf, ax = ax_rf)
			ax_rf.axvline(cf_ix, color = 'r', lw = 1.5)
			Spikes.plot_sorted_raster(rf_rast, rf_stimparams, ax = ax_rfrast)

			''' calcualte and plot RRTFs for CF and noise stimuli '''
			rr_ix = unitinfo['stimtype'].index('RR')
			
			f_rr = h5py.File(unitinfo['fpath'][rr_ix], 'r')
			rr_rast = f_rr['rast'].value
			rr_stimparams = f_rr['stimID'].value
			f_rr.close()

			# find the played CF
			rr_ufreqs = np.unique(rr_stimparams[:, 0])
			urrs = np.unique(rr_stimparams[:, 1])
			rr_freq, rr_ufreq_ix, _ = misc.closest(rr_ufreqs, cf, log = True)

			ax_rf.axvline(RF.freq2ix(rr_freq), color = 'g', lw = 1.5)
			# calculate the PSTHs for each repetition rate
			tmp = Spikes.calc_psth_by_stim(rr_rast, rr_stimparams)
			rrtf_cfpsth = tmp[0][rr_ufreq_ix, :, :]
			rrtf_noisepsth = tmp[0][0, :, :]

			# plot repetition rate PSTHs
			for i in xrange(nrrs):
				RR.plot_rrtf(t_rrtf, rrtf_noisepsth[i, :], urrs[i], int(4*urrs[i]), onset = 0.05, duration = 0.025, ax = ax_noiserr[i])
				RR.plot_rrtf(t_rrtf, rrtf_cfpsth[i, :], urrs[i], int(4*urrs[i]), onset = 0.05, duration = 0.025, ax = ax_cfrr[i])

			ax_noiserr[0].set_title('Noise RRTFs')
			ax_cfrr[0].set_title('CF RRTFs (%.0f kHz)' % (cf/1000))
			plt.show();

				
		
		
		



		
		