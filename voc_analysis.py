import os, glob, h5py, re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import Spikes, RF, RR
import fileconversion

# studydir = '/Volumes/BOB_SAGET/Fmr1_voc'
studydir = '/Users/robert/Desktop/Fmr1_voc'
sessions = glob.glob(os.path.join(studydir, 'Sessions', '[A-Za-z]*'))

p_int = re.compile('(\d+).h5')
p_str = re.compile('([A-Za-z]+)\d+.h5')

nrrs = 7

def rr_make_contactsheets():
	'''
	loop through all the sessions and plot the rrtfs
	'''
	fig = plt.figure();
	ax_rf = fig.add_subplot(4, 4, 8);
	ax_rfrast = fig.add_subplot(4, 4, 12);
	ax_rfpsth = fig.add_subplot(4, 4, 16);
	ax_cfrr = [fig.add_subplot(nrrs, 2, i) for i in np.arange(2, 3*nrrs, 3)]
	ax_noiserr = [fig.add_subplot(nrrs, 2, i) for i in np.arange(1, 3*nrrs, 3)]

	for session in sessions:

		unitinfos = fileconversion.get_session_unitinfo(session, onlycomplete = ('RF', 'RR', 'VOC'))
			
		for unitkey in unitinfos.iterkeys():
			
			unitinfo = unitinfos[unitkey]

			rf_ix = unitinfo['stimtype'].index('RF')
			f_rf = h5py.File(unitinfo['fpath'][rf_ix], 'r')
			rf_rast = f['rast'].value
			rf_stimparams = f['stimID'].value
			f_rf.close()
			
			''' calculate and plot RF, psth, and sorted raster'''
			rf = RF.calc_rf(rf_rast, rf_stimparams)
			rf_psth = Spikes.calc_psth(rf_rast)
			RF.plot_rf(rf, ax = ax_rf)
			Spikes.plot_sorted_raster(rf_rast, rf_stimparams, ax = ax_rfrast)

			''' calcualte and plot RRTFs for CF and noise stimuli '''
			rr_ix = unitinfo['stimtype'].index('RF')
			f_rr = h5py.File(unitinfo['fpath'][rr_ix], 'r')
			rr_rast = f_rr['rast'].value
			rr_stimparams = f_rr['stimID'].value

			rrtfs = RR.calc_rrtf_all(rr_rast, rr_stimparams)

				
		
		
		



		
		