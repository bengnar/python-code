from matplotlib import pyplot as plt
import h5py, os, glob, re
import numpy as np
import RF
import misc
from sp_spectral import welch


basedir = '/Volumes/BOB_SAGET/Fmr1_voc'
# basedir = '/Users/robert/Desktop'

def calc_coherence_all(experiment, stim_sep = False, pplot = False):
	
	dtype = np.dtype([('gen', 'S2'), ('exp', 'S3'), ('sess', 'S20'), ('unit', 'i4'), ('stimID', 'i4'), ('Psignal', '129f4'), ('Pnoise', '129f4'), ('Pnoise_jn' ,'129f4'), ('coh', '26f4'), ('coh_jn', '26f4'), ('F', '129f4'), ('info_snr', 'f4'), ('info_snr_jn', 'f4')])

	p = re.compile('(\d+)')
	
	_, gen, exp, _ = experiment.split('_')
	
	fpaths = glob.glob(os.path.join(basedir, experiment, 'fileconversion', 'VOC*.h5'))
	coherence_info = np.empty(0, dtype = dtype)

	for fpath in fpaths:

		print fpath
		
		absol, relat = os.path.split(fpath)
		penname, _ = os.path.splitext(relat)
		pennum = np.int32(p.findall(penname)[0])

		f = h5py.File(fpath, 'r')
		rast = f['rast'].value
		stimparams = f['stimID'].value
		stimparams = stimparams[:, 0][..., np.newaxis]
		f.close()
		
		if stim_sep is True:
			ustim = np.unique(stimparams[:, 0])
			nstim = ustim.size
		elif stim_sep is False:
			ustim = None
			nstim = 1

		for i in range(nstim):
			Psignal, Pnoise, Pnoise_jn, f, coh_snr, coh_snr_jn, info_snr, info_snr_jn = calc_coherence(rast, stimparams, ustim = [ustim[i]])

			coherence_info.resize(coherence_info.size+1)
			coherence_info[-1] = np.array((gen, exp, experiment, pennum, ustim[i], Psignal, Pnoise, Pnoise_jn, coh_snr, coh_snr_jn, f, info_snr, info_snr_jn), dtype = dtype)

	
			if pplot:
				outpath = os.path.join(basedir, 'Analysis', '%s_%s_stim%i_SNR.png' % (experiment, penname, ustim[i]))
				fig = plot_coherence(Psignal, Pnoise, Pnoise_jn, f, coh_snr, coh_snr_jn, info_snr, info_snr_jn, outpath)
				plt.close(fig)
		
	return coherence_info

def calc_coherence(rast, stimparams, ustim = None):
	'''Calculate the noise using a signal generated from all trials as well as a signal that does not incldue the trial for wich you calculate the noise (jack knifing)
	'''

	if ustim is None:
		ustim = np.unique(stimparams[:, 0])
	ustim = np.array(ustim)

	nstim = ustim.size
	npsth = 16000.
	ntrials = [(stimparams[:, 0]==x).sum() for x in ustim]
	assert np.unique(ntrials).size==1
	ntrial = ntrials[0]

	noise = np.zeros((npsth, ntrial, nstim))
	orig_signal = np.zeros((npsth, ntrial, nstim)) # raster
	noise_jn = np.zeros((npsth, ntrial, nstim))
	signal = np.zeros((npsth, nstim)) #signal ave across all the trials
	
	specsamplerate = 20
	## Loop on stimuli
	for k in range(nstim):
		
		ix = RF.get_trials(stimparams, [k+1])
		resp = rast[ix, :16000]
		ntrial, npsth = resp.shape  # Number of trials for this stim
		binsize = 1000./specsamplerate	#Sampling rate ms.
		ndur = np.round(npsth/binsize) # Length for this stim-response pair
	
		# Signal is estimated as the average of all trials for that stimulus
		signal[:, k] = resp.mean(0)	
		for itrial in range(ntrial):
		
			# Calculate the two estimates of the noise
			ix = np.arange(ntrial)!=itrial
			noise_jn[:, itrial, k] = resp[itrial, :] - resp[ix, :].mean(0) #jackknifing (exclude current trial)
			noise[:, itrial, k] = resp[itrial, :] - signal[:, k]
	
	## Calculate the noise and signal psd obtained by averaging before and after jackknifing
	nfft = 256
	fs = 1000.
	
	noise_ravel = noise.ravel(order = 'F')
	noise_jn_ravel = noise_jn.ravel(order = 'F')
	signal_ravel = signal.ravel(order = 'F')
	assert (noise_ravel[:npsth]==noise[:, 0, 0]).all()
	assert (noise_jn_ravel[:npsth]==noise_jn[:, 0, 0]).all()
	assert (signal_ravel[:npsth]==signal[:, 0]).all()
	
	f, Pnoise = welch(noise_ravel, fs = fs, nfft = nfft)
	f, Pnoise_jn = welch(noise_jn_ravel, fs = fs, nfft = nfft)
	f, Psignal = welch(signal_ravel, fs = fs, nfft = nfft)

	## Calculate and plot the coherence calculated from the SNR
	infobound = 100 # upper bound on information in Hz
	ix = f<=infobound
	coh_snr = Psignal[ix] / (Pnoise[ix] + Psignal[ix])
	df = f[1] - f[0]
	info_snr = (-np.log2(1-coh_snr)).sum()*df
	coh_snr_jn = Psignal[ix] / (Pnoise_jn[ix] + Psignal[ix])
	info_snr_jn = (-np.log2(1-coh_snr_jn)).sum()*df
	
	return Psignal, Pnoise, Pnoise_jn, f, coh_snr, coh_snr_jn, info_snr, info_snr_jn, 

def plot_coherence(Psignal, Pnoise, Pnoise_jn, f, coh_snr, coh_snr_jn, info_snr, info_snr_jn, outpath, infobound = 100):
		
	ix = f<=infobound
	## Plot SNR power across frequencies
	fig = plt.figure()
	ax1 = fig.add_subplot(211)
	ax1.plot(f, 10*np.log10(Pnoise_jn), 'r')
	ax1.plot(f, 10*np.log10(Psignal), 'k')
	ax1.legend(['Noise','Signal'])
	ax1.set_xlabel('Frequency (Hz)')
	ax1.set_ylabel('Power (dB)')
	ax1.set_title('Estimated power for Signal and Noise across Freq')

	ax2 = fig.add_subplot(212)
	ax2.plot(f[ix], coh_snr_jn, 'b')
	ax2.set_ylabel('Coherence')
	ax2.set_xlabel('Frequency (Hz)')
	#leg('Coherence SNR','Coherence SNR D1')
	ax2.legend(['Coherence'])
	ax2.set_title('Total information: UB %.f  LB %.0f bits/s' % (info_snr, info_snr_jn))
	
	fig.savefig(outpath)
	plt.show()
	return fig


def plot_info_bounds(coh_info, ax = None):

	## Plots the information upper and lower bound per penetration, averaging over stimuli
	coh_info = coh_info[coh_info['info_snr'].argsort()]
	
	ax, fig = misc.axis_check(ax)

	ax.plot(np.arange(coh_info.size), coh_info['info_snr'], 'b.') #plots upper bound
	ax.plot(np.arange(coh_info.size), coh_info['info_snr_jn'],'r.') #plots lower bound
	#line(info_mat(:,2),info_mat(:,3))
	ax.set_xticks(np.arange(coh_info.size))
	ax.set_xticklabels(coh_info['unit'])
	ax.set_xlabel('Penetration')
	ax.set_ylabel('Information Bound')
	ax.legend(['Upper Bound', 'Lower Bound'], loc = 'upper left')
	ax.set_title('Info Bounds per penetration averaging across stimuli')




