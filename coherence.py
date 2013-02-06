import os, glob, h5py, re
import numpy as np
import nitime as nt
import nitime.viz as viz
from nitime.analysis import SNRAnalyzer

basedir = '/Volumes/BOB_SAGET/Fmr1_RR'

def calc_coherence_all(experiment, lb = 0, ub = 100, prefix = 'VOC'):
	
	_, gen, exp, _ = experiment.split('_')
	
	fpaths = glob.glob(os.path.join(basedir, experiment, 'fileconversion', '%s*.h5' % prefix))

	npsth = 16000

	p = re.compile('(\d+)')
	colnames = ['gen', 'exp', 'sess', 'unit', 'stim', 'Psignal', 'Pnoise', 'snr', 'coh', 'F', 'info_snr']
	df = pd.DataFrame(columns = colnames)

	# dtype = np.dtype([('gen', 'S2'), ('exp', 'S3'), ('sess', 'S20'), ('unit', 'i4'), ('stim', 'i4'), ('Psignal', '1600f4'), ('Pnoise', '1600f4'), ('snr', '1600f4'), ('coh', '1600f4'), ('F', '1600f4'), ('info_snr', 'f4')])
	# coh_info = np.empty(0, dtype = dtype)

	# loop through penetrations
	for fpath in fpaths[:1]:

		absol, relat = os.path.split(fpath)
		penname, _ = os.path.splitext(relat)
		pennum = np.int32(p.findall(penname)[0])

		f = h5py.File(fpath, 'r')
		rast = f['rast'].value
		stimparams = f['stimID'].value
		stimparams = stimparams[:, 0][..., np.newaxis]
	
		ustim = np.unique(stimparams[:, 0])
		nstim = ustim.size
		
		# loop through stimuli
		for i in range(nstim):

			ix = stimparams[:, 0]==ustim[i]
			rast_ = rast[ix, :npsth]
			t1 = nt.TimeSeries(rast_, sampling_interval = 0.001)
			snr1 = SNRAnalyzer(t1)
		
			F = snr1.mt_frequencies
			freq_ix = F<ub
			F = F[freq_ix]
			info = snr1.mt_information[freq_ix]
			pnoise = snr1.mt_noise_psd[freq_ix]
			psignal = snr1.mt_signal_psd[freq_ix]
			psnr = snr1.mt_snr[freq_ix]
			coh = snr1.mt_coherence[freq_ix]
			fstep = np.diff(F[:2])[0]
			info_snr = (-np.log2(1-coh)).sum()*fstep
	
			# coh_info.resize(coh_info.size+1)
			# coh_info[-1] = np.array((gen, exp, experiment, pennum, ustim[i], psignal, pnoise, psnr, coh, F, info_snr), dtype = dtype)
			data = [gen, exp, experiment, pennum, ustim[i], psignal, pnoise, psnr, coh, F, info_snr]
			df_ = dict(zip(colnames, data))
			df = df.append(df_, ignore_index = True)

	return coh_info