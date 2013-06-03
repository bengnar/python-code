import os, glob, h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import misc

studydir = '/Volumes/BOB_SAGET/TNFalpha/tinnitus/unilateralstim'

def make_contactsheets():

	fpaths = glob.glob(os.path.join(studydir, 'Sessions', 'Pre', 'both', '*.h5'))
	for fpath in fpaths:

		absol, relat = os.path.split(fpath)
		fname, _ = os.path.splitext(relat)
		print fname

		fin = pd.HDFStore(fpath, 'r')
		df = fin['df']
		fin.close()

		ntrials = 100
		duration = 0.333
		nattens = np.unique(df.atten).size

		gp = df.groupby(('atten', 'hemi', 'relhemi'))
		means = gp.lfp.apply(np.mean)
		errs = gp.lfp.apply(np.std) / float(np.sqrt(ntrials))

		t = np.linspace(0, duration, df.lfp[0].size)

		fig = plt.figure()
		ax = []
		for i, ((k, m), (k, er)) in enumerate(zip(means.iterkv(), errs.iterkv())):

			ax.append(fig.add_subplot(nattens, 4, i+1))
			misc.errorfill(t[:-1], m[:-1], er[:-1], ax = ax[-1])
			ax[-1].set_title(k)
			ax[-1].set_xlim([0, duration])

		misc.sameyaxis(ax)
		[a.set_xticklabels('') for a in ax[:-1]]
		[a.set_yticklabels('') for a in ax[1:]]
		fig.tight_layout()
		fig.savefig(os.path.join(studydir, 'Sheets', '%s_sheet.png' % fname))
		plt.close(fig)

def combine():

	fpaths = glob.glob(os.path.join(studydir, 'Sessions', 'Pre', 'fileconversion', '*.h5'))
	subjIDs = np.unique([os.path.splitext(os.path.split(fpath)[1])[0].split('_')[0] for fpath in fpaths])

	hemis = 'rl'
	ears = 'rl'

	for subjID in subjIDs:

		LFP = []; ATTEN = []; EAR = []; HEMI = []; RELHEMI = [];

		for ear in ears:

			fpath = glob.glob(os.path.join(studydir, 'Sessions', 'Pre', 'fileconversion', '%s_b*_%s.h5' % (subjID, ear.upper())))[0]

			fin = h5py.File(fpath, 'r')
			lfp = fin['lfp'].value
			stimparams = fin['stimID'].value
			fin.close()

			nchan, ntrials, nsamp = lfp.shape

			for i in xrange(nchan):

				# LFP.extend(lfp[i, ...].tolist())
				[LFP.append(lfp[i, j, :]) for j in xrange(ntrials)]
				ATTEN.extend(stimparams[:, 1].tolist())
				EAR.extend([ear]*ntrials)
				HEMI.extend([hemis[i]]*ntrials)
				if ear==hemis[i]:
					relhemi = 'ipsi'
				elif ear!=hemis[i]:
					relhemi = 'contra'
				RELHEMI.extend([relhemi]*ntrials)

		d = dict(lfp = LFP, atten = ATTEN, ear = EAR, hemi = HEMI, relhemi = RELHEMI)
		df = pd.DataFrame(d)

		fout = pd.HDFStore(os.path.join(studydir, 'Sessions', 'Pre', 'both', '%s_both.h5' % subjID))
		fout['df'] = df
		fout.close()


