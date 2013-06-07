import os, glob, h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import misc
import itertools

studydir = '/Volumes/BOB_SAGET/TNFalpha/tinnitus/unilateralstim'
subjIDs = ['Red102', 'Green102', 'Blue102', 'Red101']
ears = 'rl'
hemis = 'rl'
epochs = ['Pre', 'Post']
colors = dict(Pre = 'b', Post = 'r')

def post_lesion_comparison():

	fig = plt.figure(figsize = (15, 8));
	ax = [fig.add_subplot(2, 2, i) for i in range(4)]
	for subjID in subjIDs:

		df = combine_pre_post(subjID)
		gp = df.groupby(['ear', 'hemi', 'epoch'])
		lfp_mean = gp.lfp.apply(np.mean)

		for i, (hemi, ear) in enumerate(itertools.product(hemis, ears)):
			for epoch in epochs:
				
				l, = ax[i].plot(lfp_mean.ix[(ear, hemi, epoch)], color = colors[epoch], label = epoch)
				ax[i].set_title('%s ear / %s hemi' % (ear, hemi))

		ax[0].legend()
		fig.savefig(os.path.join(studydir, 'Analysis', '%s_post_lesion_comparison.png' % subjID))
		[a.cla() for a in fig.get_axes()]

	plt.close(fig)

def combine_pre_post(subjID):

	epochs = ['Pre', 'Post']

	df_ = []
	for epoch in epochs:
		sesspaths = glob.glob(os.path.join(studydir, 'Sessions', epoch, '*'))
		for sesspath in sesspaths:
			fpath = os.path.join(sesspath, 'both', '%s_both.h5' % subjID)
			d = pd.HDFStore(fpath, 'r')

			df_.append(d['df'])
			
			d.close()
			print sesspath, len(df_[-1])
			df_[-1]['epoch'] = epoch

	df = pd.concat(df_)
	df.index = range(len(df))

	return df

def make_contactsheets():

	subjIDs = ['Red102', 'Green102', 'Blue102', 'Red101']
	epochs = ['Pre', 'Post']

	for subjID in subjIDs:
		for epoch in epochs:
			sesspaths = glob.glob(os.path.join(studydir, 'Sessions', epoch, '*'))
			for sesspath in sesspaths:
				# fpaths = 
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

def combine(epoch = 'Pre', stimtype = 'noise'):

	sesspaths = glob.glob(os.path.join(studydir, 'Sessions', epoch, '*'))
	for sesspath in sesspaths:
		absol, sessID = os.path.split(sesspath)
		fpaths = glob.glob(os.path.join(sesspath, 'fileconversion', '*.h5'))
		subjIDs = np.unique([os.path.splitext(os.path.split(fpath)[1])[0].split('_')[0] for fpath in fpaths])

		hemis = 'rl'
		ears = 'rl'

		for subjID in subjIDs:

			LFP = []; ATTEN = []; EAR = []; HEMI = []; RELHEMI = []; SESSID = [];

			for ear in ears:

				fpath = glob.glob(os.path.join(sesspath, 'fileconversion', '%s_%s_%s_b*.h5' % (subjID, ear.upper(), stimtype)))[0]

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
					SESSID.extend([sessID]*ntrials)

			d = dict(lfp = LFP, atten = ATTEN, ear = EAR, hemi = HEMI, relhemi = RELHEMI, sess = SESSID)
			df = pd.DataFrame(d)

			fout = pd.HDFStore(os.path.join(sesspath, 'both', '%s_both.h5' % subjID))
			fout['df'] = df
			fout.close()


