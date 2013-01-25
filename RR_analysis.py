import numpy as np
import h5py, glob, os, shutil
import Spikes; reload(Spikes)
import RF; reload(RF)
import RR; reload(RR)
import matplotlib.pyplot as plt
import misc; reload(misc);
import scipy.stats as st
import re
import shutil
import itertools

# experiments = ['ko_20120113', 'ko_20120116', 'ko_20120202', 'ko_20120206', 'ko_20120209', \
# 			'ko_20120305', 'ko_20120306', \
#  			'wt_20120123', 'wt_20120125', 'wt_20120126', 'wt_20120127', 'wt_20120203', 'wt_20120208']

# experiments = ['ko_20120206', 'ko_20120209', 'ko_20120305', 'ko_20120306', \
# 			'wt_20120126', 'wt_20120203', 'wt_20120208']

experiments = ['wt_exp_20120314', 'wt_nai_20120208', \
'wt_nai_20120410', 'wt_exp_20120320', 'wt_nai_20120412', \
'ko_exp_20120315', 'ko_nai_20120202', 'ko_nai_20120305', \
'ko_exp_20120321', 'ko_nai_20120206', 'ko_nai_20120306', \
'ko_exp_20120326', 'ko_exp_20120420', 'ko_nai_20120209']

basedir = '/Volumes/BOB_SAGET/Fmr1_RR/Sessions'
# basedir = '/Users/robert/Desktop/Fmr1_RR/Sessions'
nbins = 4200
npips = 6
onset = 0.05
wind = np.hamming(10); wind = wind / wind.sum()
ix2freq = 1000 * 2**(np.arange(0, 64)/10.)

# onset and offset times of first pip
pip_start = np.int32((onset + 0.005) * 1000)
pip_end = pip_start + 20

# dtype = np.dtype([('group', 'S2'), ('sess', 'S20'), ('unit', 'i8'), ('cf', 'f4'), \
# 	('vs_noise', 'f4'), ('vs_tone', 'f4'), ('p_noise', 'f4'), ('p_tone', 'f4'), \
# 	('dist_il_noise', 'f4'), ('dist_il_tone', 'f4'), \
# 	('dist_in_noise', 'f4'), ('dist_in_tone', 'f4'), \
# 	('power_noise', 'f4'), ('power_tone', 'f4'), \
# 	('rrtf_noise', 'f4'), ('rrtf_tone', 'f4'), \
# 	('rr', 'f4'), ('base_mean', 'f4'), ('evoked_mean', 'f4')])
	
dtype = np.dtype([('gen', 'S2'), ('exp', 'S3'), ('sess', 'S20'), ('unit', 'i8'), \
	('cf', 'f4'), ('bw', '8f4'), ('thresh', 'f4'), ('coord', '2f4'), \
	('rrtf_noise', '7f4'), ('rrtf_tone', '7f4'), \
	('rrtf_noise_norm', '7f4'), ('rrtf_tone_norm', '7f4'), \
	('rr', '7f4'), ('base_mean', 'f4'), ('evoked_mean', 'f4'), \
	('ev_psth', '333f4')])
	
def characterize(experiments = experiments, pplot = True, verbose = False):

	if type(experiments) == str:
		experiments = [experiments]
		
	# set up figure
	if pplot:
		figsize = (19, 12)
		fig = plt.figure(figsize = figsize)

	# loop through experiments
	for experiment in experiments:

		DB = np.empty(0, dtype = dtype)
		
		print '%s\n%s\n\n' % (experiment, '-'*50)
		
		# build the output directory path
		savedir = os.path.join(basedir, experiment, 'analysis')
		if not os.path.exists(savedir):
			os.mkdir(savedir)
			
		# WT or KO / CTL or EXP
		gen, exp, date = experiment.split('_')

		# find the repetition rate blocks
		rr_pens = glob.glob(os.path.join(basedir, experiment, 'fileconversion', 'RR*.h5'))

		# load the cfs for this experiment
		cfs = np.loadtxt(os.path.join(basedir, experiment, 'cfs.txt'), ndmin = 1)
		
		# loop through blocks in this experiment
		for rr_pen in rr_pens:
			
			# set up the figure for plotting RF and psth
			if pplot:
				rf_ax = fig.add_subplot(3, 4, 4)
				# build axis for psth
				pos_ = rf_ax.get_position().bounds
				pos = (pos_[0], pos_[1]-0.15, pos_[2], 0.15)
				psth_ax = fig.add_axes(pos)

			# find the matching RF file
			absol, rel = os.path.split(rr_pen)
			p = re.compile('RR')
			rf_rel = p.sub('RF', rel)
			rf_pen = os.path.join(absol, rf_rel)
			blockname = os.path.splitext(os.path.split(rr_pen)[1])[0]
			p = re.compile('(\d+).h5')
			unitnum = np.int32(p.findall(rel))[0] # unit number
			if verbose:
				print rr_pen
				print rf_pen
			
			# load the repetition rate block
			f_rr = h5py.File(rr_pen, 'r')
			spktimes = f_rr['spktimes'].value
			if not np.isnan(spktimes[0]):

				'''--------RF--------'''
				# load the RF block to get the RF
				f_rf = h5py.File(rf_pen, 'r')
				rf = f_rf['rf'].value
				rf_rast = f_rf['rast'].value
				rf_spktimes = f_rf['spktimes'].value
				rf_stimparams = f_rf['stimID'].value
				rf_spktrials = f_rf['spktrials'].value
				coord = f_rf['coord'].value
				ntrials = f_rf['rast'].shape[0]
				f_rf.close()
				rf_psth = Spikes.calc_psth(rf_spktimes, nbins = 334, ntrials = ntrials, normed = True)
				
				# plot the RF psth
				if pplot:
					psth_ax.plot(np.convolve(rf_psth, wind, 'same'))
					# psth_ax.set_ylim([0, rf_psth.max()+0.1*rf_psth.max()])
					misc.vline(psth_ax, [50, 75], color = 'r', ls = '--')

				'''------RepRate-------'''
				# finish loading the reprate block
				spktrials = f_rr['spktrials'].value
				stimparams = f_rr['stimID'].value
				ntrials_rr = f_rr['stimID'].value.shape[0]
				rast = f_rr['rast'].value
				f_rr.close()
				

				# calculate the psths
				psths, ustimparams = Spikes.calc_psth_by_stim(rast, stimparams, normed = True)
				psths = psths[:, :7, :]
				if (ustimparams[1]==0).sum() > 0:
					print 'HOLF DYDHSISHIIIITTT!!1'
				
				# get stimulus parameters
				urrs = np.empty(7) * np.nan
				urrs_ = np.unique(ustimparams[1][:7])
				urrs[:urrs_.size] = urrs_
				nrrs = 7
			
				# calculate the baseline firing rate
				base_mean = np.mean(psths[:, :, :np.int32(onset*1000)])
				# calculate the evoked firing rate
				evoked_mean = np.mean(psths[:, :, 58:88])

				# # initialize analysis measurements
				# vs_all = np.zeros((2, nrrs)) * np.nan
				# p_all = np.zeros((2, nrrs)) * np.nan 
				

				rrtf_noise = np.zeros(nrrs) * np.nan
				rrtf_tone = np.zeros(nrrs) * np.nan
				rrtf_noise_norm = np.zeros(nrrs) * np.nan
				rrtf_tone_norm = np.zeros(nrrs) * np.nan
				if pplot:
					axs = list()
					axs2 = list()

				# plot psths (noise)
				if ustimparams[0][0] == 0:
				
					# calculate the mean response to the first pip (used to calculate the RRTF)
					ix = stimparams[:, 0] == 0
					resp_1st = rast[ix, pip_start:pip_end].mean()
					# RepRate loop
					
					'''
					Loop through PSTHs for each of the noise RRs
					'''
					for i, psth in enumerate(psths[0, :, :]):
				
						# # Vector stregth NOISE
						spktimes_, spktrials_ = RF.get_spktimes(spktimes, spktrials, stimparams, [0., urrs[i]])
						# vs_all[0, i], p_all[0, i] = RR.calc_vs(spktimes_, urrs[i], npips, onset, analysistype = 'mean', rem_first = True)
						# # vs_all[0, i] = RR.calc_r_from_psth(psth, urrs[i], npips, onset, remove_first = True, base_mean = base_mean)

						# RRTF NOISE
						ix = RF.get_trials(stimparams, [0, urrs[i]])
						rast_ = rast[ix, :]
						rrtf_noise[i] = RR.calc_rrtf(rast_, urrs[i], resp_1st)
						rrtf_noise_norm[i] = RR.calc_rrtf_norm(rast_, urrs[i], resp_1st)
						
						psth_smoo = np.convolve(psth, wind, 'same')
						# # FFT following NOISE
						# power_noise[i] = RR.calc_fft_following(psth_smoo, urrs[i], npips)
						
						# van rossum distance (intertrain) NOISE
						trials = np.unique(spktrials_)
						dist_in_noise[i] = RR.calc_vR_dist_intertrain(spktimes_, spktrials_, trials, urrs[i], nbins)
						
						# van rossum distance (intertrial) NOISE
						dist_il_noise[i] = RR.calc_vR_dist_intertrial(spktimes_, spktrials_, trials, urrs[i], nbins)
						
						if pplot:					
							ax = fig.add_subplot(nrrs, 3, 3*i+1);
							ax.plot(psth_smoo);
							RR.plot_tone_pips(ax, urrs[i], npips, onset, 0.025)
							ax.set_xlim([50-200/urrs[i], 50+(1000*npips/urrs[i])])
							axs.append(ax)
							if i == 0:
								ax.set_title(unitnum)
							pos_ = ax.get_position().bounds
							pos = (pos_[0]+pos_[2]-0.03, pos_[1], pos_[2]/3., pos_[3])
							ax2 = fig.add_axes(pos, polar = True);
							ax2.set_xticks([]);
							ax2.set_yticks([])
							RR.circ_psth(ax2, spktimes_, urrs[i], npips, onset, bins = 20, rem_first = True)
							axs2.append(ax2)
					
							# if i == psths.shape[1]-1:
							# 						
							# 	# check for speaker noise
							# 	ax3 = fig.add_subplot(4, 4, 16);
							# 	# ax3.cla();
							# 	apsth = RR.aligned_psth(spktimes_, urrs[i], npips, onset)
							# 	ax3.plot(np.convolve(apsth, wind, 'same'))
							# 	ax3.axvline(5, color = 'r', ls = '--')
							# 	ax3.axvline(30, color = 'r', ls = '--')
				
					if pplot:
						misc.sameyaxis(axs)
						misc.sameyaxis(axs2)

				ix = cfs[:, 0] == unitnum
				
				# if the unit has a manually labeled CF
				if ix.sum() > 0:
					cf = cfs[ix, 1][0]
					# plot RF
					if pplot:
						rf_ax.imshow(rf, aspect = 'auto', interpolation = 'nearest')
						rf_ax.axvline(cf, color = 'w', lw = 5)

					ufreqs = np.unique(stimparams[:, 0])
					if verbose:
						print ufreqs
					# which of the RR frequencies corresponds to this unit's CF?
					freq_played, closest_f_ix = misc.closest(ufreqs[1:], ix2freq[20:][cf], log = True)
					freq_played = np.round(freq_played, 1)
					
					if verbose:
						print freq_played
					closest_f_ix =  closest_f_ix[0] + 1
		
					_, freq_played_ix = misc.closest(ix2freq[20:], freq_played, log = True)
		
					# make sure the marked cf and played freq don't differ by more than a half octave
					if np.abs(freq_played_ix - cf) < 5:
		
						if pplot:
							rf_ax.set_title('Tone @ %i' % freq_played[0])
							rf_ax.axvline(freq_played_ix, color = 'k', lw = 5)

							# plot psths (tone)
							axs = list()
							axs2 = list()
						
						# calculate the mean response to the first pip (used to calculate the RRTF)
						ix = stimparams[:, 0] == freq_played
						resp_1st = rast[ix, pip_start:pip_end].mean()
				
						'''RR -- TONES'''
						for i, psth in enumerate(psths[closest_f_ix, :, :]):
					
							# # VS TONES
							spktimes_, spktrials_ = RF.get_spktimes(spktimes, spktrials, stimparams, [freq_played, urrs[i]])
							# vs_all[1, i], p_all[1, i] = RR.calc_vs(spktimes_, urrs[i], npips, onset, analysistype = 'mean', rem_first = True)
							# # vs_all[1, i] = RR.calc_r_from_psth(psth, urrs[i], npips, onset, remove_first = True, base_mean = base_mean)

							# RRTF TONES
							ix = RF.get_trials(stimparams, [freq_played, urrs[i]])
							rast_ = rast[ix, :]
							rrtf_tone[i] = RR.calc_rrtf(rast_, urrs[i], resp_1st)
							rrtf_tone_norm[i] = RR.calc_rrtf_norm(rast_, urrs[i], resp_1st)

							psth_smoo = np.convolve(psth, wind, 'same')
							# # FFT following TONES
							# power_tone[i] = RR.calc_fft_following(psth_smoo, urrs[i], npips)

							# van rossum distance (intertrain) TONE
							trials = RF.get_trials(stimparams, [freq_played, urrs[i]])
							dist_in_tone[i] = RR.calc_vR_dist_intertrain(spktimes_, spktrials_, trials, urrs[i], nbins)
						
							# van rossum distance (intertrial) TONE
							dist_il_tone[i] = RR.calc_vR_dist_intertrial(spktimes_, spktrials_, trials, urrs[i], nbins)

							# plot VS
							if pplot:
					
								# plot PSTHs
								ax = fig.add_subplot(nrrs, 3, 3*i+2)
								ax.plot(psth_smoo, 'g')
								RR.plot_tone_pips(ax, urrs[i], npips, onset, 0.025)
								ax.set_xlim([50-200/urrs[i], 50+(1000*npips/urrs[i])])
								axs.append(ax)
					
								# plot circular histogram for tones
								pos_ = ax.get_position().bounds # build axis
								pos = (pos_[0]+pos_[2]-0.03, pos_[1], pos_[2]/3., pos_[3]) # build axis
								ax2 = fig.add_axes(pos, polar = True); # set axis
								ax2.set_xticks([]);
								ax2.set_yticks([])
								RR.circ_psth(ax2, spktimes_, urrs[i], npips, onset, bins = 20, color = 'g', rem_first = True)
								axs2.append(ax2)
					
						# plot VS/RRTF for this unit
						if pplot:
							# set psth and circular psth axes to the same range
							misc.sameyaxis(axs)
							misc.sameyaxis(axs2)
						
							ax = fig.add_subplot(2, 4, 8);
							# plot VS
							# ax.plot(urrs, vs_all.T, 'o-')
							# plot RRTF
							ax.plot(urrs, rrtf_noise_norm, 'bd-')
							ax.plot(urrs, rrtf_tone_norm, 'gd-')
						
							ax.plot(urrs, rrtf_noise, 'bd--')
							ax.plot(urrs, rrtf_tone, 'gd--')
							ax.set_ylim([0., 1.3])
							ax.axhline(1.0, color = 'r', ls = '--')
				
						# for i in range(urrs.size):
						DB.resize(DB.size + 1)
						DB[-1] = np.array((gen, exp, experiment, unitnum, cf, bw, thresh, coord, \
							dist_il_noise, dist_il_tone, dist_in_noise, dist_in_tone, \
							rrtf_noise, rrtf_tone, rrtf_noise_norm, rrtf_tone_norm, urrs, base_mean, evoked_mean, \
							ev_psth[:333]), dtype = dtype)

				if pplot:
					figpath = os.path.join(savedir, blockname + '.png')
					fig.savefig(figpath);
					fig.clf()
				
				# end unit loop
		
		np.savez(os.path.join(basedir, experiment, experiment + '_DB.npz'), DB = DB)
		if verbose:
			print '\n'*4
		# end experiment loop		
		
			# except IOError:
			# 	print 'Error!'
	
def analyze_DB(DB, p_val = 0.01):

	fig = plt.figure();
	ax1 = fig.add_subplot(221);
	ax2 = fig.add_subplot(223);
	ax3 = fig.add_subplot(222);
	ax4 = fig.add_subplot(224);
	ax1.set_title('Tone @ BF')
	ax2.set_title('Noise')
	ax2.set_ylabel('Vector strength')
	ax2.set_xlabel('Repetition rate (pips per second)')
	ax4.set_ylabel('Proportion significant responders (p<%3.3f)' % p_val)
	
	groups = np.unique(DB['group'])
	ngroups = len(groups)
	rrs = np.unique(DB['rr'])
	nrrs = rrs.size
	clrs = 'bg'
	xdat = np.zeros((ngroups, nrrs)); xstd = np.zeros((ngroups, nrrs))
	xdat_noise = np.zeros((ngroups, nrrs)); xstd_noise = np.zeros((ngroups, nrrs))
	nresp_tone = np.zeros((ngroups, nrrs)); n_tone = np.zeros((ngroups, nrrs))
	nresp_noise = np.zeros((ngroups, nrrs)); n_noise = np.zeros((ngroups, nrrs))

	for g, group in enumerate(groups):
		print group
		for r, rr in enumerate(rrs):
			ix = np.logical_and(DB['rr'] == rr, DB['group'] == group)
			dat = DB['vs_tone'][ix]
			xdat[g, r] = st.nanmean(dat.flatten())
			xstd[g, r] = st.sem(dat.flatten())
			dat = DB['vs_noise'][ix]
			xdat_noise[g, r] = st.nanmean(dat.flatten())
			xstd_noise[g, r] = st.sem(dat.flatten())
	
			dat = DB['p_tone'][ix]
			nresp_tone[g, r] = (dat < p_val).sum()
			n_tone[g, r] = dat.size
			dat = DB['p_noise'][ix]
			nresp_noise[g, r] = (dat < p_val).sum()
			n_noise[g, r] = dat.size
			
		ax1.errorbar(rrs+g/6., xdat[g, :], yerr = xstd[g, :], color = clrs[g], marker = 'o')
		ax2.errorbar(rrs+g/6., xdat_noise[g, :], yerr = xstd_noise[g, :], color = clrs[g], marker = 'o')
		ax3.plot(rrs, nresp_tone[g, :] / n_tone[g, :], color = clrs[g], marker = 'o')
		ax4.plot(rrs, nresp_noise[g, :] / n_noise[g, :], color = clrs[g], marker = 'o')

	
	# # calc stats (vector strength)
	# pval_tone= np.zeros(nrrs)
	# pval_
	
	# calc stats (proportion of responders)
	pval_resp_tone = np.zeros(nrrs)
	pval_resp_noise = np.zeros(nrrs)
	for r, rr in enumerate(rrs):
		contingency_tone = [[nresp_tone[0, r], nresp_tone[1, r]], [n_tone[0, r], n_tone[1, r]]]
		contingency_noise = [[nresp_noise[0, r], nresp_noise[1, r]], [n_noise[0, r], n_noise[1, r]]]
		[_, p] = st.fisher_exact(contingency_tone, 'two-sided')
		pval_resp_tone[r] = p
		[_, p] = st.fisher_exact(contingency_noise, 'two-sided')
		pval_resp_noise[r] = p
	
	sig = np.int32(rrs[(pval_resp_tone<p_val).nonzero()[0]])
	ax3.plot(sig, np.ones(sig.size)*0.97, 'r*', ms = 8)
	misc.vline(ax3, sig, ls = '--', color = 'r')
	sig = np.int32(rrs[(pval_resp_noise<p_val).nonzero()[0]])
	ax4.plot(sig, np.ones(sig.size)*0.97, 'r*', ms = 8)
	misc.vline(ax4, sig, ls = '--', color = 'r')
	
	ax2.legend(groups)
	misc.sameyaxis([ax1, ax2])
	ax3.set_ylim([0, 1])
	ax4.set_ylim([0, 1])

	plt.show()
	
	return xdat, xstd, rrs
	
def load_DB(experiments = experiments):
	
	DB = np.empty(0, dtype = dtype)
	for e in experiments:
		print e
		db = np.load(os.path.join(basedir, e, e + '_DB.npz'))['DB']
		for db_ in db:
			DB.resize(DB.size + 1)
			DB[-1] = db_
			
	return DB

def load_DBs(groups, suffix = 'DB', DB = None, v = False):
	
	if DB is None:
		DB = np.empty(0, dtype = dtype)
	
	if groups == 'all': # uses shell-style wildcards
		groups = ['[!_]']
		
	pat = '_'.join(groups) + '*' + '/'
	dirs = glob.glob(os.path.join(basedir, pat))
	for d in dirs:
		_, sess = os.path.split(os.path.split(d)[0])
		try:
			fname = os.path.join(d, sess + '_' + suffix + '.npz')
			db = np.load(fname)['DB']
			for db_ in db:
				DB.resize(DB.size + 1)
				DB[-1] = db_
			nunits = db.size
		except IOError:
			nunits = 0
		if v:
			print '%s\t%i' % (fname, nunits)
	return DB
		

def check_noise(fpath):
	
	f = h5py.File(fpath, 'r')

	npips = 6.
	onset = 0.03

	spktimes = f['spktimes'].value
	spktrials = f['spktrials'].value
	stimparams = f['stimID'].value
	spkwaveform = f['spkwaveform'].value
	rast = f['rast'].value
	rrs = stimparams[:, 1]
	urrs = np.unique(rrs)
	urrs = urrs[urrs>2]
	nrrs = urrs.size


	axs = list()
	ix1 = (stimparams[:, 0] == 0)
	
	fig = plt.figure();
	for i in range(nrrs):
		ix2 = (stimparams[:, 1] == urrs[i])
		ix = np.logical_and(ix1, ix2)
		psth = rast[ix, :].sum(0)
		aligned_psth = RR.aligned_psth(psth, urrs[i], npips, onset, 50)
		ax = fig.add_subplot(nrrs, 2, 2*i+1)
		ax.plot(aligned_psth)
		ax.axvline(20, color = 'r', ls = '-')
		if i == 0:
			ax.set_title(fpath)
		axs.append(ax)
	
	misc.sameyaxis(axs)

	ix1 = ix1.nonzero()[0]
	ix_ = [i for (i, y) in enumerate(spktrials) if y in ix1]
	spkwav = spkwaveform[ix_, :]

	x = Spikes.pca(spkwav, npcs = 2)

	ax = fig.add_subplot(2, 2, 2)
	ax.scatter(x[::5, 0], x[::5, 1])
	xlim = np.array([np.min(x[:, 0]), np.max(x[:, 0])])*1.1
	ylim = np.array([np.min(x[:, 1]), np.max(x[:, 1])])*1.1
	ax.set_xlim(xlim); ax.set_ylim(ylim);
	
	ax = fig.add_subplot(2, 2, 4)
	ax.plot(spkwav[::5, :].T)

	plt.show();


def glm(DB):
	
	cfs = DB['cf']
	group = DB['group'] == 'ko'
	rr = DB['rr']
	db = DB['vs_tone']
	predictors = ['Group', 'RR', 'GxR', 'CF']
	X = np.hstack((group[:, np.newaxis], rr[:, np.newaxis], (group*rr)[:, np.newaxis], cfs[:, np.newaxis]))
	y = db
	Q = sm.GLM(y, X)
	Q.fit()
	print 'Data: Vector strength\n%s\t%s\t%s\n%s' % ('Predictor', 'B', 'p-value', '-'*40)
	for predictor, param, pvalue in zip(predictors, Q.results.params, Q.results.pvalues):
		print '%s\t\t%3.3f\t%3.3f' % (predictor, param, pvalue)


def compare_firing_rates(DB):
	
	groups = ['wt', 'ko']
	ngroups = len(groups)
	base_mean = np.zeros(ngroups)
	base_sem = np.zeros(ngroups)
	evoked_mean = np.zeros(ngroups)
	evoked_sem = np.zeros(ngroups)
	
	fig = plt.figure(figsize = (5, 6));
	ax1 = fig.add_subplot(211); ax1.set_title('Baseline firing rate')
	ax2 = fig.add_subplot(212); ax2.set_title('Evoked firing rate')

	dat1 = DB[DB['group'] == 'wt']['base_mean'] * 1000 / 50
	dat2 = DB[DB['group'] == 'ko']['base_mean'] * 1000 / 50
	dat3 = DB[DB['group'] == 'wt']['evoked_mean'] * 1000 / 30
	dat4 = DB[DB['group'] == 'ko']['evoked_mean'] * 1000 / 30
	ax1.boxplot([dat1, dat2])
	ax2.boxplot([dat3, dat4])
	
	ax1.set_xticklabels(groups)
	ax2.set_xticklabels(groups)

	[_, p_base] = st.ttest_ind(dat1, dat2)
	[_, p_ev] = st.ttest_ind(dat3, dat4)
	print p_base
	print p_ev
	ax1.text(1, 0, 'mean = %3.3f' % np.mean(dat1))
	ax1.text(2, 0, 'mean = %3.3f' % np.mean(dat2))
	ax2.text(1, 0, 'mean = %3.3f' % np.mean(dat3))
	ax2.text(2, 0, 'mean = %3.3f' % np.mean(dat4))
	ax2.set_ylabel('Firing rate (spikes / sec)')
	ax2.set_xlabel('Group')




	
	
	
	