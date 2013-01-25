import os
import numpy as np
import matplotlib.pyplot as plt
import world
params = world.params
import load_session

spktimes = f['spktimes'].value
spktrials = f['spktrials'].value
stimparams = f['stimID'].value
rast = f['rast'].value


basedir = '/Volumes/BOB_SAGET/Fmr1_RR'
e = 'ko_20120127'
for x_ in x:
	rel = os.path.split(x_)[1]
	newrel = rel[1:]
	shutil.move(os.path.join(basedir, e, 'fileconversion', newrel), os.path.join(basedir, e, 'fileconversion', rel))
	try:
		shutil.move(os.path.join(basedir, e, 'fileconversion', 'RR' + newrel[2:]), os.path.join(basedir, e, 'fileconversion', '_RR' + newrel[3:]))
	except:
		pass

fig = plt.figure(figsize = (10, 10));
ax1 = fig.add_subplot(211);
ax2 = fig.add_subplot(212);

for e in experiments:
	ax1.cla(); ax2.cla();
	RF.look_at_map(e, onlygood = True, ax = ax1)
	RF.look_at_map(e, onlygood = False, ax = ax2)
	
	fig.savefig(os.path.join(basedir, 'map' + e + '.png'))


rast = f['rast'].value

for e in experiments:
	files = glob.glob(os.path.join(basedir, e, 'fileconversion', 'RR*.h5'))
	print e
	for file_ in files:
		f = h5py.File(file_, 'r')
		stimparams = f['stimID'].value
		f.close()
		print file_
		if np.unique(stimparams[:, 1])[0] == 0:
			print '^^^^^^  AFUUFFFFFCCCCCCCCFFFF1!!!  ^^'



'''---Spike train distances---'''
dist_in = []
dist_il = []
experiment = []
unit = []
ufreqs = []
uattens = []
for e in experiments:
	files = glob.glob(os.path.join(basedir, e, 'fileconversion','RR*.h5'))
	for file_ in files:
		f = h5py.File(file_, 'r')
		dist_in.append(RR.calc_vR_dist_intertrain(f['spktimes'].value, f['spktrials'].value, f['stimID'].value))
		dist_il.append(RR.calc_vR_dist_intertrial(f['spktimes'].value, f['spktrials'].value, f['stimID'].value))
		ufreqs.append(np.unique(f['stimID'].value[:, 0]))
		uattens.append(np.unique(f['stimID'].value[:, 1]))
		f.close()
		experiment.append(e)
		unit.append(file_)


group = [e[:2] for e in experiment]
unitnum = [np.int32(os.path.split(u)[1][2:5]) for u in unit]

# get unit CFs
cf = []
for u, e in zip(unitnum, experiment):
	cfs = np.loadtxt(os.path.join(basedir, e, 'cfs.txt'))
	cf_ = cfs[cfs[:, 0] == u, 1]
	cf.append(cf_[0])


dtype = np.dtype([('rr', 'f4'), ('group', 'S2')])
data = np.empty(0, dtype = dtype)
for g, e, il_, in_, ufreq, cf_ in zip(group, experiment, dist_il, dist_in, ufreqs, cf):
	if il_.shape[1] == 10:
		rrs = np.array([2, 4, 6, 8, 10, 12, 14, 16, 18, 20.])
	elif il_.shape[1] == 6:
		rrs = np.array([1.5, 2, 3, 4, 6, 8])
	elif il_.shape[1] == 7:
		rrs = np.array([1.5, 2, 3, 4, 6, 8, 12])
	elif il_.shape[1] == 8:
		rrs = np.array([1.5, 2, 3, 4, 6, 8, 12, 16])
	[_, ix] = misc.closest(ufreq, cf_, log = True)

	for k, (il_n, il_t, in_n, in_t) in enumerate(zip(il_[0], il_[ix[0]], in_[0], in_[ix[0]])):
		data.resize(data.size + 1)
		data[-1] = np.array((rrs[k], il_n, il_t, in_n, in_t, g), dtype = dtype)


indeps = ('dist_il_noise', 'dist_il_tone', 'dist_in_noise','dist_in_tone')
fig = plt.figure();
for i, indep in enumerate(indeps):
	ax = fig.add_subplot(2, 2, i+1)
	analysis.bar_by_dep_by_group(indep, 'rr', data, bins = np.array([1.5, 2, 3, 4, 6, 8, 12]), ax = ax)





nspks = np.logspace(1, 4, 20)
N = 1000
duration = 2. # in seconds
nbins = duration * 1000.
dist = np.zeros((nspks.size, N))
for i in range(nspks.size):
	for n in range(N):
		spktimes_1 = np.random.random(nspks[i]) * duration
		spktimes_2 = np.random.random(nspks[i]) * duration
		dist[i, n] = Spikes.vR_dist(spktimes_1, spktimes_2, nbins = nbins)
	
dist_mean = np.mean(dist, 1)
dist_sem = st.sem(dist, 1)


fig = plt.figure();
deps = ('power_noise', 'power_tone')
axs = []
for i, dep in enumerate(deps):
	ax = fig.add_subplot(2, 1, i+1)
	analysis.bar_by_dep_by_group('rr', dep, data, bins = np.array([1.5, 2, 3, 4, 6, 8, 12]), iscategorical = True, ax = ax)
	axs.append(ax)

misc.sameyaxis(axs)


# remove the units that don't follow noise at the lowest train rate
d = db.db(data)

q = d.find('rr < 2.1', 'p_noise > 0.05')

badunits = np.unique(zip(q['sess'], q['unit']))

for sess in np.unique(badunits[:, 0]):
	np.savetxt(os.path.join(basedir, sess, 'nonfollowing.txt'), np.int32(badunits[badunits[:, 0] == sess, 1]))

# remove units that don't have evoked firing rate of at least 2 times baseline firing rate
ix = d.data['evoked_mean'] < 2.*d.data['base_mean']
badunits = np.unique(zip(d.data[ix]['sess'], d.data[ix]['unit']))
for sess in np.unique(badunits[:, 0]):
	np.savetxt(os.path.join(basedir, sess, 'lowevoked.txt'), np.int32(badunits[badunits[:, 0] == sess, 1]))





# for each distance measure, plot distance vs CF (RRs done separately)
groups = ['wt', 'ko']
deps = ['dist_il_noise', 'dist_il_tone', 'dist_in_noise', 'dist_in_tone']
rrs = np.array([1.5, 2, 3, 4, 6, 8, 12])
nrrs = rrs.size
for dep in deps:
	axs = []
	fig = plt.figure(figsize = (11, 9));
	fig.suptitle(dep)
	for i, rr in enumerate(rrs):
		ax = fig.add_subplot(4, 2, i+1)
		axs.append(ax)
		data_ = data[data['rr'] == rr]
		data_ = data_[data_['p_noise'] < 0.05]
		analysis.bar_by_dep_by_group('cf', dep, data_, bins = np.arange(0, 41, 10), iscategorical = False, ax = ax)
		
	misc.sameyaxis(axs)
	fig.savefig(os.path.join(basedir, dep + '_by_exp_following.png'))
	

# for each distance measure, plot distance vs RR (CFs done separately)
groups = ['wt', 'ko']
deps = ['dist_il_noise', 'dist_il_tone', 'dist_in_noise', 'dist_in_tone']
cfs = np.arange(0, 41, 10)
ncfs = cfs.size-1
for dep in deps:
	axs = []
	fig = plt.figure(figsize = (11, 9));
	fig.suptitle(dep)
	for i in range(ncfs):
		ax = fig.add_subplot(4, 1, i+1)
		axs.append(ax)
		data_ = data[misc.isbetween(data['cf'], cfs[i], cfs[i+1])]
		data_ = data_[data_['p_noise'] < 0.05]
		analysis.bar_by_dep_by_group('rr', dep, data_, bins = np.array([1.5, 2, 3, 4, 6, 8, 12]), iscategorical = True, ax = ax)

	misc.sameyaxis(axs)
	fig.savefig(os.path.join(basedir, dep + '_by_cf_following.png'))
	


clrs = 'bg'
for i, dep in enumerate(deps):
	fig = plt.figure();
	ax = fig.add_subplot(111);
	for group in groups:
		q = d.find('group == ' + group)
		x = q[dep]
		X = np.vstack((x, np.ones(x.size))).T
		y = q['power_noise']
		Q = boots_lr.boots_lr(X, y, 2000)
		Q.fit()
		Q.plot(ax)



deps = ['power_tone', 'power_noise']
ndeps = len(deps)
fig = plt.figure();
for i, dep in enumerate(deps):
	ax = fig.add_subplot(ndeps, 1, i+1)
	analysis.bar_by_dep_by_group('rr', dep, data, bins = np.array([1.5, 2, 3, 4, 6, 8, 12]), grouptype = 'sess', ax = ax)
	
deps = ['power_noise', 'power_tone']
cfs = np.arange(5, 46, 20)
ncfs = cfs.size-1
for dep in deps:
	axs = []
	fig = plt.figure(figsize = (11, 9));
	fig.suptitle(dep)
	for i in range(ncfs):
		ax = fig.add_subplot(ncfs, 1, i+1)
		axs.append(ax)
		data_ = data[misc.isbetween(data['cf'], cfs[i], cfs[i+1])]
		data_ = data_[data_['p_noise'] < 0.05]
		analysis.bar_by_dep_by_group('rr', dep, data_, bins = np.array([1.5, 2, 3, 4, 6, 8, 12]), iscategorical = True, ax = ax, grouptype = 'sess')

	misc.sameyaxis(axs)
	# fig.savefig(os.path.join(basedir, dep + '_by_cf_following.png'))




e = 'ko_20120306'
files = glob.glob(os.path.join(basedir, e, 'fileconversion', 'RR*.h5'))
nfiles = len(files)
freq = 0
rrtf = np.empty((nfiles, 7))
for i, file_ in enumerate(files):
	f = h5py.File(file_, 'r')
	rast = f['rast'].value
	stimparams = f['stimID'].value
	blockid = f['blockID'].value
	f.close()
	rrs = np.unique(stimparams[:, 1])
	rrtf[i, :] = RR.calc_rrtf_all(rast, stimparams, freq, rrs)


freq = 0
rr = 1.5
onset = 0.05

ix = stimparams[:, 0] == freq
nspks_1st = rast[ix, 50:100].mean()

rast_ = rast[RF.get_trials(stimparams, [freq, rr]), :]

npips = 6.
pip_start = np.int32((onset + np.arange(npips) / rr) * 1000)
pip_end = pip_start + 50
response = np.empty(npips) * np.nan
fig = plt.figure();
for i in np.arange(npips):
	rast_chunk = rast_[:, pip_start[i]:pip_end[i]]
	response[i] = rast_chunk.sum()
	ax = fig.add_subplot(npips, 1, i+1)
	ax.imshow(1-rast_chunk, aspect = 'auto', interpolation = 'nearest', cmap = 'gray')
	
plt.show()
	
	
	
	
for e in experiments:
	
	files = glob.glob(os.path.join(basedir, e, 'fileconversion', 'RR*.h5'))
	f = h5py.File(files[0], 'r')
	stimparams = f['stimID'].value
	f.close()
	print e
	print np.unique(stimparams[:, 1])
	



	
rast = f['rast'].value	
rf = f['rf'].value
stimparams = f['stimID'].value

nattens, nfreqs = rf.shape
ufreqs = np.unique(stimparams[:, 0])
uattens = np.unique(stimparams[:, 1])
rf_psth = np.empty((nfreqs, nattens, rast.shape[1]))
rf_psth_mask = np.empty((nfreqs, nattens))
wind = np.hamming(10); wind = wind / np.sum(wind)
base_psth = np.convolve(rast[:50, :].mean(0), wind, 'same')
baselinestd = base_psth.std()
baselinemean = base_psth.mean()
thresh = baselinemean + 3*baselinestd
for i in range(nfreqs):
	for j in range(nattens):
		ix = RF.get_trials(stimparams, (ufreqs[i], uattens[j]))
		rf_psth_ = rast[ix, :].mean(0)
		rf_psth_mask[i, j] = rf_psth_[55:80].mean() > thresh
		rf_psth[i, j, :] = rf_psth_

B = RF.findmaxcluster(rf_psth_mask)
rf_psth_mask = np.ma.make_mask(B[0])
ev_psth = rf_psth[rf_psth_mask].mean(0)
fig = plt.figure();
ax = fig.add_subplot(111);
ax.plot(ev_psth)
plt.show()




fig = plt.figure();
ax1 = fig.add_subplot(211);
rf = RF.calc_rf(rast, stimparams)
ax1.imshow(rf);
ax2 = fig.add_subplot(212);
lfp_rf = RF.calc_lfp_rf(lfp, stimparams)
ax2.imshow(lfp_rf);
plt.show();


fig = plt.figure();
indeps = ['freq', 'rate']
for i, indep in enumerate(indeps):
	ax = fig.add_subplot(1, 2, i+1)
	analysis.startle_ratio_by_indep(indep, db, ax = ax, show_all = True)

X = db['startle'][:100, :]
nobs, nsamp = X.shape

X_mean = X.mean(1)
X_std = X.std(1)
X -= np.tile(X_mean, (nsamp, 1)).T
X /= np.tile(X_std, (nsamp, 1)).T

X_ = ica.fit(X.T).transform(X.T)

fig = plt.figure();
q = np.linspace(0, 5, 10)
Q = np.tile(q, (X_.shape[0], 1)).T
for i in range(4):
	ax = fig.add_subplot(1, 4, i+1)
	plot(X_[:, (i*10):((i+1)*10)]+Q.T)
	
plt.show()


fig = plt.figure();
indeps = ['freq', 'rate']
for i, indep in enumerate(indeps):
	ax = fig.add_subplot(1, 2, i+1)
	analysis.startle_ratio_by_indep(indep, db, ax = ax, show_all = True)

sesss = np.unique(DB['sess'])
for sess in sesss:
	DB_ = DB[DB['sess']==sess]
	dats = [DB_, DB_[DB_['animal']<4], DB_[DB_['animal']>3]]
	titles = ['All', 'Old', 'Young']
	fig = plt.figure(figsize = (10, 10));
	axs = []
	for i, (dat, title) in enumerate(zip(dats, titles)):
		ax = fig.add_subplot(3, 1, i+1);
		axs.append(ax)
		analysis.startle_ratio(dat, title = title, ax = ax)
		ax.set_xlim([-1, ax.get_xlim()[1]])
		if i!=len(titles):
			ax.set_xlabel('')
		# ax.legend(l, loc = 'lower right')
		# fig.suptitle(sess)
		plt.show()
	misc.sameyaxis(axs)




sesss = np.unique(DB['sess'])
animals = np.unique(DB['animal'])
i = 0
fig = plt.figure();
axs = []
for sess in sesss:
	for a, animal in enumerate(animals):
		i += 1
		ax = fig.add_subplot(sesss.size, animals.size, i)
		axs.append(ax)
		startle = DB[np.c_[DB['sess']==sess, DB['animal']==animal].all(1)]['startle']
		startle = np.abs(startle)
		misc.errorfill(startle, ax = ax)
		if a==0:
			ax.set_ylabel(sess)
misc.sameyaxis(axs)

sesss = ['wt_nai_20120208', 'wt_nai_20120410', \
'wt_exp_20120314', 'wt_exp_20120320', 'wt_exp_20120508', \
'wt_exp_20120509', 'wt_nai_20120412', \
'ko_exp_20120315', 'ko_nai_20120202', 'ko_nai_20120305', \
'ko_exp_20120321', 'ko_nai_20120206', 'ko_nai_20120306', \
'ko_exp_20120326', 'ko_exp_20120420', 'ko_nai_20120209', \
'ko_w2_20120521', 'ko_w2_20120523']

dep_keys = ['rrtf_noise_norm', 'rrtf_tone_norm']
fig = plt.figure();
ax[0] = fig.add_subplot(211)
ax[1] = fig.add_subplot(212)
for sess in sesss[:1]:
	for i, dep_key in enumerate(dep_keys):
		data = DB[DB['sess']==sess]
		analysis.bar_by_indep_by_group_2d(dep_key, 'rr', group_key = 'cf', data, ax = ax[i])
		

'''do it all'''
for sess in sesss:
	
	'''1. fileconversion'''
	fileconversion.fileconversion(sess)
	
	'''2. add_bf_man'''
	RF.add_bf_man(sess)
	
	'''3. look_at_map'''
	RF.look_at_map(sess, makechanges = True)
	
	
	'''4. remove units'''
	RF.remove_no_cf_units(sess)
	RF.remove_units(sess, 'nonA1')

	'''5. characterize'''
	RR_analysis.characterize(sess)
	
	'''6. look at contact sheets'''
	'''7. create nonfollowing.txt'''
	
	'''8. remove nonfollowing units'''
	try:
		RF.remove_units(sess, 'nonfollowing')
		RF.remove_from_DB(sess, 'nonfollowing')
	except:
		pass



'''print maps'''
ncols = np.ceil(np.sqrt(len(sesss)))
fig = plt.figure()
ax = fig.add_subplot(111);
for s, sess in enumerate(sesss):
	# try:
	ax.cla()
	RF.look_at_map(sess, ax = ax)
	ax.set_title(sess)
	fig.savefig(os.path.join(basedir, 'maps', sess + '_map.png'))
	except:
		print sess
	
	
fig = plt.figure();
ax = fig.add_subplot(111);
clrs = 'brgymck'*10
for sess in sesss:

	if sess[:2] == 'ko':
		clr = 'r'
	elif sess[:2] == 'wt':
		clr = 'b'
	if sess[3:6] == 'nai':
		ls = '-'
	elif sess[3:6] == 'exp':
		ls = '--'
	analysis.bar_by_indep_2d('rrtf_noise_norm', 'rr', DB[DB['sess']==sess], ax = ax, c = clr, ls = ls)
	
	
'''rearing effect'''
gens = ['wt', 'ko']
deps = ['rrtf_noise', 'rrtf_tone', 'rrtf_noise_norm', 'rrtf_tone_norm']
ax = []
for gen in gens:
	fig = plt.figure();
	fig.suptitle(gen)
	for i, dep in enumerate(deps):
		ax_ = fig.add_subplot(2, 2, i+1)
		analysis.bar_by_indep_by_group_2d(dep, 'rr', DB[DB['gen']==gen], group_key = 'exp', ax = ax_)
		ax.append(ax_)
misc.vline(ax, [8], color = 'r', ls = '--')
misc.sameyaxis(ax)
plt.show()

'''genotype effect'''
exps = ['nai', 'exp']
deps = ['rrtf_noise', 'rrtf_tone', 'rrtf_noise_norm', 'rrtf_tone_norm']
ax = []
for exp in exps:
	fig = plt.figure();
	fig.suptitle(exp)
	for i, dep in enumerate(deps):
		ax_ = fig.add_subplot(2, 2, i+1)
		analysis.bar_by_indep_by_group_2d(dep, 'rr', DB[DB['exp']==exp], group_key = 'gen', ax = ax_)
		ax.append(ax_)
misc.vline(ax, [8], color = 'r', ls = '--')
misc.sameyaxis(ax)
plt.show()

fig = plt.figure()
ax = fig.add_subplot(111)
ix = np.vstack([DB['exp']=='nai', DB['cf']>30, DB['gen']=='ko']).all(0)
plot_x = np.unique(DB['rr'][ix])
plot_x = plot_x[~np.isnan(plot_x)]
ax.plot(plot_x, DB['rrtf_noise_norm'][ix].mean(0))

ix = np.vstack([DB['exp']=='nai', DB['cf']>30, DB['gen']=='wt']).all(0)
plot_x = np.unique(DB['rr'][ix])
plot_x = plot_x[~np.isnan(plot_x)]
ax.plot(plot_x, DB['rrtf_noise_norm'][ix].mean(0))

plt.show()


rrs = np.array([1.5, 2, 3, 4, 6, 8, 12])

sesss = [['wt_exp_20120314', 'wt_exp_20120320'], ['wt_exp_20120508', 'wt_exp_20120509']]
fig = plt.figure()
ax = fig.add_subplot(111)
for sess in sesss:
	DB = RR_analysis.combine_DB(sess)
	dat_mean = np.mean(DB['rrtf_tone_norm'], 0)
	dat_sem = st.sem(DB['rrtf_tone_norm'], 0) 
	ax.errorbar(rrs, dat_mean, yerr = dat_sem)


DB['rrtf_noise']


dtype = np.dtype([('gen', 'S3'), ('exp', 'S3'), ('sess', 'S20'), ('cf', 'f4'), ('ev_psth', '333f4')])
DB = np.empty(0, dtype = dtype)
for sess in sesss:
	gen, exp, _ = sess.split('_')
	cfs = np.loadtxt(os.path.join(basedir, sess, 'cfs.txt'))
	files = glob.glob(os.path.join(basedir, sess, 'fileconversion', 'RF*.h5'))
	for file_ in files:
		unitnum = np.int32(re.findall('(\d+)', os.path.splitext(os.path.split(file_)[-1])[0])[0])
		cf = cfs[cfs[:, 0] == unitnum, 1][0]
		f = h5py.File(file_, 'r')
		rast = f['rast'].value
		stimparams = f['stimID'].value
		f.close()
		ev_psth = RF.calc_evoked_psth(rast, stimparams)
		DB.resize(DB.size+1)
		DB[-1] = np.array((gen, exp, sess, cf, ev_psth[:333]), dtype = dtype)




fig = plt.figure()
ax = fig.add_subplot(111)
bins = np.arange(0, 50, 10)
colors = 'bgmcry'
cf_bin = misc.bin(DB['cf'], bins)
dat_mean = DB['rrtf_tone_norm'].mean(1)
med = np.median(dat_mean)
dat_bin = dat_mean > med
bins = [True, False]
for i, bin in enumerate(bins):
	db = DB[dat_bin == bin]
	data = db['ev_psth']
	misc.errorfill(data, ax, color = colors[i])
	
	
	
	
	
# check number of reps for each RR
for sess in sesss:
	fnames = glob.glob(os.path.join(basedir, sess, 'fileconversion', 'RR*.h5'))
	f = h5py.File(fnames[0], 'r')
	stimparams = f['stimID'].value
	f.close();
	print '%s\t%i' % (sess, RF.get_trials(stimparams, [0, 2]).size)
	


ix = (DB['thresh'] <= 2).nonzero()[0]
for i in ix:
	print DB['sess'][i], DB['unit'][i]

cf_bins = np.arange(0, 50, 10)
cf_binned = misc.bin(DB['cf'], bins = cf_bins)
for cf_bin in cf_bins:
	cf_binned = cf_bin
	

# check BW and Threshold calculations
# plot RF and thresholded RF side-by-side
fig = plt.figure();
for sess in sesss:
	fig.clf()
	fnames = glob.glob(os.path.join(basedir, sess, 'fileconversion', 'RF*.h5'))
	nfiles = len(fnames)
	ncols = 2
	nrows = np.ceil(nfiles/2.)
	for i, fname in enumerate(fnames):
		ax = fig.add_subplot(nrows, ncols, i+1)
		f = h5py.File(fname, 'r')
		rast = f['rast'].value
		stimparams = f['stimID'].value
		rf = f['rf'].value; rf = rf / rf.max()
		f.close()
		ev_psth, B = RF.calc_evoked_psth(rast, stimparams)
		RF.plot_RF(np.hstack((rf, B.T)), ax = ax)
		ax.set_title(os.path.split(fname)[-1])
	fig.savefig(sess+'_RF_thresh.png')



# make file structure
import glob, os, re

basedir = '/Users/robert/Desktop/Fmr1_RR/Sessions'
fnames = glob.glob(os.path.join(basedir, '*[wt|ko]*'))
for fname in fnames:
	p = [os.path.join(basedir, fname, 'experimentfiles'), \
		os.path.join(basedir, fname, 'fileconversion'), \
		os.path.join(basedir, fname, 'data')]
	for p_ in p:
		if not os.path.exists(p_):
			os.mkdir(p_)

'''
rename all *_R.mat to *_R1.mat
'''
dirs = glob.glob(os.path.join(basedir, '[wt|ko]*'))
for d in dirs:
	files = glob.glob(os.path.join(basedir, d, 'data', '*_R.mat'))
	for f in files:
		f_abs, f_rel = os.path.split(f)
		f_name, f_ext = os.path.splitext(f_rel)
		f_name2 = f_name + '1'
		os.rename(f, os.path.join(f_abs, f_name2 + f_ext))

x = np.empty((0, 3), dtype = 'str')
for d in dirs:
	g, w, _ = d.split('_')
	x = np.vstack((x, np.hstack((d, g, w))))
	


'''
fix the _R problem
(rename all R## < max(R##) to _R## and Rmax## to b/r**.mat)
'''
dirs = glob.glob(os.path.join(basedir, '[wt|ko]*'))
for d in dirs:
	# f_abs, sess = os.path.split(d)
	# try:
	# 	fileconversion.fileconversion(sess)
	# except:
	# 	print sess
	fileconversion.sort_redos(d)



'''
recreate coordinate text files from an image of the map
'''
sess = 'ko_nai_20120113'
sesss = glob.glob('/Users/robert/Desktop/aaa/*.jpg')
sess = os.path.splitext(os.path.split(sesss[9])[1])[0]
print sess
imgdir = '/Users/robert/Desktop/aaa/'

I = misc.loadimg(os.path.join(imgdir, sess + '.jpg'))
fig = plt.figure()
ax = fig.add_subplot(111)
ax.imshow(I)
ax.set_xlim([-20, ax.get_xlim()[1]])
plt.show()
i = 0
coords = np.empty((0, 3), dtype = float32)
while x > 0:
	i += 1
	ax.set_title(i)
	plt.show()
	(x, y) = zip(*plt.ginput(1, timeout = 0))
	coords = np.vstack((coords, np.hstack((i, x, y))))
	if x > 0:
		np.savetxt(os.path.join(imgdir, sess + '.txt'), coords, '%4.4f', delimiter = '\t')


dirs = glob.glob(os.path.join('/Users/robert/Desktop/robert/', '[wt|ko]*'))
for d in dirs:
	f_abs, sess = os.path.split(d)
	dest_dir = os.path.join(basedir, sess)
	if len(glob.glob(os.path.join(dest_dir, 'data', 'b*.mat'))) == 0:
		print os.path.join(d, 'data'), dest_dir
		os.removedirs(os.path.join(dest_dir, 'data'))
		shutil.move(os.path.join(d, 'data'), dest_dir)
	
	

sesss = ['ko_exp_20120420', 'ko_exp_20120420', 'ko_nai_20120113', 'ko_nai_20120116', 'ko_nai_20120127', 'ko_nai_20120202', 'ko_nai_20120206', 'ko_nai_20120209', 'ko_nai_20120229', 'ko_nai_20120305', 'ko_nai_20120306', 'wt_exp_20120314', 'wt_nai_20120123', 'wt_nai_20120125', 'wt_nai_20120126', 'wt_nai_20120127', 'wt_nai_20120203', 'wt_nai_20120208', 'wt_nai_20120223', 'wt_nai_20120225', 'wt_nai_20120226', 'wt_nai_20120227', 'wt_w3_20120703']

basedir = '/Volumes/BOB_SAGET/Fmr1_RR/'
fnames = glob.glob(os.path.join(basedir, 'Sessions', 'wt_w[1|2]*'))


for fname in fnames:
	sess = 
	I = misc.loadimg(os.path.join(sess, 'experimentfiles', ))

# quick look at how many fileconversion files there are for each session
fnames = glob.glob(os.path.join(basedir, '[wt|ko]*'))
for fname in fnames:

	print os.path.splitext(os.path.split(fname)[-1])[0]
	sess = os.path.splitext(os.path.split(fname)[-1])[0]
	nb = len(glob.glob(os.path.join(fname, 'fileconversion', 'RF*.h5')))
	nr = len(glob.glob(os.path.join(fname, 'fileconversion', 'RR*.h5')))
	print '.'*nb
	print '*'*nr
	print

# which sessions have _DB files?
fnames = glob.glob(os.path.join(basedir, 'Sessions', '[wt|ko]*'))
for fname in fnames:

	sess = os.path.splitext(os.path.split(fname)[-1])[0]
	
	print sess
	if os.path.exists(os.path.join(basedir, 'DB', sess + '_DB.npz')):
		print '*'*20
	print
	
basedir = '/Volumes/BOB_SAGET/Fmr1_RR'
	
# which sessions have cfs.txt?
fnames = glob.glob(os.path.join(basedir, 'Sessions', '[wt|ko]*'))
for fname in fnames:

	sess = os.path.splitext(os.path.split(fname)[-1])[0]
	
	print sess
	if os.path.exists(os.path.join(fname, 'cfs.txt')):
		print '*'*20
	print


sesss = ['wt_w1_20120625', 'wt_w1_20120626', 'wt_w1_20120628', 'wt_w1_20120629', 'wt_w2_20120622']
	
	
basedir = '/Volumes/BOB_SAGET/Fmr1_RR/Sessions'
destdir = '/Users/robert/Desktop/Fmr1_RR/Sessions'
dirs = glob.glob(os.path.join(basedir, '*[wt|ko]*'))
for d in dirs:
	f_abs, f_rel = os.path.split(d)
	dest_dir = os.path.join(destdir, f_rel)
	if not os.path.exists(dest_dir):
		os.mkdir(dest_dir)
	shutil.copytree(os.path.join(d, 'experimentfiles'), os.path.join(dest_dir, 'experimentfiles'))



basedir = '/Volumes/BOB_SAGET/Fmr1_RR/'
fnames = glob.glob(os.path.join(basedir, '[wt|ko]*'))
for fname in fnames:
	_, f_rel = os.path.split(fname)
	sess = f_rel.split('_DB.npz')[0]
	shutil.move(fname, os.path.join(basedir, 'Sessions', sess, f_rel))


basedir = '/Volumes/BOB_SAGET/Fmr1_RR/Sessions'
for s in sesss:
	if os.path.exists(os.path.join(basedir, s, 'nonfollowing.txt')):
		RF.remove_from_DB(s, 'nonfollowing')
	if os.path.exists(os.path.join(basedir, s, 'nonA1.txt')):
		RF.remove_from_DB(s, 'nonA1')
		

# rerun all fileconversion, RR_analysis.characterize, remove bad unit
sesss = glob.glob(os.path.join(basedir, '[!_]*'))
for sess in sesss:
	sess = os.path.basename(sess)
	fileconversion.fileconversion(sess, v = True)
	RF.remove_no_cf_units(sess)
	RF.remove_units(sess, 'nonfollowing')
	RF.remove_units(sess, 'nonA1')
	RR_analysis.characterize(sess)
	
''' discover the data type, endianness, and header size of a binary file'''
# endians = '!><@='
endians = '>'
offsets = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15]
# dtypes = 'hHiIlLqQfd'
dtypes = 'f'
sampwidths = [4]
nsamps = 3200000
i = 0
for (e, o, (d, s)) in itertools.product(endians, offsets, zip(dtypes, sampwidths)):
	print e, o, d, s
	f.seek(o)
	xs = f.read(nsamps*s)
	try:
		x = struct.unpack('%s%i%s' % (e, nsamps, d), xs)
		fig = plt.figure(1)
		ax = fig.add_subplot(6, 6, np.mod(i, 36)+1)
		# ax.set_yticks([])
		# ax.set_xticks([])
		ax.plot(x[1000000:1500000])
		ax.set_title('%s%i%s%i' % (e, o, d, s))
		i += 1
	except:
		pass
plt.show()



f = h5py.File('/Users/robert/Desktop/Vocalization/KO4/P13/KO4_11_7_30_2012.h5', 'r')
Pxx = f['Pxx'].value
F = f['F'].value
fs = f['fs'].value
f.close()

F_range = (30000, np.inf)

basedir = '/Users/robert/Desktop/Vocalization'
cage = 'KO4'
sess = 'P13'
fnames = glob.glob(os.path.join(basedir, cage, sess, cage + '_*.h5'))
for fname in fnames[:10]:
	f = h5py.File(fname, 'r')
	Pxx = f['Pxx'].value
	F = f['F'].value
	fs = f['fs'].value
	f.close()
	
	blobs = vocalization.power_search(Pxx, F, F_range)
	figs = vocalization.plot_power_search(blobs, Pxx, F)
	for i, fig in enumerate(figs):
		figname = os.path.join(basedir, 'analysis', os.path.splitext(os.path.split(fname)[1])[0] + '_%2.2i.png' % i)
		fig.savefig(figname)
		plt.close(fig)

fnames = glob.glob(os.path.join(basedir, cage, sess, cage + '_*.h5'))
for fname in fnames[10:20]:
	f = h5py.File(fname, 'r')
	Pxx = f['Pxx'].value
	F = f['F'].value	
	f.close()
	
	F_ix = np.vstack((F>30000, F<np.inf)).all(0)
	Pxx_ = Pxx[F_ix]
	F_ = F[F_ix]
	npts = Pxx.shape[1]
	fs = npts / 60.

	P_thresh = Pxx_.mean(1) + Pxx_.std(1)
	t_width = np.round(0.25*fs)
	t_step = np.round(0.125*fs)

	t = np.arange(0, npts, t_step)
	t_hits = []
	for t_ in t:
		pxx = Pxx_[:, t_:t_+t_width].mean(1)
		if (pxx>P_thresh).any(0):
			t_hits.append(t_)
	
	if len(t_hits)>0:
		fig = plt.figure(figsize = (12, 8))
		for i, t_hit in enumerate(t_hits[:9]):
			ax = fig.add_subplot(3, 3, i+1)
			vocalization.plot_spec(Pxx_[:, t_hit:t_hit+t_width]**0.1, ax = ax)
			# ax.set_xticklabels('')
			ax.set_yticklabels('')
		fig.suptitle(fname)
		fig.tight_layout()
		fig.savefig(os.tempnam(basedir, 'spec_'))
	
	print fname, len(t_hits)
	





# RAT VOCALIZATIONS

basedir = '/Volumes/BOB_SAGET/for_shaowen'
stimnames = np.loadtxt('/Volumes/BOB_SAGET/for_shaowen/stimulus/stimnames.txt', 'S')
for i, stimname in enumerate(stimnames):
	stimnames[i] = os.path.splitext(stimname.split('\\')[-1])[0]
	
nstims = stimnames.size


# make giant contact sheets of psths and rasters
plt.close('all')
voc_fnames = glob.glob('VOC*.h5')
onset_offset = np.loadtxt(os.path.join(basedir, 'sounds', 'onset_offset.txt'))
fig = plt.figure(figsize = (16, 20))
for voc_fname in voc_fnames:
	
	f = h5py.File(voc_fname, 'r')
	rast = f['rast'].value
	stimparams = f['stimID'].value
	f.close()

	rast_by_stim, ustimparams = Spikes.calc_rast_by_stim(rast, stimparams)
	psth, ustimparams = Spikes.calc_psth_by_stim(rast, stimparams)

	ax = []
	for i in range(nstims):
		
		onset = 50 + 1000*onset_offset[i, 0]
		offset = 50 + 1000*onset_offset[i, 1]
		ax_ = fig.add_subplot(10, 4, i+1)
		Spikes.plot_raster(rast_by_stim[:, i, 0, :], ax_)
		ax_.plot(5*Spikes.hamming_smoo(psth[i, 0, :])-5)
		ax_.set_title(stimnames[i])
		ax_.set_xlim([onset-50, offset+50])
		ax_.axvline(50+1000*onset_offset[i, 0], color = 'r', ls = '--')
		ax_.axvline(50+1000*onset_offset[i, 1], color = 'r', ls = '--')
		ax.append(ax_)
		if i>0:
			ax_.set_xticklabels('')
			ax_.set_yticklabels('')

	misc.sameyaxis(ax)
	fig.suptitle(voc_fname)
	# fig.tight_layout()
	plt.show()
	fig.savefig(os.path.join('/Volumes/BOB_SAGET/for_shaowen/analysis', '%s_psth_rast.png' % os.path.splitext(voc_fname)[0]))
	fig.clf();


stim_fnames = glob.glob('/Volumes/BOB_SAGET/for_shaowen/sounds/*.wavtxt')
for stim_fname in stim_fnames:
	x = np.loadtxt(stim_fname)
	h5_fname = os.path.splitext(stim_fname)[0] + '.h5'
	f = h5py.File(h5_fname, 'w')
	f.create_dataset('s', data = x, compression = 'gzip')
	f.create_dataset('fs', data = 200000)
	f.close()
	


# add spectrogram to sound databases
stim_fnames = glob.glob('/Volumes/BOB_SAGET/for_shaowen/sounds/*.h5')
NFFT = 512
noverlap = 256
for stim_fname in stim_fnames:
	f = h5py.File(stim_fname, 'r+')
	P, F, T = specgram(f['s'].value, Fs = f['fs'].value, NFFT = NFFT, noverlap = noverlap)
	f.create_dataset('P', data = P, compression = 'gzip')
	f.create_dataset('F', data = F, compression = 'gzip')
	f.create_dataset('T', data = T, compression = 'gzip')
	f.close()
	
# mark sound onset times
ax, fig = misc.axis_check(None)
# onset_offset = np.zeros((40, 2), np.float32)

for i, stimname in enumerate(stimnames):
	f = h5py.File(os.path.join(basedir, 'sounds', stimname + '.h5'), 'r')
	P = f['P'].value
	F = f['F'].value
	T = f['T'].value
	f.close()
	
	voc.plot_spec(P**0.3, F, T, ax = ax)
	ax.set_title(stimname)
	fig.savefig(os.path.join(basedir, 'sounds', stimname+'.png'))
	ax.cla()
	# x = plt.ginput(2, timeout = -1)
	# onset_offset[i, :] = np.array([x[0][0], x[1][0]])
	

straight_series = np.arange(29, 40)
c_series = np.arange(2, 11)
def calc_tuning_curve(psth, ixs):
	resp = np.zeros(ixs.size)
	for i, ix in enumerate(ixs):
		onset = 58 + 1000*onset_offset[ix, 0]
		# offset = 68 + 1000*onset_offset[ix, 1]
		offset = onset + 30
		resp[i] = psth[ix, 0, onset:offset].mean()
	
	return resp


# plot tuning curves for various stimulus series
ixs = c_series
voc_fnames = glob.glob('VOC*.h5')
# voc_fnames = ['VOC017.h5', 'VOC029.h5']
onset_offset = np.loadtxt(os.path.join(basedir, 'sounds', 'onset_offset.txt'))
plt.close('all')
fig = plt.figure()
ax = []
for i, voc_fname in enumerate(voc_fnames):
	
	print voc_fname
	f = h5py.File(voc_fname, 'r')
	rast = f['rast'].value
	stimparams = f['stimID'].value
	f.close()

	rast_by_stim, ustimparams = Spikes.calc_rast_by_stim(rast, stimparams)
	psth, ustimparams = Spikes.calc_psth_by_stim(rast, stimparams)

	ax_ = fig.add_subplot(4, 2, i+1)
	ax_.plot(calc_tuning_curve(psth, ixs))
	ax_.set_title(voc_fname)
	ax.append(ax_)
	
fig.suptitle('synthesized morph series')
fig.savefig(os.path.join(basedir, 'analysis', 'synth_morph_series.png'))

# make tuning sheet for one neuron
psth_all = np.zeros(0)
for i, ix in enumerate(ixs):
	onset = 58 + 1000*onset_offset[ix, 0]
	# offset = 68 + 1000*onset_offset[ix, 1]
	offset = onset + 30

	psth_ = psth[ix, 0, onset:offset]
	psth_all = np.hstack((psth_all, psth_, np.zeros(100-psth_.size)))


tuning_curve = calc_tuning_curve(psth, ixs)
x = np.arange(15, 15+9*100, 100)

fig = plt.figure(figsize = (10, 8))
ax1 = fig.add_axes([0.1, 0.5, 0.8, 0.4])
ax1.set_title(os.path.splitext(voc_fname)[0])
ax1.plot(x, tuning_curve, 'ro-')
ax1.set_ylabel('Mean firing rate (spikes / second)')
ax1.set_xticks(x)
ax1.set_xticklabels('')
ax1.set_xlim([x[0]-30, x[-1]+30])
ax2 = fig.add_axes([0.1, 0.1, 0.8, 0.4])
ax2.plot(Spikes.hamming_smoo(psth_all, windlen = 10))
ax2.set_xlabel('Stimulus ID')
ax2.set_ylabel('Firing rate (spikes / second)')
ax2.set_ylim([0, 1])
ax2.set_xticks(x)
ax2.set_xticklabels(stimnames[ixs])
ax2.set_xlim([x[0]-30, x[-1]+30])


# plot RFs
plt.close('all');
fig = plt.figure()
rf_fnames = ['RF017.h5', 'RF018.h5', 'RF019.h5', 'RF020.h5']
for i, rf_fname in enumerate(rf_fnames):
	f = h5py.File(rf_fname, 'r')
	rast = f['rast'].value
	stimparams = f['stimID'].value
	f.close()
	rf = RF.calc_rf(rast, stimparams)
	ax = fig.add_subplot(4, 2, i+1)
	RF.plot_RF(rf, ax = ax)
	ax.set_title(rf_fname)



dtype = np.dtype([('gen', 'S2'),  ('exp', 'S3'), ('sess', 'S20'), ('cf', 'f4'), ('rrtf', '7f4')])
db = np.empty(0, dtype = dtype)
for sess in f.iterkeys():
	gen, exp, _ = sess.split('_')
	for unit in f[sess].iterkeys():
		print f[sess][unit]['rrs']
		rrtf = f[sess][unit]['rrtf'].value
		if rrtf.size != 7:
			rrtf = np.append(rrtf, np.nan)
		
		db.resize(db.size+1)
		db[-1] = np.array((gen, exp, sess, f[sess][unit]['cf'].value, rrtf), dtype = dtype)

clrs = {'wt' : {'nai' : 'b-', 'exp' : 'b--'}, 'ko' : {'nai' : 'r-', 'exp' : 'r--'}}
gens = ['wt', 'ko']
exps = ['nai', 'exp']
rrs = np.array([1.5, 2, 3, 4, 6, 8, 12])
fig = plt.figure()
ax = fig.add_subplot(111)
leg = []
for gen, exp in itertools.product(gens, exps):
	clrs_ = clrs[gen][exp]
	data = db[np.vstack((db['gen']==gen, db['exp']==exp)).all(0)]['rrtf']
	ax.errorbar(rrs, data.mean(0), yerr = st.sem(data), color = clrs_[0], linestyle = clrs_[1:])
	leg.append(' '.join((gen, exp)))

ax.set_title('Repetition rate transfer function')
ax.set_xlabel('Repetition rate (pips/sec)')
ax.set_ylabel('RRTF')
ax.legend(leg)
plt.show()

ax = fig.add_subplot(122)
leg = []
for gen, exp in itertools.product(gens, exps):
	clrs_ = clrs[gen][exp]
	data = db[np.vstack((db['gen']==gen, db['exp']==exp, db['cf']>30)).all(0)]['rrtf']
	ax.errorbar(rrs, data.mean(0), yerr = sem(data), color = clrs_[0], linestyle = clrs_[1:])
	leg.append(' '.join((gen, exp)))

ax.set_title('RRTF for units outside of exposure range (4-32 kHz)')
ax.set_xlabel('Repetition rate (pips/sec)')
ax.set_ylabel('RRTF')
ax.legend(leg)
plt.show()

lfp = f['lfp'].value
rast = f['rast'].value
stimparams = f['stimID'].value
f.close()

lfp_stim = calc
fs_lfp = lfp.shape[1] / 4.2
P_lfp = 1./fs_lfp
lfp_t = np.arange(0, 4.2, P_lfp)

fig = plt.figure()
for i in range(psth_stim.shape[0]):
	for j in range(7):
		ax = fig.add_subplot(4, 7, (i-1)*7 + j)
		ax.plot(lfp_t, lfp_stim[i, j, :])
		ax.set_ylabel('%f %f' % (ustimparams[0][i], ustimparams[1][j]))

plt.show()
	
# do statistics by animal
sesss = np.unique(db['sess'])
rrtfs = {'wt_nai' : np.empty((0, 7)), 'wt_exp' : np.empty((0, 7)), \
'ko_nai' : np.empty((0, 7)), 'ko_exp' : np.empty((0, 7)),}
for sess in sesss:
	db_ = db[db['sess'] == sess]
	gen_ = db_[0]['gen']
	exp_ = db_[0]['exp']
	rrtf_ = db_['rrtf'].mean(0)
	key = gen_+'_'+exp_
	rrtfs[key] = np.vstack((rrtfs[key], rrtf_[np.newaxis, ...]))

dtype2 = np.dtype([('gen', 'S2'), ('exp', 'S3'), ('rr', 'f4'), ('rrtf', 'f4')])
db2 = np.empty(0, dtype = dtype2)
rrs = np.array([1.5, 2, 3, 4, 6, 8, 12])
for key in rrtfs:
	gen, exp = key.split('_')
	for rrtf in rrtfs[key]:
		for rr, rrtf_ in zip(rrs, rrtf):
			db2.resize(db2.size+1)
			db2[-1] = np.array((gen, exp, rr, rrtf_), dtype = dtype2)
	

clrs = {'wt' : {'nai' : 'b-', 'exp' : 'b--'}, 'ko' : {'nai' : 'r-', 'exp' : 'r--'}}
gens = ['wt', 'ko']
exps = ['nai', 'exp']

fig = plt.figure()
ax = fig.add_subplot(111)
leg = []
for gen, exp in itertools.product(gens, exps):
	clrs_ = clrs[gen][exp]
	data = rrtfs[gen+'_'+exp]
	ax.errorbar(rrs, data.mean(0), yerr = st.sem(data), color = clrs_[0], linestyle = clrs_[1:])
	leg.append(' '.join((gen, exp)))

ax.set_title('Repetition rate transfer functions (RRTFs)')
ax.set_xlabel('Repetition rate (pips/sec)')
ax.set_ylabel('RRTF')
ax.legend(leg)
path = misc.get_path()
fig.savefig(path, format = 'eps')

plt.show()


clrs = {'wt' : {'nai' : 'b-', 'exp' : 'b--'}, 'ko' : {'nai' : 'r-', 'exp' : 'r--'}}
exp = 'nai'
leg = []
fig = plt.figure()
ax = fig.add_subplot(111)

for gen in gens:
	clrs_ = clrs[gen][exp]
	data = db[db['gen']==gen]['rrtf']
	ax.errorbar(rrs, data.mean(0), yerr = st.sem(data), color = clrs_[0], linestyle = clrs_[1:])
	leg.append(gen)

ax.set_title('Repetition rate transfer functions (RRTFs)')
ax.set_xlabel('Repetition rate (pips/sec)')
ax.set_ylabel('RRTF')
ax.legend(leg)
path = misc.get_path()
fig.savefig(path, format = 'eps')

plt.show()

b1 = 2.
b2 = 3.
b12 = 4.
b3 = 5.

x0 = 10.
x1 = np.random.rand(100)
x2 = np.random.rand(100)
c = np.random.rand(100)

y = x0 + b1*x1 + b2*x2 + b12*x1*x2 + b3*c + np.random.rand(100)
f = open('/Users/robert/Desktop/test.csv', 'w')
f.write(','.join(('y', 'x1', 'x2', 'c')))
f.write('\n')
for i in range(100):
	f.write('%4.4f,%4.4f,%4.4f,%4.4f\n' % (y[i], x1[i], x2[i], c[i]))
f.close()
	

clrs = {'wt' : {'nai' : 'b-', 'exp' : 'b--'}, 'ko' : {'nai' : 'r-', 'exp' : 'r--'}}
gens = ['wt', 'ko']
exps = ['nai', 'exp']
leg = []
t = np.arange(-50, 333-50)
fig = plt.figure()
ax = fig.add_subplot(111)
for gen in gens:
	data = db_[db_['gen']==gen]['ev_psth']
	ax.errorbar(t, data.mean(0), yerr = st.sem(data, 0), color = clrs[gen]['nai'][0], ls = clrs[gen]['nai'][1:])
	leg.append(gen+'-nai')
ax.set_xlim((5, 35))
plt.show()


stimparams = f['stimID'].value
for i in range(16):
	rast = f['site%2.2i' % i]['rast'].value
	fig = plt.figure()
	ax = fig.add_subplot(111)
	Spikes.plot_sorted_raster(rast, stimparams, ax = ax)


dtype = np.dtype([('gen', 'S2'), ('exp', 'S3'), ('animal', 'i4'), ('freq', 'f4'), ('atten', 'f4'), ('peaks', '5f4')])

d2 = np.empty(0, dtype = dtype)
for d_ in d:
	gen = d_['gen']
	exp = d_['exp']
	animal = d_['animal']
	freq = d_['freq']
	for a in range(d_['attens'].size):
		atten = d_['attens'][a]
		peaks = d_['peaks'][a]
		d2.resize(d2.size+1)
		d2[-1] = np.array((gen, exp, animal, freq, atten, peaks), dtype = dtype)


savepath = os.path.join(basedir, fname+'.txt')

assert not os.path.exists(savepath)

f = open(savepath, 'w')
line = ''
dat = data[0]
for key_ in dat.dtype.descr:
	key = key_[0]
	newfield = dat[key]
	if newfield.size==1:
		newfield = [newfield]
	for i, newfield_ in enumerate(newfield):
		header = key
		if len(newfield)>1:
			header = '%s%2.2i' % (header, i+1)
		line += '%s\t' % header
line += '\n'
f.write(line)

for dat in data:
	line = ''
	for key_ in dat.dtype.descr:
		key = key_[0]
		newfield = dat[key]
		if newfield.size==1:
			newfield = [newfield]
		for newfield_ in newfield:
			if not type(newfield_) in [str, np.string_]:
				newfield_ = '%.5e' % newfield_
			line += '%s\t' % newfield_
	line += '\n'
	f.write(line)

f.close()

plot_opts = {'wt' : {'color' : 'b'}, 'ko' : {'color' : 'r'}, 'nai' : {'ls' : '-'}, 'exp' : {'ls' : '--'}}

peak_no = np.arange(1, 6)

gens = np.unique(d['gen'])
exps = np.unique(d['exp'])

fig = plt.figure()
ax = [[]]*2
ax[0] = fig.add_subplot(211)
ax[1] = fig.add_subplot(212)
for g, e in itertools.product(gens, exps):
	d_ = d[np.vstack((d['gen']==g, d['exp']==e)).all(0)]
	for a in range(2):
		dat = d_['peaks'][:, i, :]
		y = dat.mean(0)
		yerr = st.sem(dat, 0)
		ax[a].errorbar(peak_no, y, yerr = yerr, color = plot_opts[g]['color'], ls = plot_opts[e]['ls'])

plt.show()


ons = []
for db_ in db:
	on, off = Spikes.calc_on_off(db_['ev_psth'])
	ons.append(on)

basedir = '/Volumes/BOB_SAGET/fmr1_forRobert/naive/'
sess_fnames = glob.glob(os.path.join(basedir, '*'))
for sess_fname in sess_fnames:
	vor_fname = glob.glob(os.path.join(sess_fname, 'VOR2*.txt'))[0]
	vor = np.loadtxt(vor_fname)
	vor = vor[vor[:,0]>0]
	cf = vor[:, :3]
	np.savetxt(os.path.join(sess_fname, 'cf.txt'), cf)


f = h5py.File('/Users/robert/Documents/Work/Bao/Fmr1_Heesoo/RFs.h5')
basedir = '/Volumes/BOB_SAGET/fmr1_forRobert/naive/'
sess_fnames = glob.glob(os.path.join(basedir, '*'))
for sess_fname in sess_fnames:
	print sess_fname
	absol, sess = os.path.split(sess_fname)
	s = f.create_group(sess)
	pen_fnames = glob.glob(os.path.join(sess_fname, 'fileconversion', 'RF*.h5'))
	for pen_fname in pen_fnames:
		print pen_fname
		absol, pen = os.path.split(pen_fname)
		pen = os.path.splitext(pen)[0]
		frf_ = h5py.File(pen_fname, 'r')
		stimparams = frf_['stimID']
		# print stimparams.shape
		frf = frf_['site00']
		funit = s.create_group(pen)
		RF.add_rf_analysis(frf, funit, stimparams)
		frf_.close()

f.close()


fig = plt.figure()
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)
basedir = '/Volumes/BOB_SAGET/fmr1_forRobert/naive/'
for s_ in f:
	s = f[s_]
	for u_ in s:
		u = s[u_]
		RF.plot_RF(u['rf'], ax = ax1)
		RF.plot_RF(u['rf_clust'], ax = ax2)
		title = '_'.join((s_, u_))
		ax1.set_title(title)
		fig.savefig(os.path.join(basedir, 'analysis', title + '.png'), format = 'png', dpi = fig.dpi)
		ax1.cla()
		ax2.cla()
				

fig = plt.figure(figsize = (16, 4))
ax = fig.add_subplot(111)
t = np.arange(333)
dat = db[db['gen']=='wt']['ev_psth']
misc.errorfill(dat, color = 'b', ax = ax)
dat = db[db['gen']=='ko']['ev_psth']
misc.errorfill(dat, color = 'r', ax = ax)

y = db['ev_psth'][:, 60:100]
y = y / np.tile(y.max(1), (40, 1)).T

pca = PCA(n_components = 2)
y_ = pca.fit_transform(y)

km = KMeans(k = 2)
km.fit(y_)

fig = plt.figure()
ax = fig.add_subplot(111)
t = np.arange(333)
dat = db[km.labels_==0]['ev_psth']
misc.errorfill(dat, color = 'b', ax = ax)
dat = db[km.labels_==1]['ev_psth']
misc.errorfill(dat, color = 'r', ax = ax)

peak_time = np.argmax(db['ev_psth'], 1)

fig = plt.figure()
ax = fig.add_subplot(111)
ix = np.argsort(db['resp_on'])
for i in ix[:20]:
	if db[i]['gen']=='wt':
		color = 'b'
	if db[i]['gen']=='ko':
		color = 'r'
	y = db[i]['ev_psth']
	y = y / y.max()
	wind = np.hanning(10)
	wind = wind / wind.sum()
	y = np.convolve(y, wind, 'same')
	ax.plot(y, color = color)


for i in ix[:3]:
	y = db_ko[i]['ev_psth']
	y = y / y.max()
	wind = np.hanning(10)
	wind = wind / wind.sum()
	y = np.convolve(y, wind, 'same')
	plot(y)



fig = plt.figure()
ax = fig.add_subplot(111)
for db_ in db:
	RF.plot_RF(db_['rf_clust'], bw_lr = db_['bw_lr'], ax = ax)
	title = '_'.join((db_['sess'], db_['unit']))
	ax.set_title(title)
	fig.savefig(os.path.join('/Users/robert/Documents/Work/Bao/Fmr1_RR/analysis/RFs', title+'.png'))
	ax.cla()
	




fig = plt.figure()
ax = fig.add_subplot(111)
for s_ in f:
	s = f[s_]
	for u_ in s:
		u = s[u_]
		RF.plot_RF(u['rf_clust'].value, bw_lr = u['bw_lr'].value, ax = ax)
		ax.plot(u['cf'].value, u['thresh'].value-2, '*', color = 'r', ms = 10)
		xy = plt.ginput()
		bw.append((s_, u_, xy))


bw = []
fig = plt.figure()
ax = fig.add_subplot(111)
for (s_, u_) in su:
		u = f[s_][u_]
		RF.plot_RF(u['rf_clust'].value, bw_lr = u['bw_lr'].value, ax = ax)
		ax.plot(u['cf'].value, u['thresh'].value-2, '*', color = 'r', ms = 10)
		plt.show()
		xy_ = plt.ginput(2)
		xy = (xy_[0][0], xy_[1][0])
		bw.append((s_, u_, xy))
		ax.cla()


for (s_, u_, bw_lr) in bw:
	u = f[s_][u_]
	u['bw'][2] = np.abs(bw_lr[1]-bw_lr[0])


















	

