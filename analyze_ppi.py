import os, glob, shutil
import numpy as np
from scipy import stats as st
import pandas as pd
import Gap; reload(Gap)
import matplotlib as mpl
from matplotlib import pyplot as plt
import misc

color_cycle = mpl.rcParams['axes.color_cycle']
colors = dict(control = 'k', salicylate = 'r', prenihl = 'k', postnihl = 'g', postnihl5d = 'y', thalid = 'r', vehicle = 'b', \
	preinjection = 'k', pre = 'k', tnfa = 'r', thalidwashout = 'y', vehiclewashout = 'm')
lss = dict(control = '-', salicylate = '--', prenihl = '-', postnihl = '--', postnihl5d = '--', thalid = '--', vehicle = '--', \
	preinjection = '-', pre = '-', tnfa = '--', thalidwashout = '--', vehiclewashout = '--')

global studydir, animalIDs
basedir = '/Volumes/BOB_SAGET/TNFalpha/tinnitus/behavior'

def compare_condition_diffs(df = None, control = 'preinjection'):

	df = load_experiment()
	df_control = df[df.condition==control]
	animalgp_control = df_control.groupby(('animalID', 'condition'))
	animalmeans_control = animalgp_control.agg(dict(gapratio = np.mean))

	df_manip = df[df.postdate1>0]
	postdategp = df_manip.groupby(('animalID', 'postdate1'))
	postdate_manip = postdategp.agg(dict(gapratio = np.mean, condition = misc.get_first, postdate1 = misc.get_first))
	
	# add a column for gap detection difference for each animal/postdate
	gapdiff = []
	for ((animalID, _), (postdate1, condition, gapratio)) in postdate_manip.iterrows():
		gapdiff.append(gapratio-animalmeans_control.ix[animalID]['gapratio'][0])
	postdate_manip['gapdiff'] = gapdiff

	datecondgp = postdate_manip.groupby(('postdate1', 'condition'))
	datecondmean = datecondgp.agg(dict(gapdiff = np.mean)).unstack(level = 'condition')
	datecondsem = datecondgp.agg(dict(gapdiff = st.sem)).unstack(level = 'condition')

	fig = plt.figure();
	ax = fig.add_subplot(111);
	x = np.arange(len(datecondmean))
	xticklabels = datecondmean.index
	colors = dict(tnfa = 'gray', vehicle = 'white')
	for i, (((_, condition), y), (k, yerr)) in enumerate(zip(datecondmean.iteritems(), datecondsem.iteritems())):
		ax.bar(x+(-i*0.4), y, yerr = yerr, width = 0.4, color = colors[condition], ecolor = 'k', label = condition)
	ax.legend(loc = 'best')
	ax.set_xticks(x)
	ax.set_xticklabels(xticklabels)
	ax.set_xlabel('Days post injection')
	ax.set_ylabel('Change in gap detection ratio')
	ax.set_xlim([-0.5, len(x)-0.5])
	ax.axhline(0, color = 'k')
	ax.set_title(os.path.split(studydir)[1])

	fig.savefig(os.path.join(studydir, 'Analysis', 'ppi_compare_condition_diffs.png'))

def compare_conditions_pairwise3(control = 'preinjection'):
	'''
	Performs a pairwise comparison between a control condition and several "manipulated"
	conditions where the manipulations were performed on DIFFERENT ANIMALS.
	
	This script plots the paired pre-post difference for each animal, and the manipulation
	group is indicated by different marker type/line styles
	'''
	df = load_experiment()
	df = df[np.vstack((5000<df.freq, df.freq<40000)).all(0)]
	ucond = np.unique(df.condition).values
	ncond = len(ucond)
	# df = pd.concat([df[df.condition==control], df[df.postdate1>6]])
	animalgp = df.groupby(('animalID', 'condition'))
	animalmeans = animalgp.agg(dict(gapratio = np.mean))
	condgp = animalmeans.groupby(level = 'condition')
	condmeans = condgp.gapratio.apply(np.mean)
	condsems = condgp.gapratio.apply(st.sem)
	condmean_control = condmeans[control]
	condsem_control = condsems[control]
	condsem_manip = condsems.select(lambda x: x != control)
	condmean_manip = condmeans.select(lambda x: x != control)
	condsem_manip = condsems.select(lambda x: x != control)
	animalIDs = zip(*animalmeans.index)[0]

	styles = dict(tnfa = dict(ls = '-', marker = '^', mfc = [0.2, 0.2, 0.2]), \
		vehicle = dict(ls = '--', marker = 'o', mfc = 'w'))
	fig = plt.figure(figsize = (4, 8)); ax = fig.add_subplot(111)
	for animalID in animalIDs:

		animalmean = animalmeans.ix[animalID]
		comparemean = animalmean.select(lambda x: x != control)
		compare = comparemean.index[0]

		ycontrol = animalmean.ix[control].values[0]
		ycompare = comparemean.ix[compare].values[0]

		ax.plot([0.3, 0.7], [ycontrol, ycompare], '.-', color = 'k', \
			linestyle = styles[compare]['ls'], marker = styles[compare]['marker'], \
			mfc = styles[compare]['mfc'], markersize = 10)
	
	barwidth = 0.2
	ax.bar(0, condmean_control, yerr = condsem_control, width = barwidth, color = 'w', ecolor = 'k')
	for j, ((i, y), (i, yerr)) in enumerate(zip(condmean_manip.iteritems(), condsem_manip.iteritems())):
		ax.bar(0.8+j*barwidth, y, yerr = yerr, width = barwidth, color = 'w', ecolor = 'k', facecolor = styles[i]['mfc'], label = i)

	ax.set_xticks([]); #ax.set_xticklabels([control[:5], 'manip'])
	ax.set_ylabel('Gap ratio')
	ax.set_xlim([-0.1, 1.3]); ax.set_ylim([0, 1])
	ax.legend()

	fig.subplots_adjust(left = 0.2)
	figpath = os.path.join(studydir, 'Analysis', 'ppi_compare_conditions_pairwise3_%s.png' % control)
	fig.savefig(figpath)

def compare_conditions_pairwise2(control = 'preinjection'):

	df = load_experiment()
	# df = pd.concat([df[df.condition==control], df[df.postdate1>6]])
	df = df[np.vstack((5000<df.freq, df.freq<40000)).all(0)]

	uconds = np.unique(df.condition).values

	animalgp = df.groupby(('animalID', 'condition'))
	animalmeans = animalgp.agg(dict(gapratio = np.mean))
	gapratio = animalmeans.unstack('condition').gapratio

	for ucond in uconds:

		gapratio_ = gapratio[~np.isnan(gapratio[ucond])]
		
		groupmeans = gapratio_.apply(np.mean)
		groupsems = gapratio_.apply(st.sem)

		fig = plt.figure(figsize = (3.5, 6))
		ax = fig.add_subplot(111)

		ax.bar([0, 0.8], [groupmeans[control], groupmeans[ucond]], yerr = [groupsems[control], groupsems[ucond]], width = 0.2, color = 'w', ecolor = 'k')
		for i, g in gapratio_.iterrows():
			ax.plot([0.3, 0.7], [g[control], g[ucond]], '.-', color = 'k', markersize = 15)
		
		ax.set_xlim([-0.1, 1.1])
		ax.set_xticks([0.1, 0.9])
		ax.set_xticklabels([control, ucond])
		ax.set_xlabel('Condition')

		ax.set_ylim([0, 1])
		ax.set_ylabel('Gap ratio')

		fig.subplots_adjust(left = 0.2)
		figpath = os.path.join(studydir, 'Analysis', 'ppi_compare_conditions_%s_vs_%s.png' % (control, ucond))
		fig.savefig(figpath)


def compare_conditions_by_postdate1():

	df = load_experiment()
	df_pre = df[df.postdate1<0]
	df = df[df.postdate1>0]
	df['postdate1'][df.postdate1>7] = 7

	freqgp = df_pre.groupby('freq')
	y_pre = freqgp.gapratio.apply(np.mean)
	yerr_pre = freqgp.gapratio.apply(st.sem)


	udate = np.unique(df.postdate1).values
	ndates = len(udate)
	ncols = np.ceil(np.sqrt(ndates))

	sessgp = df.groupby(('postdate1', 'condition', 'freq'))
	# x = 
	y = sessgp.gapratio.apply(np.mean).unstack('condition')
	yerr = sessgp.gapratio.apply(st.sem).unstack('condition')

	colors = dict(tnfa = 'r', vehicle = 'b')
	fig = plt.figure(figsize = (10, 8.5))

	for j, date in enumerate(udate):

		ax = fig.add_subplot(1, 3, j+1)

		y_ = y.ix[date]
		yerr_ = yerr.ix[date]

		x = range(len(y_.index))

		ax.errorbar(x, y_pre, yerr = yerr_pre, color = 'k', ecolor = 'k', label = 'preinjection')

		for ((i, yc), (i, yerrc)) in zip(y_.iteritems(), yerr_.iteritems()):
			ax.errorbar(x, yc, yerr = yerrc, color = colors[i], label = i)
			ax.set_title('%u days post injection' % date)
			ax.set_xticklabels((y_.index / 1000.).astype(int))
			# ax.set_xlim([-1, len(x)])
			format_axis(ax)

	ax.legend()
	figpath = os.path.join(studydir, 'Analysis', 'ppi_compare_postdate1.png')
	fig.savefig(figpath)

def timed_effect(animalID, startdate):

	if type(startdate) is str:
		startdate = misc.str2date(startdate)

	dfs = []
	fpaths = glob.glob(os.path.join(studydir, 'pairpulse', '%s*nihl*.csv' % animalID))
	for fpath in fpaths:
		df_ = pd.read_csv(fpath)
		df.append(df_)

	df = pd.concat(dfs)

	ufreqs = np.unique(df.freq).values
	x = range(len(ufreqs))
	sessgp = df.groupby('sess')

	groups = sessgp['sess'].groups.keys()
	groups.sort()


	fig = plt.figure();
	ax = fig.add_subplot(111);
	for group in groups:
		date = misc.str2date(group, format = 'YYMMDD')
		g = sessgp.get_group(group)
		ax.plot(x, g.gapratio, marker = '.', ls = lss[g.condition[0]], label = '%ud' % (date-startdate).days)

	ax.legend()
	ax.set_xticks(x)
	ax.set_xticklabels(ufreqs)
	ax.set_ylabel('Gap ratio')
	ax.set_title('NIHL effect: %s' % animalID)

	figpath = os.path.join(studydir, 'Analysis', 'ppi_hearinglesion_%s.png' % animalID)
	fig.savefig(figpath)
	plt.close(fig)

def compare_conditions_by_day(conditions = ('tnfa', 'vehicle')):

	nconditions = len(conditions)
	fig = plt.figure(figsize = (10, 8))

	for k, condition in enumerate(conditions):
		df = load_experiment(condition)

		ufreqs = np.unique(df.freq).values
		x = range(len(ufreqs))
		sessgp = df.groupby(('condition', 'sess', 'freq'))
		y = sessgp.agg(dict(gapratio = np.mean))
		yerr = sessgp.agg(dict(gapratio = st.sem))

		ax = fig.add_subplot(1, nconditions, k+1);
		for j, (((i, y_), (i, yerr_))) in enumerate(zip(y.unstack().iterrows(), yerr.unstack().iterrows())):
			misc.errorfill(x, y_, yerr_, label = i, ax = ax, marker = 'o', color = color_cycle[j])

		ax.legend()
		ax.set_xticks(x)
		ax.set_xlim([0, len(x)])
		ax.set_xticklabels(np.int32(ufreqs/1000))
		ax.set_xlabel('Frequency (kHz)')
		format_axis(ax)

	figpath = os.path.join(studydir, 'Analysis', 'ppi_compare_conditions_by_day_%s.png' % '_'.join(conditions))
	fig.savefig(figpath)


def compare_conditions(df = None, conditions = None, ax = None):

	if df is None:
		df = load_experiment(conditions = conditions)

	# df_pre = pd.concat((df[df.condition=='postnihl'], df[df.condition=='prenihl']))
	# ix = ['1' in i for i in df.condition]
	# df_ = df[ix]
	# df = pd.concat((df_pre, df_))
	df['condition'] = [filter(lambda x: x.isalpha(), i) for i in df.condition]
	animalgp = df.groupby(('animalID', 'condition', 'freq'))
	animaldf = animalgp.agg(dict(gapratio = np.mean))
	condgp = animaldf.groupby(level = ('condition', 'freq'))

	ufreqs = np.unique(df.freq).values
	x = np.arange(ufreqs.size)
	y = condgp.agg(np.mean)
	yerr = condgp.agg(st.sem)

	fig = plt.figure(figsize = (10, 8))
	ax = fig.add_subplot(111)
	for ((i, y_), (i, yerr_)) in zip(y.unstack().iterrows(), yerr.unstack().iterrows()):
		# ax.errorbar(x, y_, yerr_, label = i, marker = 'o', color = colors[i], ls = lss[i])
		misc.errorfill(x, y_, yerr_, label = i, ax = ax, marker = 'o', color = colors[i], ls = lss[i])
	
	ax.set_xticks(x)
	ax.set_xlim([0, len(x)])
	ax.set_xticklabels(np.int32(ufreqs/1000))
	ax.set_xlabel('Frequency (kHz)')
	ax.set_ylabel('Gap ratio')
	ax.set_ylim([0, 1.5])
	ax.axhline([1.0], color = 'r', ls = '--')
	ax.legend()
	ax.set_title('Start age: %u' % df.age.min())

	figdir = os.path.join(studydir, 'Analysis')
	if not os.path.exists(figdir):
		os.mkdir(figdir)
	figpath = os.path.join(figdir, 'ppi_compare_conditions.png')
	fig.savefig(figpath)


def single_subject_conditionmeans():

	fig = plt.figure(figsize = (10, 8));
	ax = fig.add_subplot(111);
	figdir = os.path.join(studydir, 'Analysis')
	if not os.path.exists(figdir):
		os.mkdir(figdir)
	for animalID in animalIDs:

		fpaths = glob.glob(os.path.join(studydir, 'pairpulse', '%s*.csv' % animalID))
		dfs = []
		for fpath in fpaths:
			df_ = pd.read_csv(fpath)
			dfs.append(df_)
		df = pd.concat(dfs)

		condgp = df.groupby(('condition', 'freq'))

		ufreqs = np.unique(df.freq).values
		x = np.arange(ufreqs.size)
		y = condgp.agg(dict(gapratio = np.mean))
		yerr = condgp.agg(dict(gapratio = st.sem))

		for ((i, y_), (i, yerr_)) in zip(y.unstack().iterrows(), yerr.unstack().iterrows()):
			# ax.errorbar(x, y_, yerr_, label = i, marker = 'o', color = colors[i], ls = lss[i])
			misc.errorfill(x, y_, yerr_, label = i, ax = ax, marker = 'o', color = colors[i], ls = lss[i])

		ax.set_xticklabels(np.int32(ufreqs/1000))
		ax.set_xlabel('Frequency (kHz)')
		ax.set_ylabel('Gap ratio')
		ax.set_ylim([0, 1.5])
		ax.axhline([1.0], color = 'r', ls = '--')
		ax.legend()

		figpath = os.path.join(figdir, 'conditionmeans_%s.png' % animalID)
		fig.savefig(figpath)
		ax.cla();

def single_subject_dailyresults(condition = 'all', cond_color = True):

	fig = plt.figure(figsize = (10, 8));
	ax = fig.add_subplot(111);

	if condition=='all':
		cond_patt = ''
	else:
		cond_patt = '_'+condition+'_'
	
	for animalID in animalIDs:
		print animalID
		fpaths = glob.glob(os.path.join(studydir, 'pairpulse', '%s*%s*.csv' % (animalID, cond_patt)))
		dfs = []
		for fpath in fpaths:
			df_ = pd.read_csv(fpath)
			dfs.append(df_)
		try:
			df = pd.concat(dfs)
		except:
			continue
		df = df[df.freq<40000]

		ufreqs = np.unique(df.freq).values
		x = range(len(ufreqs))
		sessgp = df.groupby('sess')

		for j, (i, g) in enumerate(sessgp):
			if cond_color:
				color = colors[g.condition[0]]
			else:
				color = np.array(color_cycle)[(j%len(color_cycle))]
			ax.plot(x, g.gapratio, marker = '.', color = color, label = i)

		ax.set_title('%s: %s' % (animalID, condition))
		ax.set_xticks(x)
		ax.set_xticklabels(ufreqs)
		ax.legend()
		format_axis(ax)

		figpath = os.path.join(studydir, 'Analysis', 'ppi_dailyresults_%s_%s.png' % (condition, animalID))
		fig.savefig(figpath)
		ax.cla();

def single_subject_startleampl():

	fig = plt.figure()
	animalIDs = pd.read_csv(os.path.join(studydir, 'dobs.csv'))['animalID'].values
	for animalID in animalIDs:
		fnames = glob.glob(os.path.join(studydir, 'data', 'PPI', Gap.reverse_animalID(animalID), '*.txt'))
		nfnames = len(fnames)
		nrows, ncols = misc.get_subplot_grid(nfnames)
		for i, fname in enumerate(fnames):
			ax = fig.add_subplot(nrows, ncols, i+1)
			df = Gap.txt2pd(fname)
			Gap.plot_startleampl(df, ax)
			ax.set_ylim([0, 8])
			ax.set_title(os.path.split(fname)[1])
		figpath = os.path.join(studydir, 'Analysis', 'ppi_startleampl_%s.png' % animalID)
		fig.savefig(figpath)
		fig.clf();


def load_experiment(conditions = None):

	if conditions is None:
		cond_patts = ['']
	else:
		if type(conditions) is str:
			conditions = [conditions]

		cond_patts = []
		for condition in conditions:	
			cond_patts.append('_'+condition+'_')

	fnames = []
	for cond_patt in cond_patts:
		patt = os.path.join(studydir, 'pairpulse', '[A-Za-z]*%s*' % cond_patt)
		fnames.extend(glob.glob(patt))
	
	dfs = []
	if len(fnames)==0:
		print 'No files found with pattern\n%s' % patt
	else:
		for fname in fnames:
			df_ = pd.read_csv(fname)
			dfs.append(df_)
		df = pd.concat(dfs)
		return df

def format_axis(ax):
	ax.set_xlabel('Frequency (kHz)')
	ax.set_ylabel('Gap ratio')
	ax.set_ylim([0, 1.5])
	ax.axhline([1.0], color = 'r', ls = '--')
	# ax.legend()

def select_study(studyID, gen):
	global studydir, animalIDs
	studydir = os.path.join(basedir, studyID, gen)
	animalIDs = pd.read_csv(os.path.join(studydir, 'dobs.csv'))['animalID'] 
	return studydir

def analyze():
	single_subject_dailyresults(cond_color = False)
	single_subject_conditionmeans()
	single_subject_startleampl()
	try:
		compare_conditions()
		compare_conditions_by_day()
	except:
		pass
