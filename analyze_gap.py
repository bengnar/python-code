import os, glob, shutil
import numpy as np
from scipy import stats as st
import pandas as pd
import Gap; reload(Gap)
import matplotlib as mpl
from matplotlib import pyplot as plt
import misc
import pdb

color_cycle = mpl.rcParams['axes.color_cycle']
colors = dict(control = 'k', salicylate = 'r', prenihl = 'k', postnihl = 'g', postnihl5d = 'y', thalid = 'r', vehicle = 'b', \
	preinjection = 'k', preinjectiontnfa = 'k', preinjectionvehicle = 'k', pre = 'k', tnfa = 'r', thalidwashout = 'y', vehiclewashout = 'm')
lss = dict(control = '-', salicylate = '--', prenihl = '-', postnihl = '--', postnihl5d = '--', thalid = '--', vehicle = '--', \
	preinjection = '-', preinjectiontnfa = '-', preinjectionvehicle = '-', pre = '-', tnfa = '--', thalidwashout = '--', vehiclewashout = '--')



class GapAnalysis(object):

	def __init__(self, studyID=None, cageID=None):

		# self.basedir = '/Volumes/BOB_SAGET/TNFalpha/tinnitus/behavior'
		self.basedir = '/Users/robert/Desktop'
		if (studyID is not None) and (cageID is not None):
			self.select_study(studyID, cageID)

			
	def load_experiment(conditions = None, onlygood = True):

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
			if onlygood:
				patt = os.path.join(studydir, 'gapdetection', '[A-Za-z]*%s*' % cond_patt)
			else:
				patt = os.path.join(studydir, 'gapdetection', '*%s*' % cond_patt)
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


	def select_study(self, studyID, cageID):

		studydir = os.path.join(self.basedir, studyID, cageID)
		animalIDs = pd.read_csv(os.path.join(studydir, 'dobs.csv'))['animalID'] 
		animalIDs = [i for i in animalIDs if i[0]!='_']
		
		self.animalIDs = animalIDs
		self.studydir = studydir
		self.studyID = studyID
		self.cageID = cageID


def compare_condition_diffs_by_freq(df = None):
	if df is None:
		df = load_experiment()

	df = df[df.freq<40000]
	# select the final day's results as the control
	df_pre = df[df.postdate1<0]
	df_pre.groupby('animalID')
	df_control = []
	for animalID, df_ in df_pre.groupby('animalID'):
		last_sess = df_.postdate1.max()
		df_control.append(df_[df_.postdate1==last_sess])
	df_control = pd.concat(df_control)
	df_control.pivot('freq', 'animalID', 'gapratio').plot()

	df_control2 = df_control.pivot('freq', 'animalID', 'gapratio')

	df_manip = df[df.postdate1>0]
	df_manip.pivot('freq', 'animalID', 'gapratio')
	for (animalID, postdate), df_ in df_manip.groupby(('animalID', 'postdate1')):
		
		df_.gapratio-df_control2[animalID]
		pass

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

	fig.savefig(os.path.join(studydir, 'Analysis', 'compare_condition_diffs.png'))

def compare_conditions_pairwise_by_freq(control = 'prenihl'):

	df = load_experiment(onlygood=True)
	df = df[df.freq<40000]
	postdate = df.filter(regex='post*').values
	df_control = df[postdate<0]
	freqgp_control = df_control.groupby(('animalID', 'freq'))
	freqmean_control = freqgp_control.agg({'gapratio': np.mean})

	df_manip = df[postdate>0]
	freqgp_manip = df_manip.groupby(('animalID', 'condition', 'freq'))
	freqmean_manip = freqgp_manip.agg({'gapratio': np.mean})
	
	gapdiff = []
	for (animalID, condition, freq), df_ in freqmean_manip.iterrows():
		gapdiff.append((df_['gapratio'] - freqmean_control.ix[(animalID, freq)]).values[0])
	freqmean_manip['gapdiff'] = gapdiff

	ufreqs = np.unique(df.freq).values
	x = np.arange(len(ufreqs))
	condgp = freqmean_manip.unstack('freq').groupby(level='condition')
	condmean = condgp.agg()
	for cond, df_ in condgp:
		plot(x, df_['gapdiff'])

def compare_conditions_pairwise(df = None, control = 'preinjection'):
	'''
	Performs a pairwise comparison between a control condition and several "manipulated"
	conditions where the manipulations were performed on DIFFERENT ANIMALS.
	
	This script plots the paired pre-post difference for each animal, and the manipulation
	group is indicated by different marker type/line styles
	'''
	if df is None:
		df = load_experiment()
	# remove below 5k and 40k
	df = df[np.vstack((5000<df.freq, df.freq<40000)).all(0)]
	ucond = np.unique(df.condition).values
	ncond = len(ucond)

	'''
	ANIMAL-WISE
	First we'll make groups, one for each animal/condition
	Red123-preinjection, Red123-postinjection
	Blue123-preinjection, Blue123-postinjection, etc...

	Take animal-wise means for each condition.

	CONDITION-WISE (COMBINE ANIMALS)
	'''
	# make animal/condition groups
	animalgp = df.groupby(('animalID', 'condition'))
	# take means of gap performance for those groups
	animalmeans = animalgp.agg(dict(gapratio = np.mean))

	# make condition-wise groups
	condgp = animalmeans.groupby(level = 'condition')
	condmeans = condgp.gapratio.apply(np.mean)
	condsems = condgp.gapratio.apply(st.sem)
	# since we'll compare each of the groups to the "control" group
	# let's separate control from manipulations
	condmean_control = condmeans[control]
	condsem_control = condsems[control]
	condmean_manip = condmeans.select(lambda x: x != control)
	condsem_manip = condsems.select(lambda x: x != control)
	
	animalIDs = zip(*animalmeans.index)[0]

	styles = dict(tnfa=dict(ls='-', marker='^', mfc=[0.2, 0.2, 0.2]),
		vehicle=dict(ls='--', marker='o', mfc='w'),
		prenihl={'ls': '-', 'marker': '.', 'mfc': 'lightgray'},
		postnihl={'ls': '--', 'marker': '.', 'mfc': 'darkgray'},
		thalid1={'ls': '-', 'marker': '.', 'mfc': 'gray'},
		thalid2={'ls': '-', 'marker': '.', 'mfc': 'gray'},
		vehicle1={'ls': '-', 'marker': '.', 'mfc': 'gray'},
		vehicle2={'ls': '-', 'marker': '.', 'mfc': 'gray'},
		thalidwashout1={'ls': '-', 'marker': '.', 'mfc': 'gray'},
		thalidwashout2={'ls': '-', 'marker': '.', 'mfc': 'gray'},
		vehiclewashout1={'ls': '-', 'marker': '.', 'mfc': 'gray'},
		vehiclewashout2={'ls': '-', 'marker': '.', 'mfc': 'gray'})

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
	figpath = os.path.join(studydir, 'Analysis', 'compare_conditions_pairwise_%s.png' % control)
	fig.savefig(figpath)

def compare_conditions_by_postdate1():

	df = load_experiment()
	df = df[df.freq<40000]
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

	colors = dict(tnfa='r', vehicle='b')
	fig = plt.figure(figsize=(10, 8.5))

	for j, date in enumerate(udate):

		ax = fig.add_subplot(2, 2, j+1)

		y_ = y.ix[date]
		yerr_ = yerr.ix[date]

		x = range(len(y_.index))

		ax.errorbar(x, y_pre, yerr=yerr_pre, color='k', ecolor='k', label='preinjection')

		for ((i, yc), (i, yerrc)) in zip(y_.iteritems(), yerr_.iteritems()):
			ax.errorbar(x, yc, yerr=yerrc, color=colors[i], label=i)
			ax.set_title('%u days post injection' % date)
			ax.set_xticklabels((y_.index / 1000.).astype(int))
			# ax.set_xlim([-1, len(x)])
			format_axis(ax)

	ax.legend()
	figpath = os.path.join(studydir, 'Analysis', 'compare_postdate1.png')
	fig.savefig(figpath)

def compare_conditions_by_day(conditions=('tnfa', 'vehicle')):

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

	figpath = os.path.join(studydir, 'Analysis', 'compare_conditions_by_day_%s.png' % '_'.join(conditions))
	fig.savefig(figpath)

	return fig

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
	ax.axhline([1.0], color='r', ls='--')
	ax.legend(loc='best')
	ax.set_title('Start age: %u' % df.age.min())

	figdir = os.path.join(studydir, 'Analysis')
	if not os.path.exists(figdir):
		os.mkdir(figdir)
	figpath = os.path.join(figdir, 'compare_conditions.png')
	fig.savefig(figpath)

	return fig

def single_subject_conditionmeans():

	fig = plt.figure(figsize = (10, 8));
	ax = fig.add_subplot(111);
	figdir = os.path.join(studydir, 'Analysis')
	if not os.path.exists(figdir):
		os.mkdir(figdir)
	for animalID in animalIDs:

		fpaths = glob.glob(os.path.join(studydir, 'gapdetection', '%s*.csv' % animalID))
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

	return fig

def single_subject_dailyresults(df=None, condition='all', cond_color=True):

	if condition=='all':
		cond_patt = ''
	else:
		cond_patt = '_'+condition+'_'
	
	if df is None:
		df = load_experiment(onlygood=False)

	animalIDs = np.unique(df.animalID).values
	for animalID in animalIDs:
		fig = plt.figure(figsize = (6, 6));
		ax = fig.add_subplot(111);
		fpaths = glob.glob(os.path.join(studydir, 'gapdetection', '%s*%s*.csv' % (animalID, cond_patt)))
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
		ax.legend(loc = 'upper left')
		format_axis(ax)

		figpath = os.path.join(studydir, 'Analysis', 'dailyresults_%s_%s.png' % (condition, animalID))
		fig.savefig(figpath)
		# ax.cla();

def single_subject_startleampl():

	fig = plt.figure()
	animalIDs = pd.read_csv(os.path.join(studydir, 'dobs.csv'))['animalID'].values
	for animalID in animalIDs:
		fnames = glob.glob(os.path.join(studydir, 'data', 'Gap', Gap.reverse_animalID(animalID), '*.txt'))
		nfnames = len(fnames)
		nrows, ncols = misc.get_subplot_grid(nfnames)
		for i, fname in enumerate(fnames):
			ax = fig.add_subplot(nrows, ncols, i+1)
			df = Gap.txt2pd(fname)
			Gap.plot_startleampl(df, ax)
			ax.set_ylim([0, 8])
			ax.set_title(os.path.split(fname)[1])
		figpath = os.path.join(studydir, 'Analysis', 'startleampl_%s.png' % animalID)
		fig.savefig(figpath)
		fig.clf();

	return fig

def format_axis(ax):
	ax.set_xlabel('Frequency (kHz)')
	ax.set_ylabel('Gap ratio')
	ax.set_ylim([0, 1.5])
	ax.axhline([1.0], color = 'r', ls = '--')
	# ax.legend()



def analyze():
	try:
		fig = compare_conditions()
		plt.close(fig)
		fig = compare_conditions_by_day()
		plt.close(fig)
	except:
		pass
	fig = single_subject_conditionmeans()
	plt.close(fig)
	fig = single_subject_startleampl(); plt.close(fig)
	single_subject_dailyresults(cond_color = False)
	df = load_experiment()
	# useful for determining which of the initial sessions should be included
	for key, df_ in df.groupby(('animalID', 'condition')):
		single_subject_dailyresults(df_, condition = key[1], cond_color = False)
