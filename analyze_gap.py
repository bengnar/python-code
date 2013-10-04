import os, glob, shutil
import numpy as np
from scipy import stats as st
import pandas as pd
import Gap; reload(Gap)
import matplotlib as mpl
from matplotlib import pyplot as plt
import misc; reload(misc)
import pdb

color_cycle = mpl.rcParams['axes.color_cycle']
colors = dict(control = 'k', salicylate = 'r', prenihl = 'k', postnihl = 'g', postnihl5d = 'y', thalid = 'r', vehicle = 'b', \
	preinjection = 'k', preinjectiontnfa = 'k', preinjectionvehicle = 'k', pre = 'k', tnfa = 'r', thalidwashout = 'y', vehiclewashout = 'm')
lss = dict(control = '-', salicylate = '--', prenihl = '-', postnihl = '--', postnihl5d = '--', thalid = '--', vehicle = '--', \
	preinjection = '-', preinjectiontnfa = '-', preinjectionvehicle = '-', pre = '-', tnfa = '--', thalidwashout = '--', vehiclewashout = '--')

def agg_nanmean(ser):
	return np.mean(ser[~pd.isnull(ser)])

def agg_nansem(ser, axis=0):
	return st.sem(ser[~pd.isnull(ser)])

class GapAnalysis(object):

	def __init__(self, studyID=None, cageID=None):

		# self.basedir = '/Volumes/BOB_SAGET/TNFalpha/tinnitus/behavior'
		self.basedir = '/Users/robert/Desktop'
		# studydir = os.path.join(self.basedir, studyID, cageID)
		if (studyID is not None) and (cageID is not None):
			self.select_study(studyID, cageID)

		self.figdir = os.path.join(self.studydir, 'Analysis')
		if not os.path.exists(self.figdir):
			os.mkdir(self.figdir)

		self.load_experiment()
		self.postdate_bins = [-100, 0, 2, 4, 6, 8, 10, 12]
		self.figformat = 'png'

	def load_experiment(self, conditions = None, onlygood = True):
		'''

		'''
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
				patt = os.path.join(self.studydir, 'gapdetection', '[A-Za-z]*%s*' % cond_patt)
			else:
				patt = os.path.join(self.studydir, 'gapdetection', '*%s*' % cond_patt)
			fnames.extend(glob.glob(patt))

		dfs = []
		if len(fnames)==0:
			print 'No files found with pattern\n%s' % patt
		else:
			for fname in fnames:
				df_ = pd.read_csv(fname)
				dfs.append(df_)
			df = pd.concat(dfs)

			df = df[np.vstack((5000<df.freq, df.freq<40000)).all(0)]
			self.df = df

	def select_study(self, studyID, cageID):

		studydir = os.path.join(self.basedir, studyID, cageID)
		animalIDs = pd.read_csv(os.path.join(studydir, 'dobs.csv'))['animalID'] 
		animalIDs = [i for i in animalIDs if i[0]!='_']
		
		self.animalIDs = animalIDs
		self.studydir = studydir
		self.studyID = studyID
		self.cageID = cageID


	def compare_condition_diffs_by_freq(self):

		df = self.df
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

	def condition_means_by_animal_by_day(self, bins=None):
		'''
		Plots the raw gap detection ratio as a function of frequency 
		averaged over a couple days (specify bins)
		'''
		df = self.df

		df = df[df.freq<28000]
		if bins is None:
			bins = [-100, 0, 2, 4, 6, 8, 10]

		df['postdate_bin'] = pd.cut(df.postdate1, bins=bins)
		gp = df.groupby(('animalID', 'freq', 'postdate_bin'))
		
		condmean = gp.gapratio.agg(agg_nanmean).unstack('postdate_bin')
		condmean = gp.agg(dict(gapratio=agg_nanmean, condition=misc.get_first)).unstack('postdate_bin')
		condsem = gp.agg(dict(gapratio=agg_nansem, condition=misc.get_first)).unstack('postdate_bin')

		animals = np.unique(df.animalID)
		fig, ax = plt.subplots()
		for animal in animals:
			condmean.ix[animal].plot(kind='line', ax=ax)
			cond = df.ix[df['animalID']==animal].iloc[-1]['condition']
			ax.set_title('%s (%s)' % (animal, cond))
			fig.savefig(os.path.join(self.figdir, '%s_%s_means_by_animal_by_day.%s' % (animal, cond, self.figformat)))
			ax.cla()

		plt.close(fig)

	def condition_means_by_day(self, bins=None):

		pass

	def pairwise_compare_condition_by_day_by_freq(self):

		df = self.df

		df_control = df[df.condition==control]
		animalgp_control = df_control.groupby(('animalID', 'freq'))
		animalmeans_control = animalgp_control.agg(dict(gapratio=np.mean))

		df_manip = df[df.postdate1>0]
		if bins is None:
			bins = np.arange(df_manip.postdate1.min()-1, df_manip.postdate1.max()+1, 2)
		df_manip['postdate_bin'] = pd.cut(df_manip.postdate1, bins=bins)
		postdategp = df_manip.groupby(('animalID', 'freq', 'postdate_bin'))
		postdate_manip = postdategp.agg(dict(gapratio=np.mean)).unstack('postdate_bin')
		
		gapdiff = pd.DataFrame(index=animalmeans_control.index)
		for (key, df_) in postdate_manip.iteritems():
			gapdiff[key[1]] = df_ - animalmeans_control.gapratio

		condition_by_animal = df_manip.groupby('animalID').condition.first()

		cond = []
		for (key, value) in gapdiff.iterrows():
			cond.append(condition_by_animal[key[0]])
		gapdiff['condition'] = cond

		# add frequency column (from index)
		gapdiff['freq'] = gapdiff.index.get_level_values('freq')

		condgp = gapdiff.groupby(('condition', 'freq'))
		condmean = condgp.agg(agg_nanmean)
		condsem = condgp.agg(agg_nansem)


	def pairwise_compare_condition_by_day(self, control = 'pre', bins=None):
		'''
		Plots change in gap ratio pairwise for each animal, pre- and post-manipulation
		
		The gap ratio is calculated as the average for each this animal, each frequency
		and then averaged for all frequencies.
		Change in gap ratio is the difference between this value for the pre- and post-
		manipulation days.
		
		bins : Post-manipulation data is averaged in bins of 2 days, by default, but arbitrary bin
			sizes can be used
		'''

		df = self.df
		df_control = df[df.condition==control]
		animalgp_control = df_control.groupby(('animalID', 'freq'))
		animalmeans_control = animalgp_control.agg(dict(gapratio=np.mean))

		df_manip = df[df.postdate1>0]
		if bins is None:
			bins = np.arange(df_manip.postdate1.min()-1, df_manip.postdate1.max()+1, 2)
		df_manip['postdate_bin'] = pd.cut(df_manip.postdate1, bins=bins)
		postdategp = df_manip.groupby(('animalID', 'freq', 'postdate_bin'))
		postdate_manip = postdategp.agg(dict(gapratio=np.mean)).unstack('postdate_bin')
		
		gapdiff = pd.DataFrame(index=animalmeans_control.index)
		for (key, df_) in postdate_manip.iteritems():
			gapdiff[key[1]] = df_ - animalmeans_control.gapratio

		condition_by_animal = df_manip.groupby('animalID').condition.first()

		cond = []
		for (key, value) in gapdiff.iterrows():
			cond.append(condition_by_animal[key[0]])
		gapdiff['condition'] = cond

		condgp = gapdiff.groupby('condition')
		condmean = condgp.agg(agg_nanmean)
		condsem = condgp.agg(agg_nansem)

		styles = dict(tnfa=dict(color='r', hatch='///'), vehicle=dict(color='b', hatch=None))
		x = 2*np.arange(len(condmean.columns))
		width=0.8
		fig, ax = plt.subplots()
		for i,  ((key, y), (_, yerr)) in enumerate(zip(condmean.iterrows(), condsem.iterrows())):
			ax.bar(x-width+(i*width), y, yerr=yerr, width=width, fill=False, hatch=styles[key]['hatch'], ecolor='k', label=key)

		ax.legend(loc = 'best')
		ax.set_xticks(x)
		ax.set_xticklabels(condmean.columns)
		ax.set_xlabel('Days post injection')
		ax.set_ylabel('Change in gap detection ratio')
		ax.axhline(0, color = 'k')
		ax.set_title(os.path.split(self.studydir)[1])

		fig.savefig(os.path.join(self.figdir, 'pairwise_compare_condition_by_day.%s' % self.figformat))

	def pairwise_compare_conditions_by_freq(self, control = 'prenihl'):

		df = self.df

		ucond = np.unique(df.condition).values
		ncond = len(ucond)

		# make animal/condition/freq groups
		animalgp = df.groupby(('animalID', 'condition', 'freq'))
		# take means of gap performance for those groups
		animalmeans = animalgp.agg(dict(gapratio = np.mean))

		# put the control and manipulation columns side-by-side
		condstack = animalmeans.unstack('condition')
		
		# add a column indicating which condition this animal is
		animalcondition = []
		for key, value in condstack.iterrows():
			if not pd.isnull(value['gapratio']['tnfa']):
				animalcondition.append('tnfa')
			else:
				animalcondition.append('vehicle')
		condstack['animalcondition'] = animalcondition

		# make a column that is the result after manipulation
		# (combines multiple columns for each manipulation)
		# subtract the control from the manipulation values (for each frequency)
		manip = condstack.gapratio.vehicle.copy()
		manip.update(condstack.gapratio.tnfa)
		condstack['manip'] = manip
		condstack['manipdiff'] = condstack.manip - condstack.gapratio[control]

		# add a frequency column (turns frequency index into column)
		condstack['freq'] = condstack.index.get_level_values('freq')

		# group changes in gap detection by condition and frequency
		condgp = condstack.groupby(('animalcondition', 'freq'))
		condmeans = condgp.manipdiff.apply(np.mean).unstack('animalcondition')
		condsems = condgp.manipdiff.apply(st.sem).unstack('animalcondition')

		# line plot (change in gap detection VS frequency)
		ufreqs = np.unique(condstack.freq).values
		x = range(len(ufreqs))
		fig, ax = plt.subplots()
		styles = dict(tnfa=dict(color='r', hatch='///'), vehicle=dict(color='b', hatch=None))
		for (key, y), (_, yerr) in zip(condmeans.iteritems(), condsems.iteritems()):
			misc.errorfill(x, y, yerr, ax=ax, label=key, color=styles[key]['color'])
		ax.set_xticks(x)
		ax.set_xticklabels((ufreqs/1000.).astype(int))
		ax.legend()

	def pairwise_compare_conditions(self, control = 'preinjection'):
		'''
		Performs a pairwise comparison between a control condition and several "manipulated"
		conditions where the manipulations were performed on DIFFERENT ANIMALS.
		
		This script plots the paired pre-post difference for each animal, and the manipulation
		group is indicated by different marker type/line styles
		'''

		df = self.df

		ucond = np.unique(df.condition).values
		ncond = len(ucond)

		# make animal/condition/freq groups
		animalgp = df.groupby(('animalID', 'condition', 'freq'))
		# take means of gap performance for those groups
		animalmeans = animalgp.agg(dict(gapratio = np.mean))

		# put the control and manipulation columns side-by-side
		condstack = animalmeans.unstack('condition')
		
		# add a column indicating which condition this animal is
		animalcondition = []
		for key, value in condstack.iterrows():
			if not pd.isnull(value['gapratio']['tnfa']):
				animalcondition.append('tnfa')
			else:
				animalcondition.append('vehicle')
		condstack['animalcondition'] = animalcondition

		# make a column that is the result after manipulation
		# (combines multiple columns for each manipulation)
		# subtract the control from the manipulation values (for each frequency)
		manip = condstack.gapratio.vehicle.copy()
		manip.update(condstack.gapratio.tnfa)
		condstack['manip'] = manip
		condstack['manipdiff'] = condstack.manip - condstack.gapratio[control]

		# add a frequency column (turns frequency index into column)
		condstack['freq'] = condstack.index.get_level_values('freq')

		# calculate change collapsing over freq
		condgp = condstack.groupby('animalcondition')
		condmeans = condgp.manipdiff.apply(np.mean)
		condsems = condgp.manipdiff.apply(st.sem)

		# bar plot (change in gap detection, collapsing frequency)
		fig, ax = plt.subplots()
		for i, ((key, y), (_, yerr)) in enumerate(zip(condmeans.iteritems(), condsems.iteritems())):
			ax.bar(i+0.6, y, yerr=yerr, width=0.8, hatch=styles[key]['hatch'], fill=False, ecolor='k', label=key)
		ax.axhline(0, color='k')
		ax.legend(loc='upper center')
		ax.set_xlim([0, 3])

	def compare_conditions_by_postdate1(self):

		df = self.df

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
				self.format_axis(ax)

		ax.legend()
		figpath = os.path.join(self.figdir, 'compare_postdate1.%s' % self.figformat)
		fig.savefig(figpath)

		plt.close(fig)

	def group_compare_conditions_by_day(self, control='pre'):

		df = self.df

		uconds = np.unique(df.condition).values
		manip = [i for i in uconds if i!=control]

		if not hasattr(df, 'condition_postdate'):
			df = self.add_condition_postdate(df)

		udays = np.unique(df.postdate_bin).values
		ndays = len(udays)-1 # control should have its own date bin

		animalgp = df.groupby(('animalID', 'condition', 'postdate_bin', 'freq'))
		animaldf = animalgp.agg(dict(gapratio = np.mean))
		condgp = animaldf.groupby(level = ('condition', 'postdate_bin', 'freq'))

		ufreqs = np.unique(df.freq).values
		x = np.arange(ufreqs.size)
		y = condgp.agg(np.mean)
		yerr = condgp.agg(st.sem)
		ymean_control = y.ix[control]
		yerr_control = yerr.ix[control]
		ymean_manip = y.ix[manip].unstack('postdate_bin')
		yerr_manip = yerr.ix[manip].unstack('postdate_bin')

		fig, axs = plt.subplots(1, ndays, figsize=(12, 8))
		for i, ((ix, y_), (ix2, yerr_)) in enumerate(zip(ymean_manip.iteritems(), yerr_manip.iteritems())):
			assert ix==ix2
			for (cond, cond_y), (cond, cond_yerr) in zip(y_.unstack('condition').iteritems(), yerr_.unstack('condition').iteritems()):
				axs[i].errorbar(x, cond_y, yerr=cond_yerr, c=colors[cond], label=cond)
			axs[i].errorbar(x, ymean_control.values.flatten(), yerr=yerr_control.values.flatten(), c=colors[control], marker='.', label=control)
			axs[i].set_title('Postdate %s' % ix[1])

		[a.set_xticks(x) for a in axs]
		[a.set_xlim([-0.5, len(x)-0.5]) for a in axs]
		[a.set_ylim([0, 1.5]) for a in axs]
		[a.set_xticklabels(np.int32(ufreqs/1000)) for a in axs]
		axs[0].set_ylabel('Gap ratio')
		[a.set_xlabel('Frequency (kHz') for a in axs]
		axs[-1].legend(loc='best')

		figpath = os.path.join(self.figdir, 'group_compare_conditions_by_day.%s' % self.figformat)
		fig.savefig(figpath)

		plt.close(fig)

	def group_compare_conditions(self, conditions = None, ax = None):

		df = self.df

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
			misc.errorfill(x, y_, yerr_, label = i, ax = ax, marker = '.')
		
		ax.set_xticks(x)
		ax.set_xlim([0, len(x)-1])
		ax.set_xticklabels(np.int32(ufreqs/1000))
		self.format_axis(ax)
		ax.legend(loc='best')
		ax.set_title('Start age: %u' % df.age.min())

		figpath = os.path.join(self.figdir, 'compare_conditions.%s' % self.figformat)
		fig.savefig(figpath)

		plt.close(fig)


	def single_subject_compare_conditions_by_day(self):
		'''
		For each subject, plot the mean gap ratio for each condition
		binned by date. (same as single_subject_compare_conditions, but
		uses binned date as an additional condition)
		'''
		fig, ax = plt.subplots(figsize = (10, 8))

		df = self.df
		df = self.add_condition_postdate(df)

		for animalID in self.animalIDs:

			animaldf = df[df.animalID==animalID]
			condgp = animaldf.groupby(('condition_postdate', 'freq'))

			ufreqs = np.unique(animaldf.freq).values
			x = np.arange(ufreqs.size)
			y = condgp.agg(dict(gapratio = np.mean))
			yerr = condgp.agg(dict(gapratio = st.sem))

			for ((i, y_), (i, yerr_)) in zip(y.unstack().iterrows(), yerr.unstack().iterrows()):
				misc.errorfill(x, y_, yerr_, label = i, ax = ax, marker = 'o')

			ax.set_xticks(x)
			ax.set_xticklabels(np.int32(ufreqs/1000))
			self.format_axis(ax)
			ax.legend()

			figpath = os.path.join(self.figdir, 'single_subject_compare_conditions_by_day_%s.%s' % (animalID, self.figformat))
			fig.savefig(figpath)
			ax.cla();

		plt.close(fig)		

	def single_subject_compare_conditions(self):
		'''
		For each subject, plot the condition means for each condition
		on the same plot for comparison.
		'''
		fig, ax = plt.subplots(figsize = (10, 8))

		for animalID in self.animalIDs:

			animaldf = self.df[self.df.animalID==animalID]
			condgp = animaldf.groupby(('condition', 'freq'))

			ufreqs = np.unique(animaldf.freq).values
			x = np.arange(ufreqs.size)
			y = condgp.agg(dict(gapratio = np.mean))
			yerr = condgp.agg(dict(gapratio = st.sem))

			for ((i, y_), (i, yerr_)) in zip(y.unstack().iterrows(), yerr.unstack().iterrows()):
				# ax.errorbar(x, y_, yerr_, label = i, marker = 'o', color = colors[i], ls = lss[i])
				misc.errorfill(x, y_, yerr_, label = i, ax = ax, marker = 'o', color = colors[i], ls = lss[i])

			ax.set_xticks(x)
			ax.set_xticklabels(np.int32(ufreqs/1000))
			self.format_axis(ax)
			ax.legend()

			figpath = os.path.join(self.figdir, 'single_subject_compareconditions_%s.%s' % (animalID, self.figformat))
			fig.savefig(figpath)
			ax.cla();

		plt.close(fig)

	def single_subject_daily_results(self):
		'''
		For each condition, plots the daily results and a condition mean
		'''

		df = self.df
		df = self.add_condition_postdate(df)
		
		animals = np.unique(df.animalID).values
		ufreqs = np.unique(df.freq).values
		x = range(len(ufreqs))

		fig, ax = plt.subplots()
		for animal in animals:
		
			animaldf = df[df.animalID==animal]

			sessdf = animaldf.filter(['freq', 'gapratio', 'postdate1'])
			sessgapratio = sessdf.pivot(index='freq', columns='postdate1', values='gapratio')

			condgp = animaldf.groupby('condition_postdate')

			for ix, gp in condgp:
				sessgapratio = gp.pivot(index='freq', columns='sess', values='gapratio')
				for i, y in sessgapratio.iteritems():
					ax.plot(x, y, marker='.', label=i)
				ax.plot(x, st.nanmean(sessgapratio, axis=1), c='k', lw=2, marker='.')
				# misc.errorfill(x, st.nanmean(sessgapratio, 1), yerr=st.nanstd(sessgapratio, 1), c='k', lw=2, ax=ax)
				
				ax.set_title('%s: %s' % (animal, ix))
				ax.set_xlim([0, len(x)-1]); ax.set_ylim([0., 1.4])
				ax.axhline(1., c='r', ls='--')
				ax.set_xticks(x)
				ax.set_xticklabels((ufreqs/1000.).astype(int))
				ax.set_xlabel('Frequency (kHz)'); ax.set_ylabel('Gap ratio')
				ax.legend()

				# save out figure
				figname = 'single_subject_dailyresults_%s_%s.%s' % (animal, ix, self.figformat)
				fig.savefig(os.path.join(self.figdir, figname))

				ax.cla()

		plt.close(fig)

	def single_subject_startleampl(self):

		fig = plt.figure()
		animalIDs = pd.read_csv(os.path.join(self.studydir, 'dobs.csv'))['animalID'].values
		for animalID in animalIDs:
			fnames = glob.glob(os.path.join(self.studydir, 'data', 'Gap', Gap.reverse_animalID(animalID), '*.txt'))
			nfnames = len(fnames)
			nrows, ncols = misc.get_subplot_grid(nfnames)
			for i, fname in enumerate(fnames):
				ax = fig.add_subplot(nrows, ncols, i+1)
				df = Gap.txt2pd(fname)
				Gap.plot_startleampl(df, ax)
				ax.set_ylim([0, 8])
				ax.set_title(os.path.split(fname)[1])
			figpath = os.path.join(self.figdir, 'startleampl_%s.%s' % (animalID, self.figformat))
			fig.savefig(figpath)
			fig.clf();

		return fig

	def format_axis(self, ax):
		
		ax.set_xlabel('Frequency (kHz)')
		ax.set_ylabel('Gap ratio')
		ax.set_ylim([0, 1.5])
		ax.axhline([1.0], color='r', ls='--')
		# ax.legend()

	def add_condition_postdate(self, df, bins=None):

		if bins is None: bins = self.postdate_bins

		df['postdate_bin'] = pd.cut(df.postdate1, bins=[-100, 0, 2, 4, 6, 8, 10])
		df['condition_postdate'] = df.condition+df.postdate_bin
		return df

	def analyze(self):
		raise NotImplementedError
