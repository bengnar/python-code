import os, glob
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import datetime
import misc


studydir = '/Volumes/BOB_SAGET/Gap_detection'
dobs = np.loadtxt(os.path.join(studydir, 'dobs.csv'), 'S', delimiter = ',', skiprows = 1)

df = df = pd.read_csv(os.path.join(studydir, 'Processed', 'gap_detect_data.csv'))

def load_sessions():
	cagepaths = glob.glob(os.path.join(studydir, 'Cages', '*'))

	postdurs = dict(_20130301 = 50, _20130304 = 50, _20130305 = 20, _20130306 = 50, _20130307 = 20, _20130308 = 50, _20130311 = 50, _20130313 = 20, _20130314 = 50, _20130315 = 20, _20130318 = 50, _20130319 = 50, _20130320 = 20)
	
	df = []
	for cagepath in cagepaths:
		_, cage = os.path.split(cagepath)
		gen = cage[:2]
		sesspaths = glob.glob(os.path.join(cagepath, '*'))
		for sesspath in sesspaths:
			_, sess = os.path.split(sesspath)
			thisdate = datetime.date(int(sess[:4]), int(sess[4:6]), int(sess[6:]))
			animalpaths = glob.glob(os.path.join(sesspath, '%s*.txt' % gen))
			for animalpath in animalpaths:
				absol, relat = os.path.split(animalpath)
				_, color, _, _, _, _ = relat.split('_')
				ID = gen+color
				
				df_ = load_session(animalpath, postdur = postdurs['_'+sess])
				
				dob = dobs[dobs[:, 0]==ID, 1:].astype(int).flatten()
				df_['age'] = (thisdate - datetime.date(dob[0], dob[1], dob[2])).days
				df_['gen'] = gen
				df_['ID'] = ID
				df_['sess'] = sess
			
				df.append(df_)

	df = pd.concat(df)
	df.to_csv(os.path.join(studydir, 'Processed', 'gap_detect_data.csv'))
	return df

def convert_all():

	cagepaths = glob.glob(os.path.join(studydir, 'Cages', '*'))

	for cagepath in cagepaths:
		_, cage = os.path.split(cagepath)
		gen = cage[:2]
		sesspaths = glob.glob(os.path.join(cagepath, '*'))
		for sesspath in sesspaths:
			animalpaths = glob.glob(os.path.join(sesspath, '%s*.txt' % gen))
			for animalpath in animalpaths:
				absol, relat = os.path.split(animalpath)
				newpath = os.path.join(studydir, 'Processed', 'data', os.path.splitext(relat)[0]+'.h5')
				if not os.path.exists(newpath):
					print animalpath
					x = np.loadtxt(animalpath, skiprows = 1)
					f = h5py.File(newpath)
					convert_to_hdf5(x, f)
					f.close()
	
def convert_to_hdf5(x, f):

	freq = np.int32(x[:, 0]/1000) # convert to kHz
	cue = x[:, 1].astype('bool')
	gaplen = np.int32(x[:, 2]*1000) # convert to milliseconds
	holdtime = np.int32(x[:, 3]*1000) # convert to milliseconds
	startle = x[:, 4:]
	startle_amp = calc_max_startle(startle)

	f.create_dataset('freq', data = freq, compression = 'gzip')
	f.create_dataset('cue' , data = cue, compression = 'gzip')
	f.create_dataset('gaplen' , data = gaplen, compression = 'gzip')
	f.create_dataset('holdtime' , data = holdtime, compression = 'gzip')
	f.create_dataset('startle' , data = startle, compression = 'gzip')
	

def load_session(animalpath, postdur = 0.05):
	
	x = np.loadtxt(animalpath, skiprows = 1)

	freq = np.int32(x[:, 0]/1000) # convert to kHz
	cue = x[:, 1].astype('bool')
	gaplen = np.int32(x[:, 2]*1000) # convert to milliseconds
	holdtime = np.int32(x[:, 3]*1000) # convert to milliseconds
	startle = x[:, 4:]
	startle_amp = calc_max_startle(startle)
	trial = np.arange(1, x.shape[0]+1)
	
	# make data frame
	df = pd.DataFrame(data = dict(trial = trial, freq = freq, cue = cue, gaplen = gaplen, holdtime = holdtime, samplitude = startle_amp))
	
	# add post duration
	df['postdur'] = postdur

	# get startle mean for uncued in order to make startle ratio
	guncued = df[~df.cue].groupby(('freq'))
	uncued = guncued.agg(np.mean).samplitude
	
	# add column for startle ratio
	df['sratio'] = pd.Series(df['samplitude'].values / uncued[df['freq']].values, index = df.index)

	return df
			
def calc_max_startle(startle_trace):

	startle_trace = startle_trace[:, 150:300]
	maxstartle = startle_trace.max(1)
	minstartle = startle_trace.min(1)
	startle_amp = maxstartle-minstartle
	
	return startle_amp

def make_contactsheets():

	df = pd.read_csv(os.path.join(studydir, 'Processed', 'gap_detect_data.csv'))

	gp = df.groupby(('ID', 'sess'))
	for (ID, sess), df_ in gp:
		figpath = os.path.join(studydir, 'Sheets', '%s_%i_sheet.png' % (ID, sess))
		if not os.path.exists(figpath):
			f = load_startle(ID, str(sess))
			fig = make_contactsheet(f, df_)
			f.close()
			fig.savefig(figpath)
			plt.close(fig)
			
def make_contactsheet(f, df_):
	'''
	plots for one animal
	pandas dataframe containing data for this animal
	'''
	startle = f['startle'].value
	ugaplens = np.unique(df_.gaplen.values)
	ngaplens = ugaplens.size
	
	fig = plt.figure()
	fig.suptitle(df_.head(1).sess)
	# plot raw traces for each gap length
	ax1 = []
	for i in xrange(ngaplens):
		ax1.append(fig.add_subplot(3, ngaplens, i+1))
		ix = df_.gaplen==ugaplens[i]
		misc.plot_matrix(np.abs(startle[ix, 100:400].T), ax = ax1[-1])
		ax1[-1].set_title(ugaplens[i])
	
	[a.set_yticklabels('') for a in ax1[1:]]
	misc.sameyaxis(ax1)
	ax2 = fig.add_axes((0.1, 0.1, 0.8, 0.5))
	df_.groupby('gaplen').agg(np.mean).sratio.plot(kind = 'bar', ax = ax2)

	return fig
	
def load_startle(ID, sess):
	
	yr, mo, da = int(sess[:4]), int(sess[4:6]), int(sess[6:])
	gen, color, cageno = ID[:2], ID[2:3], ID[3:]
	fpath = os.path.join(studydir, 'Processed', 'data', '%s_%s%s_%i_%i_%i_0.h5' % (gen, color, cageno, mo, da, yr))
	try:
		f = h5py.File(fpath, 'r')
	except IOError:
		print 'Could not find %s!' % fpath
	
	return f