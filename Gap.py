import os, glob
import h5py
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import datetime
import misc


studydir = '/Volumes/BOB_SAGET/TNFalpha/tinnitus/behavior/Gap'
# dobs = np.loadtxt(os.path.join(studydir, 'dobs.csv'), 'S', delimiter = ',', skiprows = 1)

# df = pd.read_csv(os.path.join(studydir, 'Processed', 'gap_detect_data.csv'))

def convert_to_pd_all():
	fpaths = glob.glob(os.path.join(studydir, 'data', '*.txt'))
	for fpath in fpaths:
		convert_to_pd(fpath)

def convert_to_pd(fpath):

	header = np.loadtxt(fpath, 'S', delimiter = '\t')[0]
	startletrace_ix = len(header)
	x = np.loadtxt(fpath, 'f', skiprows = 1)
	print header

	freq = x[:, 0]
	cued = x[:, 1]
	ampl = x[:, 2]
	holdtime = x[:, 4]
	resp = [i for i in x[:, 5:]]
	# resp = x[:, 5:]
	df = pd.DataFrame(dict(freq = freq, cued = cued, ampl = ampl, holdtime = holdtime, resp = resp))

	outpath = fpath.replace('data', 'fileconversion').replace('.txt', '.h5')
	d = pd.HDFStore(outpath)
	d['df'] = df
	d.close()

def make_contact_sheet(df):

	gp1 = df.groupby(('freq', 'cued'))
	max_startle = gp1.resp.apply(calc_max_startle)
	mean_resp = gp1.resp.apply(np.mean)
	max_startle = mean_resp.apply(calc_max_startle)
	mean_startle = max_startle.apply(np.mean)
	err_startle = max_startle.apply(np.std)



	# uncued_mean = 

def calc_max_startle(gp):

	for i in xrange(len(gp)):
		x = gp[gp.index[i]][150:300]
		amp = max(x) - min(x)

	return mean(amp)
	# startle_trace = startle_trace[150:300]
	# maxstartle = max(startle_trace)
	# minstartle = min(startle_trace)
	# startle_amp = maxstartle-minstartle
	
	return startle_amp
