import os, glob, re
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt
import misc

# basedir = '/Volumes/BOB_SAGET/Fmr1_Heesoo/ABR'
basedir = '/Users/robert/Documents/Work/Bao/Fmr1_Heesoo/ABR'
dtype = np.dtype([('gen', 'S2'), ('exp', 'S3'), ('animal', 'i4'), ('freq', 'f4'), ('atten', 'f4'), ('lat', '5f4'), ('ampl', '5f4')])

plt_opts = {'wt' : {'color' : 'b', 'x_offset' : 0, 'y_offset' : 0.00001}, 'ko' : {'color' : 'r', 'x_offset' : 0.09, 'y_offset' : 0}}

def manual_threshold_all():
	
	fnames = glob.glob(os.path.join(basedir, 'fileconversion', '*.npz'))
	fnames = np.array(fnames)
	fnames = fnames[np.random.permutation(fnames.size)]
		
	fig = plt.figure(figsize = [8.3875,  8.825])
	fig.subplots_adjust(bottom = 0, top = 1, right = 1)
	ax = fig.add_subplot(111)
	
	for fname in fnames:
		d = np.load(fname)['arr_0']
		manual_threshold(d, ax = ax)
		
def manual_peaks_all():

	fig = plt.figure(figsize = [8.3875,  8.825])
	fig.subplots_adjust(bottom = 0.04, top = 0.96, right = 1)
	ax = fig.add_subplot(111)

	fnames = glob.glob(os.path.join(basedir, 'fileconversion', '*.npz'))
	for fname in fnames:
		d = np.load(fname)['arr_0']
		manual_peaks(d, ax)
		
		
def plot_abr_attens(t, abr, attens, y_offset_mag = 0.0003, ax = None, **kwargs):
	
	ax, fig = misc.axis_check(ax)
	
	nattens = abr.shape[1]
	nsamples = t.size
	# abr = abr - np.tile(abr[0, :], (nsamples, 1))
	y_offset = np.tile(np.linspace(y_offset_mag, 0, nattens), (nsamples, 1))
	
	abr = abr + y_offset

	ax.plot(t, abr, **kwargs)
	ax.legend(attens)
	
	plt.show()

def manual_threshold(d, ax = None):

	ax, fig = misc.axis_check(ax)
	savepath = os.path.join(basedir, 'thresholds.npz')
	dtype = np.dtype([('gen', 'S2'), ('exp', 'S3'), ('animal', 'i4'), ('freq', 'f4'), ('thresh', 'f4')])

	if os.path.exists(savepath):
		threshold = np.load(savepath)['arr_0']
	else:
		threshold = np.empty(0, dtype = dtype)

	ufreqs = np.unique(d['freq'])
	nsamples = 244
	duration = 0.099424
	sample_rate = nsamples / duration
	t = np.arange(0, duration, 1./sample_rate)

	for freq in ufreqs:

		d_ = d[d['freq']==freq]
		
		nattens = d_.size
		y = filter_abr(d_['data'].T)

		plot_abr_attens(t, y, d_['atten'], ax = ax)
		
		thresh = int(raw_input('Threshold--> '))
		
		ax.cla()
		
		threshold.resize(threshold.size+1)
		threshold[-1] = np.array((d_['gen'][0], d_['exp'][0], d_['animal'][0], freq, thresh), dtype = dtype)

	np.savez(savepath, threshold)

def manual_peaks(D, ax):

	savepath = os.path.join(basedir, 'peaks.npz')
	dtype = np.dtype([('gen', 'S2'), ('exp', 'S3'), ('animal', 'i4'), ('freq', 'f4'), ('attens', '2f4'), ('peaks', '(2,5)f4')])

	if os.path.exists(savepath):
		peaks = np.load(savepath)['arr_0']
	else:
		peaks = np.empty(0, dtype = dtype)

	# ufreqs = np.unique(D['freq'])
	ufreqs = [8000., 16000.]
	
	nsamples = 244
	duration = 0.099424
	sample_rate = nsamples / duration
	t = np.arange(0, duration, 1./sample_rate)
	peak_times = (np.array([0.9, 1.6, 2.3, 3.2, 4.2]) / 1000.)+0.001
	for freq in ufreqs:
		
		d_ = D[D['freq']==freq]
		
		attens_ = d_['atten']
		atten = np.sort(d_['atten'])[np.array([0, 5])]
		d_ = d_[(np.tile(d_['atten'], (2, 1)) == np.tile(atten, (d_['atten'].shape[0], 1)).T).any(0)]
		
		nattens = d_.size

		y = d_['data'].T
		y = y - np.tile(y[0, :], (244, 1))
		y_offset = np.tile(np.linspace(0.00005, 0, nattens), (244, 1))
		
		y = y + y_offset
		
		# ax = fig.add_subplot(111)
		ax.plot(t, y)
		for i in range(len(peak_times)):
			ax.axvline(peak_times[i], color = 'r', ls = '--')
		ax.legend(d_['atten'])
		plt.show()
		
		xy = plt.ginput(n = -1, timeout = -1)
		
		ax.cla()
		
		peak_ = np.array(xy)[:, 0]
		peak_ = np.array(peak_)
		peak_.resize(peak_.size/5, 5)
		peak = np.empty((2, 5)) * np.nan
		peak[:peak_.shape[0], :] = peak_
		
		attens = np.empty(2) * np.nan
		attens[:nattens] = d_['atten']
		print peak
		peaks.resize(peaks.size+1)
		peaks[-1] = np.array((d_['gen'][0], d_['exp'][0], d_['animal'][0], freq, d_['atten'], peak), dtype = dtype)

	np.savez(savepath, peaks)


def convert_all():
	
	fnames = glob.glob(os.path.join(basedir, 'data', '*.txt'))

	for fname in fnames:
		
		print fname
		convert(fname)

def convert(fname):

	relat, absol = os.path.split(fname)
	fname_, _ = os.path.splitext(absol)
	gen, exp, animal = fname_.split('_')
	
	x = open(fname)
	x = x.readlines()
	x = x[13:]
	p = re.compile('(\d+)')

	dtype = np.dtype([('gen', 'S2'), ('exp', 'S3'), ('animal', 'i4'), \
		('freq', 'f4'), ('atten', 'f4'), \
		('data', '244f4')])
	
	X = np.array([], dtype = dtype)
	done = False
	while len(x)>0:

		npts = np.int32(p.findall(x[3])[0])
		assert npts == 244
		freq = np.int32(p.findall(x[10])[0])
		atten = np.int32(p.findall(x[11])[0])
		data_ = x[12:12+npts]
		data = [float(d[:-2]) for d in data_]
	
		X.resize(X.size+1)
		X[-1] = np.array((gen, exp, animal, freq, atten, data), dtype = dtype)

		x = x[12+npts+1:]

	savepath = os.path.join(basedir, 'data', '_'.join((gen, exp, str(animal), 'ABR'))+'.npz')
	np.savez(savepath, X)
	

def filter_abr(x, duration = 0.00999424, nsamples = 244, Wn = 0.02, btype = 'high'):

	b, a = butter(10, Wn = Wn, btype = btype)
	if len(x.shape)==1:
		x_filt = filtfilt(b, a, x)
	elif len(x.shape)>1:
		x_filt = np.empty_like(x)
		for i in range(x.shape[1]):
			x_filt[:, i] = filtfilt(b, a, x[:, i])
	return x_filt

def plot_peak_gen_vs_exp(x, measure = 'ampl'):
	'''
	'''
	gens = ['wt', 'ko']
	exps = ['nai', 'exp']
	freqs = [8000., 16000.]

	plt_opts = {'wt' : {'color' : 'b', 'x_offset' : 0}, 'ko' : {'color' : 'r', 'x_offset' : 1}}

	fig = plt.figure(figsize = (14, 8))
	ax = []
	for i in range(10):
		ax.append(fig.add_subplot(2, 5, i+1))
		ax[-1].set_title('Peak %i' % ((i%5)+1))

	a = np.arange(5)
	for f, freq in enumerate(freqs):
		x2 = x[x['freq']==freq]
		attens = np.unique(x2['atten'])
		for atten in attens[:1]:
			x3 = x2[x2['atten']==atten]
		
			for gen in gens:
				ampl = np.empty((len(exps), 5))
				ampl_err = np.empty((len(exps), 5))
				for j, exp in enumerate(exps):
					x4 = x3[np.vstack((x3['gen']==gen, x3['exp']==exp)).all(0)]
					ampl[j, :] = (x4[measure]*10).mean(0)
					ampl_err[j, :] = (x4[measure]*10).std(0) / np.sqrt(x4.size)
				for i in range(ampl.shape[1]):
					ampl_ = ampl[:, i]
					ampl_err_ = ampl_err[:, i]
					ax[(f*5)+i].errorbar(np.arange(2), ampl_, yerr = ampl_err_, color = plt_opts[gen]['color'])

	misc.sameyaxis(ax)
	misc.samexaxis(ax, [-1, 2])
	for i in range(10):
		ax[i].set_xticks([0, 1])
		ax[i].set_xticklabels(exps)
	for i in range(1, 10):
		ax[i].set_yticklabels([])

	plt.show()	

def plot_threshold_gen_vs_exp(x):
	'''
	'''
	gens = ['wt', 'ko']
	exps = ['nai', 'exp']
	freqs = np.unique(x['freq'])[:-1]

	plt_opts = {'wt' : {'color' : 'b', 'x_offset' : 0}, 'ko' : {'color' : 'r', 'x_offset' : 1}}

	fig = plt.figure(figsize = (14, 4))
	ax = []
	for i, freq in enumerate(freqs):
		ax.append(fig.add_subplot(1, 4, i+1))
		ax[-1].set_title('%u kHz' % (freq/1000))
		for j, gen in enumerate(gens):
			Y = np.empty(2)
			Y_err = np.empty(2)
			for k, exp in enumerate(exps):
				x_ = x[np.vstack((x['gen']==gen, x['exp']==exp, x['freq']==freq)).all(0)]
				y = x_['thresh_calib'].mean()
				y_err = x_['thresh_calib'].std() / np.sqrt(x_.size)
				
				Y[k] = y
				Y_err[k] = y_err
		
			ax[-1].errorbar([1, 2], Y, yerr = Y_err, color = plt_opts[gen]['color'])
	
	misc.sameyaxis(ax)
	misc.samexaxis(ax, [0, 3])
	for i in range(4):
		ax[i].set_xticks([1, 2])
		ax[i].set_xticklabels(exps)
	for i in range(1, 4):
		ax[i].set_yticklabels([])

	plt.show()	

def calc_latency_delay():

	dtype = dtype = np.dtype([('gen', 'S2'), ('exp', 'S3'), ('animal', 'i4'), ('freq', 'f4'), ('atten', 'f4'), ('lat', '5f4'), ('lat_del', '5f4'), ('ampl', '5f4')])
	
	x2 = np.empty(0, dtype = dtype)
	for x_ in x:
		lat_del = np.diff(np.hstack((0, x_['lat'])))
		x2.resize(x2.size+1)
		x2[-1] = np.array((x_['gen'], x_['exp'], x_['animal'], x_['freq'], x_['atten'], x_['lat'], lat_del, x_['ampl']), dtype = dtype)

	np.savez(os.path.join(basedir, 'peaks.npz'), x2)


