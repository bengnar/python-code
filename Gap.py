import os, glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import datetime
import misc

studydir = '/Volumes/BOB_SAGET/TNFalpha/salicylate'

def combine_all():

    fpaths = glob.glob(os.path.join(studydir, 'fileconversion', '*.h5'))
    df = []
    for fpath in fpaths:
        absol, relat = os.path.split(fpath)
        animalID, gen, exp, mo, da, yr, _ = relat.split('_')
        d = pd.HDFStore(fpath, 'r')
        df_ = d['df']
        d.close()
        
        df_['animalID'] = [animalID]*len(df_)
        df_['gen'] = [gen]*len(df_)
        df_['exp'] = [exp]*len(df_)
        df_['sess'] = [relat]*len(df_)
        
        df_ = add_amplitude(df_)
        df_ = add_gapratio(df_)
        df.append(df_)

    df = pd.concat(df)
    
    return df    

def plot_indiv_results(df):

    df = df[df.cued==1]
    gp = df.groupby(['animalID', 'exp'])
    for k, v in gp:
       print k, v
       gp2 = v.groupby(['sess', 'freq'])
       gp2.gapratio.apply(np.mean)
    
def plot_group_results():
    df = df[df.cued==1]
    plt.close('all')
    gp = df.groupby(['animalID', 'freq', 'exp'])
    results_mean = gp.gapratio.apply(np.mean)
    results_err = gp.gapratio.apply(np.std)
    
    misc.pd_errorbar(results_mean, results_err)
    fig = plt.figure();
    ax = fig.add_axes([0.1, 0.3, 0.8, 0.6])
    
    results.plot(kind = 'bar', ax = ax)
    return

def add_gapratio(df):
     
    if 'ampl' not in df.keys():
        df = add_amplitude(df)

    uncued_amp = df.ampl[df.cued==0].mean()
    df['gapratio'] = df.ampl / uncued_amp
    
    return df
    
def add_amplitude(df):
    df['ampl'] = find_amplitude(df)
    return df
    
def find_amplitude(df):
    resps = []
    for i in df.index:
        if df.cued[i]==0:
            resp = df.resp[i][100:300]
        elif df.cued[i]==1:
            resp = df.resp[i][200:400]
        resps.append(resp.max()-resp.min())
    return resps

def convert_to_pd_all():
	fpaths = glob.glob(os.path.join(studydir, 'data', '*.txt'))
	for fpath in fpaths:
	    if not os.path.exists(fpath.replace('data', 'fileconversion').replace('.txt', '.h5')):
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

def calc_max_startle(gp):

	for i in xrange(len(gp)):
		x = gp[gp.index[i]][150:300]
		amp = max(x) - min(x)

	return np.mean(amp)
	# startle_trace = startle_trace[150:300]
	# maxstartle = max(startle_trace)
	# minstartle = min(startle_trace)
	# startle_amp = maxstartle-minstartle
	
	return startle_amp
