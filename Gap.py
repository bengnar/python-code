import os, glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import datetime
import misc; reload(misc);

studydir = '/Volumes/BOB_SAGET/TNFalpha/salicylate/Gap'
freqs = np.array([5, 7.1, 10, 14.1, 20, 28.3, 40])
ls = dict(salicylate = '--', control = '-')
colors = dict(salicylate = 'b', control = 'g')

def combine_all(cage = ''):

    fpaths = glob.glob(os.path.join(studydir, 'fileconversion', '*%s*.h5' % cage))
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
        df_['trial'] = range(len(df_))

        df = df[:100]
        df_ = add_amplitude(df_)
        df_ = add_gapratio(df_)
        df.append(df_)

    df = pd.concat(df)
    df.index = np.arange(len(df))
    
    return df

def plot_indiv_results(df):

    df = df[df.cued==1]
    gp = df.groupby('animalID')
    for animalID, v1 in gp:
        print animalID
        fig = plt.figure();
        ax = fig.add_subplot(111);
        gp1 = v1.groupby('exp')
        for exp, v2 in gp1:
            gp2 = v2.groupby(['freq', 'sess'])
            gapratio = gp2.gapratio.apply(np.mean)
            gapratio.unstack().plot(marker = '.', linestyle = ls[exp], ax = ax)
            #gapratio.to_csv(os.path.join(studydir, 'forshaowen', '%s_%s.csv' % (animalID, exp)))

        ax.legend(loc = 'lower left')
        ax.set_ylim([0, 1.3])
        fig.savefig(os.path.join(studydir, 'Sheets', '_'.join([animalID, 'sheet'])))
        plt.close(fig);
    
        

def get_first(df):
    x = df.ix[df.index[0]]
    return x

def apply_sem(df):
    return np.std(df) / np.sqrt(len(df))


def plot_group_results(df):

    df = df[df.cued==1]

    # group by session to get session mean gap detection
    gp_sess = df.groupby(['sess', 'freq'])
    gapratio_sess = gp_sess.agg(dict(gapratio = np.mean, freq = get_first, gen = get_first, exp = get_first, animalID = get_first))

    # combine sessions for each animal to give animal performance by condition
    gp_animal = gapratio_sess.groupby(['animalID', 'exp', 'freq'])
    gapratio_animal = gp_animal.agg(dict(gapratio = np.mean, freq = get_first, gen = get_first, exp = get_first))
    gp_exp = gapratio_animal.groupby(['freq', 'exp'])
    gapratio_mean = gp_exp.gapratio.apply(np.mean).unstack()
    gapratio_sem = gp_exp.gapratio.apply(apply_sem).unstack()
    fig = plt.figure();
    ax = fig.add_subplot(111);
    for (k, v), (kerr, verr) in zip(gapratio_mean.iteritems(), gapratio_sem.iteritems()):
        misc.errorfill(np.arange(freqs.size), v, verr, ax = ax, color = colors[k], marker = '.', ls = '-', label = k)

    ax.legend()

    maxy = ax.get_ylim()[1]
    ax.set_ylim([0, maxy])
    ax.set_ylabel('Gap ratio')
    ax.set_xlabel('Frequency (kHz)')
    ax.set_xticks(range(freqs.size))
    ax.set_xticklabels(freqs)

    return ax

def calc_gap_ratio():
    df = df[df.cued==1]
    
    sess_group = df.groupby(['sess', 'freq'])
    sess_gapratio = sess_group.agg(dict(gapratio = np.mean, exp = get_first))
    exp_group = sess_gapratio.groupby(['exp', 'freq'])
    gap_mean = exp_group.gapratio.apply(np.mean)
    plt.close('all')
    gap_mean.plot()

    gp = df.groupby(['animalID', 'freq', 'exp'])
    
    
    results_mean = gp.gapratio.apply(np.mean)
    results_err = gp.gapratio.apply(np.std)
    
    misc.pd_errorbar(results_mean, results_err)
    fig = plt.figure();
    ax = fig.add_axes([0.1, 0.3, 0.8, 0.6])
    
    results.plot(kind = 'bar', ax = ax)
    return

def add_gapratio(df):
     
    if 'startleampl' not in df.keys():
        df = add_amplitude(df)

    uncued_amp = df.startleampl[df.cued==0].mean()
    df['gapratio'] = df.startleampl / uncued_amp
    
    return df
    
def add_amplitude(df):
    df['startleampl'] = find_amplitude(df)
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
        convert_to_pd(fpath)

def convert_to_pd(fpath):

    outpath = fpath.replace('data', 'fileconversion').replace('.txt', '.h5')
    if not os.path.exists(outpath):
        print fpath
        header = np.loadtxt(fpath, 'S', delimiter = '\t')[0]
        startletrace_ix = len(header)
        x = np.loadtxt(fpath, 'f', skiprows = 1)
        freq = x[:, 0]
        cued = x[:, 1]
        ampl = x[:, 2]
        holdtime = x[:, 4]
        resp = [i for i in x[:, 5:]]
        df = pd.DataFrame(dict(freq = freq, cued = cued, ampl = ampl, holdtime = holdtime, resp = resp))
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


# # add KO and experiment type to file names
# dates = dict(salicylate = ['5_30_2013', '6_3_2013', '6_4_2013'], control = ['5_23_2013', '5_31_2013', '6_1_2013'])
# fpaths = glob.glob(os.path.join(studydir, 'data', '*.txt'))
# for fpath in fpaths:
#    absol, relat = os.path.split(fpath)
#    animalID, da, mo, yr, ext = relat.split('_')
#    for k, v in dates.iteritems():
#        for v_ in v:
#            if v_ in relat:
#                exp = k
               
#    newrelat = '_'.join([animalID, 'KO', exp, da, mo, yr, ext])
#    shutil.move(fpath, os.path.join(absol, newrelat))
    
    
    
    
    
    
    
