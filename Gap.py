import os, glob
import h5py
import numpy as np
import pandas as pd

studydir = '/Volumes/BOB_SAGET/Gap_detection'


def load_sessions():
	cagepaths = glob.glob(os.path.join(studydir, 'Cages', '*'))

	df = []
	for cagepath in cagepaths:
		_, cage = os.path.split(cagepath)
		gen = cage[:2]
		sesspaths = glob.glob(os.path.join(cagepath, '*'))
		absol, relat = os.path.split(sesspath)
		for sesspath in sesspaths:
			_, sess = os.path.split(sesspath)
			animalpaths = glob.glob(os.path.join(sesspath, '%s*.txt' % gen))
			for animalpath in animalpaths:
				absol, relat = os.path.split(animalpath)
				_, animalID, _, _, _, _ = relat.split('_')
				df_ = load_session(animalpath)
				df_['gen'] = gen
				df_['animal'] = animalID
				df_['sess'] = sess
			
				df.append(df_)

	df = pd.concat(df)
	df['ID'] = df.gen + df.animal
	df.to_csv(os.path.join(studydir, 'Processed', 'gap_detect_data.csv'))
			
def calc_max_startle(startle_trace):

	startle_trace = startle_trace[:, 150:300]
	maxstartle = startle_trace.max(1)
	minstartle = startle_trace.min(1)
	startle_amp = maxstartle-minstartle
	
	return startle_amp
	

def load_session(animalpath):
	
	x = np.loadtxt(animalpath, skiprows = 1)

	freq = x[:, 0]
	cue = x[:, 1].astype('bool')
	gaplen = x[:, 2]
	holdtime = x[:, 3]
	startle_amp = calc_max_startle(x[:, 4:])
	
	df = pd.DataFrame(data = dict(freq = freq, cue = cue, gaplen = gaplen, holdtime = holdtime, samplitude = startle_amp))

	guncued = df[~df.cue].groupby(('freq'))
	uncued = guncued.agg(mean).samplitude
	
	# add column for startle ratio
	df['sratio'] = pd.Series(df['samplitude'].values / uncued[df['freq']].values, index = df.index)

	return df			


# header = ['Frequency', 'Cue', 'Gap length', 'Hold time', 'Startle']
# for cagepath in cagepaths:
# 	_, cage = os.path.split(cagepath)
# 	gen = cage[:2]
# 	sesspaths = glob.glob(os.path.join(cagepath, '*'))
# 	absol, relat = os.path.split(sesspath)
# 	for sesspath in sesspaths:
# 		_, sess = os.path.split(sesspath)
# 		animalpaths = glob.glob(os.path.join(sesspath, '%s*.txt' % gen))
# 		for animalpath in animalpaths:
# 			
# 			absol, relat = os.path.split(animalpath)
# 			_, animalID, _, _, _, _ = relat.split('_')
# 			
# 			# open file
# 			f = open(animalpath, 'r')
# 			header = f.readline()
# 			# read all the data to a numpy array of floats
# 			x = []
# 			for line in f:
# 				x.append(line.split('\t'))
# 			y = np.asarray(x)
# 			f.close()
# 			
# 			# load 2 to end startle
# 			startle = y[:, 4:].astype('float32')
# 			# save freq, cue, gap length
# 			freq = y[:, 0].astype('float32'); cue = y[:, 1].astype('bool'); gaplen = y[:, 2].astype('float32')
# 			# fix holdtime and first startle
# 			holdtime = []
# 			firststartle = []
# 			for y_ in y[:, 3]:
# 				holdtime.append(y_[:2])
# 				firststartle.append(y_[2:])
# 			holdtime = np.asarray(holdtime).astype('float32')
# 			firststartle = np.asarray(firststartle).astype('float32')
# 			
# 			# concatenate all array parts into fixed array
# 			newx = np.hstack((freq[:, np.newaxis], cue[:, np.newaxis], gaplen[:, np.newaxis], holdtime[:, np.newaxis], firststartle[:, np.newaxis], startle))
# 			
# 			# save out file
# 			newf = open(animalpath+'2', 'w')
# 			newf.write(header)
# 			ntrials, npts = newx.shape
# 			for i in range(ntrials): # each trial
# 				newline = []
# 				for j in range(npts):
# 					newline.append(np.str(newx[i, j]))
# 				newline = '\t'.join(newline) + '\r'
# 				newf.write(newline)
# 				
# 			newf.close()
# 			
# import shutil
# for cagepath in cagepaths:
# 	_, cage = os.path.split(cagepath)
# 	gen = cage[:2]
# 	sesspaths = glob.glob(os.path.join(cagepath, '*'))
# 	absol, relat = os.path.split(sesspath)
# 	for sesspath in sesspaths:
# 		_, sess = os.path.split(sesspath)
# 		animalpaths = glob.glob(os.path.join(sesspath, '%s*.txt' % gen))
# 		animalpaths2 = glob.glob(os.path.join(sesspath, '%s*.txt2' % gen))
# 		for animalpath, animalpath2 in zip(animalpaths, animalpaths2):
# 			absol, relat = os.path.split(animalpath)
# 			shutil.move(animalpath, os.path.join(absol, '_'+relat))
# 			shutil.move(animalpath2, animalpath)




