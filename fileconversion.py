# fileconversion.py
# this should be a script that takes the blocks from
# fileconversion and makes them into python structures of a similar nature
# the blocks should be named with 
# prefix -- specifies the stimulus set used (e.g. b for RF or v for vocalizations)
# location ID -- a number indicating which location this was, continuous starting with 1
#		e.g. if you recorded two blocks in the same location but with different stimulus
#		sets, they should have different prefixes, but the same number
# .mat -- the file extension
# note: remove all other characters from the file name, e.g. b15.mat is OK, b15_R3.mat is NOT

#you should have the file0 in a folder and load in the matlab version of the file0 using pylab's
#matlab loader

import pylab
import scipy.io
import numpy as np
import os, shutil, glob, h5py, re
import Spikes; reload(Spikes)
import RF; reload(RF)

# a dictionary mapping the .mat prefix onto the .h5 prefix. for example, b##.mat files are made into RF###.h5 files
prefix = {'b' : 'RF', 'r' : 'RR', 'q' : 'QU', 'v' : 'VOC', 'i' : 'IS', 'P' : 'RF'}

# what kind of electrode was used in this experiment? dictionary keys are electrode descriptions, and items are things like the number of channels in on that electrode, the channel order (if more than one channel per electrode), 
electrodeinfo = dict(\
fourbyfour = dict(nchan=16, npenperblock=4, chanorder=np.array([[1, 7, 13, 14], [3, 4, 10, 16], [2, 8, 12, 11], [6, 5, 9, 15]])), \
tungsten = dict(nchan=4, npenperblock=4, chanorder=None), \
pen = dict(nchan=1, npenperblock=1, chanorder=None)) # a single electrode per block

# any information specific to a stimulus type
stimulusinfo = dict(\
RF = dict(nbins=333), \
RR = dict(nbins=6000), \
VOC = dict(nbins=20000))

# studydir = '/Users/robert/Desktop/tetrode_test'
# studydir = '/Volumes/BOB_SAGET/Fmr1_RR/Sessions'
studydir = '/Volumes/BOB_SAGET/Fmr1_voc'
# studydir = '/Volumes/BOB_SAGET/for_shaowen/'
# studydir = '/Volumes/BOB_SAGET/Fmr1_KO_ising/Sessions/good/'

def fileconvert(sessions, studydir = '/Volumes/BOB_SAGET/Fmr1_voc/'
, v = True, electrodetype = 'tungsten', convertprefix = 'b'):
	'''
	Converts *.mat files into python-friendly *.h5 files.
	Input:
		sessions - this can be either a list of experimental sessions or a the name of a single session as a string
		v - verbose
	Output:
		saves
	
	'''
	if type(sessions) is str:
		sessions = [sessions]
	elif sessions is None:
		sessions = glob.glob(os.path.join(studydir, '*'))
		
	# compile regexp patterns
	p_blocknum = re.compile('(?<=[a-zA-Z])[\d]+')
	
	for session in sessions:
		print session

		#get [b|r|etc]##.mat file names
		files = glob.glob(os.path.join(studydir, 'Sessions', session, 'data', '[%s]*.mat' %  convertprefix))
		nfiles = len(files)
		
		# make the destination directory if it doesn't exist (~/fileconversion/)
		if not os.path.exists(os.path.join(studydir, 'Sessions', session, 'fileconversion')):
			os.mkdir(os.path.join(studydir, 'Sessions', session, 'fileconversion'))

		blockinfo = []

		# parse block file name into
		# blockID - the block name
		# blocknum - the block number
		# blockrep - if multiple blocks were recorded at one site, which one was this
		for ff, file0 in enumerate(files):
			file0 = os.path.basename(file0)
			blockID = os.path.splitext(file0)[0] #removes extension
			tmp = blockID.split('_')
			blocknum = tmp[0]
			if len(tmp)==1:
				blockrep = None
			else:
				blockrep = int(tmp[1])
			blocknum = int(p_blocknum.findall(blocknum)[0])
			blockinfo.append((blocknum, blockrep, blockID))
		blockinfo.sort()
	
		# load cfs and coords
		cfs = load_cfs(session, v = v)
		coords = load_coords(session, v = v)
		
		#loop through penetrations
		for blocknum, blockrep, blockID in blockinfo:

			fileconvert_block(session, blocknum, blockID, blockrep, studydir = studydir, coords = coords, cfs = cfs, v = v, electrodeinfo = electrodeinfo[electrodetype], stimulusinfo = stimulusinfo[prefix[blockID[0]]])

		# end block loop

	# end session loop

def fileconvert_block(session, blocknum, blockID, blockrep, studydir = None, coords = None, cfs = None, v = True, electrodeinfo = electrodeinfo['tungsten'], stimulusinfo = stimulusinfo['RF']):
	
	chanorder = electrodeinfo['chanorder']
	nchan = electrodeinfo['nchan']
	npen = electrodeinfo['npenperblock']
	
	if cfs is None:
		cfs = load_cfs(session)
	if coords is None:
		coords = load_coords(session)
	
	tempfile = scipy.io.loadmat(os.path.join(studydir, 'Sessions', session, 'data', blockID + '.mat')) # load the .mat file0
	Data0 = tempfile['Data'] # save the Data file0 to a temporary variable
	# Data0 = tempfile['data'] # when converting already separated penetrations

	nchan = Data0['trial'][0][0][0][0]['CH'][0].size

	print '\t%s' % blockID

	for cc in range(nchan):
		if chanorder is None:
			siteno = ((blocknum-1)*electrodeinfo['npenperblock']) + cc + 1
		elif len(chanorder.shape) == 2:
			siteno =((blocknum-1)*electrodeinfo['npenperblock']) + (chanorder == cc+1).nonzero()[1][0]+1
			print cc+1, siteno

		blockstr = '%s%3.3i' % (prefix[blockID[0]], siteno)
		if blockrep is not None:
			blockstr = '%s_%1.1i' % (blockstr, blockrep)
		print '\t%s' % blockstr
	
		savepath = os.path.join(studydir, 'Sessions', session, 'fileconversion', '%s%3.3u.h5' % (prefix[blockID[0]], siteno))
		u_ = h5py.File(savepath, 'w')
		unit(u_, Data0, cc, blockID, nbins = stimulusinfo['nbins'])

		add_coords(u_, coords, siteno)
		add_cfs(u_, cfs, siteno)
		
		u_.close()

def unit(u_, Data0, cc, blockID, nbins = 500):
	
	# number of time bins to include in the LFP array
	nlfpsamp = 0
	for tt, trial in enumerate(Data0['trial'][0][0][0]):

		thislfpsamp = trial['LFP'].shape[1]
		if thislfpsamp>nlfpsamp:
			nlfpsamp = thislfpsamp
	
	
	ntrials = Data0['trial'][0][0][0].size # find number of trials
	nstimID = Data0['trial'][0][0][0][0]['Epoch_Value'][0].size
	
	# initialze LFP, spike times, spike trials, spike waveform
	lfp = np.ndarray((0, nlfpsamp), dtype = 'float32')
	spktimes = np.ndarray(0)
	spktrials = np.ndarray(0)
	spkwaveform = np.ndarray((0, 22))
	rast = np.zeros((ntrials, nbins)) 
	# initialize frequency and attenuation IDs
	stimID = np.ndarray((ntrials, nstimID), dtype = 'float32')

	for tt in range(ntrials):
		trial = Data0['trial'][0][0][0][tt]

		thisstimID = np.float32(trial['Epoch_Value'][0])


		# get the LFP for this trial and pad it with nans so it can fit in a matrix (since some of the trials have +/-1 data point for LFP)
		lfpchannel = trial['LFP'][cc]
		lfpchannel = np.concatenate((lfpchannel, np.zeros(nlfpsamp - len(lfpchannel)) * np.nan))
		lfp = np.vstack((lfp, lfpchannel))

		spktime = trial['CH'][0][cc]['latency']
		spktime = np.int32(spktime*1000)
		if (spktime>nbins-1).any():
			rm_ix = spktime>=nbins
			print 'Spike time(s) too late!'
			print spktime[rm_ix]
			spktime = spktime[~rm_ix]
			
		rast[tt, spktime] = 1
		if spktime.size > 0:
			spktrials = np.append(spktrials, np.ones(spktime.size) * tt)
			spkwaveform = np.concatenate((spkwaveform, trial['CH'][0][cc]['spkwaveform'].T), 0)
			
		# add to Epoch_Value
		stimID[tt, :] = thisstimID

			# end if valid ID

	# # stoopid TDT shift
	# stimID = stimID[1:, :]
	# rast = rast[:-1, :]
	# lfp = lfp[:-1, :]
	# spkwaveform = spkwaveform[spktrials<(rast.shape[0]-1)]


	remID = np.array([0., 0.])
	spk_mask, trial_mask = RF.make_spk_and_trial_masks(spktrials, stimID, remID)
	rast = rast[~trial_mask, :]
	# if spktime.size > 0:
	# 	spkwaveform = spkwaveform[~spk_mask]
	stimID = stimID[~trial_mask, :]
	
	# do number of spikes in raster equal number of spike waveforms saved?
	# assert rast.sum() == spkwaveform.shape[0]
	
	# save out to file
	u_.create_dataset('stimID', data = stimID)
	u_.create_dataset('chan', data = cc)
	u_.create_dataset('blockID', data = blockID)
	# add stimulus ID datasets to this stimset on this unit
	u_.create_dataset('lfp', data = lfp, compression = 'gzip')
	# u_.create_dataset('spkwaveform', data = spkwaveform, compression = 'gzip')
	u_.create_dataset('rast', data = rast, compression = 'gzip')
		
	if blockID.startswith('b'):
		rf = RF.calc_rf(rast, stimID)
		u_.create_dataset('rf', data = rf, compression = 'gzip')
	
	
def load_coords(session, v = True):
	
	try:
		coords = np.loadtxt(os.path.join(studydir, 'Sessions', session, 'experimentfiles', session + '.txt'), ndmin = 1)
		if v:
			print 'Found coordinates at %s' % os.path.join(studydir, 'Sessions', session, 'experimentfiles', session + '.txt')
	except:
		coords = np.nan
		if v:
			print 'Coordinates not found'
			
	return coords
	
def load_cfs(session, v = True):

	try:
		cfs = np.loadtxt(os.path.join(studydir, 'Sessions', session, 'cfs.txt'), 'float32', ndmin = 1)
		if v:
			print 'Found CFs at %s' % os.path.join(studydir, 'Sessions', session, 'cfs.txt')
	except:
		cfs = np.nan
		if v:
			print 'CFs not found'
	
	return cfs
	
def add_coords(u_, coords, unitnum):
	
	try:
		coord = coords[coords[:, 0] == unitnum, 1:3][0]
	except:
		coord = np.nan
	u_.create_dataset('coord', data = coord)

	
def add_cfs(u_, cfs, unitnum):
	
	try:
		ix = cfs[:, 0] == unitnum
		cf = cfs[ix, 1][0]
	except:
		cf = np.nan
	u_.create_dataset('cf', data = cf)

def remove_trials(spktimes, spktrials, spkwaveform, lfp, stimID, remID):
	
	spk_mask, trial_mask = RF.make_spk_and_trial_masks(spktrials, stimID, remID)
	
	if spk_mask.sum()>0:
		spktimes = spktimes[~spk_mask]
		spktrials = spktrials[~spk_mask]
		spkwaveform = spkwaveform[~spk_mask, :]
		lfp = lfp[~trial_mask, :]
		stimID = stimID[~trial_mask, :]
	
	return spktimes, spktrials, spkwaveform, lfp, stimID
	
def sort_redos(d, v = False):
	'''
	Input:
		d - session main directory, e.g. /Users/robert/Desktop/Fmr1_RR/Sessions/ko_nai_20120306
	'''
	# get a list of the files
	os.chdir(os.path.join(studydir, d, 'data'))
	files = glob.glob('[b|r]*.mat')
	files.append('^^^') # special treatment for the last file
	p = re.compile('_R[\d]*')
	for f in range(len(files)-1):
		blockid1 = files[f][:3]
		blockid2 = files[f+1][:3]
		if blockid1 == blockid2:
			if v:
				print 'Renamed %s to %s' % (files[f], '_' + files[f])
			os.rename(files[f], '_' + files[f])
		else:
			if files[f].find('_R') > -1:
				if v:
					print 'Renamed %s to %s' % (files[f], p.sub('', files[f]))
				os.rename(files[f], p.sub('', files[f]))
				
def unpack(f, params = ['rast', 'stimID']):
	val = []
	for param in params:
		val.append(f[param].value)
	return val

	
# def fix_redos(experiment, test = True):
# 	
# 	fpath = glob.glob(os.path.join(studydir, 'Sessions', experiment, 'data', '*.mat'))
# 	(absol, relat) = zip(*[os.path.split(i) for i in fpath])
# 	fname = [os.path.split(i)[1] for i in fpaths]
# 	blockname = [os.path.splitext(i)[0] for i in fnames]
# 	tmp = [i.split('_R') for i in blocknames]
# 	block = []
# 	repnum = []
# 	for t in tmp:
# 		block.append(t[0])
# 		if len(t)==1:
# 			repnum.append(0)
# 		elif len(t[1])==0:
# 			repnum.append(1)
# 		else:
# 			repnum.append(int(t[1]))
# 
# 	df = pd.DataFrame(dict(fpath = fpath, fname = fname, blockname = blockname, block = block, repnum = repnum))
# 	
# 	group = df.groupby(block)
# 	for label, g in group:
# 		print len(g)
# 		if len(g)==1:
# 			path0 = g.fpath.values[0]
# 			if path0.find('_R')>-1:
# 				path1 = os.path.join(studydir, 'Sessions', experiment, 'data', '%s.mat' % b[:b.find('_R])
# 				print 'Moving %s to %s' % (os.path.split(path0)[1], os.path.split(path1)[1])
# 				if not test:
# 					shutil.move(path0, path1)
# 	
# 	blocks = {blockblocks[0] : [blocknames[0]]}
# 	for blockname, blockblock in zip(blocknames[1:], blockblocks[1:]):
# 		if blocks.keys()[-1] == blockblock:
# 			blocks[blocks.keys()[-1]].append(blockname)
# 		else:
# 			blocks[blockblock] = [blockname]
# 	
# 	for key in blocks.iterkeys():
# 		block = blocks[key]
# 		if len(block)==1:
# 			b = block[0]
# 			if b.find('_R')>-1:
# 				path0 = os.path.join(studydir, 'Sessions', experiment, 'data', '%s.mat' % b)
# 				path1 = os.path.join(studydir, 'Sessions', experiment, 'data', '%s.mat' % b[:b.find('_R')])
# 				print 'Moving %s to %s' % (os.path.split(path0)[1], os.path.split(path1)[1])
# 				if not test:
# 					shutil.move(path0, path1)
# 		else:
# 			for b in block[:-1]:
# 				path0 = os.path.join(studydir, 'Sessions', experiment, 'data', '%s.mat' % b)
# 				path1 = os.path.join(studydir, 'Sessions', experiment, 'data', '_%s.mat' % b)
# 				print 'Moving %s to %s' % (os.path.split(path0)[1], os.path.split(path1)[1])
# 				if not test:
# 					shutil.move(path0, path1)
# 
# 				
# 			b = block[-1]
# 			path0 = os.path.join(studydir, 'Sessions', experiment, 'data', '%s.mat' % b)
# 			path1 = os.path.join(studydir, 'Sessions', experiment, 'data', '%s.mat' % b[:b.find('_R')])
# 			print 'Moving %s to %s' % (os.path.split(path0)[1], os.path.split(path1)[1])
# 			if not test:
# 				shutil.move(path0, path1)

	
	
	
	
	
	