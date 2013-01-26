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
import os, glob, h5py, re
import Spikes; reload(Spikes)
import RF; reload(RF)


prefix = {'b' : 'RF', 'r' : 'RR', 'q' : 'QU', 'v' : 'VOC', 'i' : 'IS', 'P' : 'RF'}

electrodetype = 'tungsten'
electrodeinfo = {'fourbyfour' : {'nchan' : 16, 'npen' : 4, 'chanorder' : np.array([[1, 7, 13, 14], [3, 4, 10, 16], [2, 8, 12, 11], [6, 5, 9, 15]]), 'npen' : 4}, 'tungsten' : {'nchan' : 4, 'npen' : 4, 'chanorder' : None}, 'pen' : {'nchan' : 1, 'npen' : 1, 'chanorder' : None}}

stimulusinfo = {'RF' : {'nbins' : 333}, 'RR' : {'nbins' : 6000}, 'VOC' : {'nbins' : 20000}}

nchanperblockdict = {'b' : 4, 'i' : 16}
# basedir = '/Users/robert/Desktop/tetrode_test'
# basedir = '/Volumes/BOB_SAGET/Fmr1_RR/Sessions'
# basedir = '/Volumes/BOB_SAGET/for_shaowen/'
# basedir = '/Volumes/BOB_SAGET/Fmr1_KO_ising/Sessions/good/'
basedir = '/Volumes/BOB_SAGET/Fmr1_voc/'

def fileconversion(experiments, v = True):
	
	if type(experiments) is str:
		experiments = [experiments]
	elif experiments is None:
		experiments = glob.glob(os.path.join(basedir, '*'))
		
	# compile regexp patterns
	p_blocknum = re.compile('(?<=[a-zA-Z])[\d]+')
	
	for experiment in experiments:
		print experiment

		#get [b|r|etc]##.mat file names
		files = glob.glob(os.path.join(basedir, experiment, 'data', '[Ptbsyrpvi]*.mat'))
		nfiles = len(files)
		
		# make the destination directory if it doesn't exist (~/fileconversion/)
		if not os.path.exists(os.path.join(basedir, experiment, 'fileconversion')):
			os.mkdir(os.path.join(basedir, experiment, 'fileconversion'))

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
		cfs = load_cfs(experiment, v = v)
		coords = load_coords(experiment, v = v)
		
		#loop through penetrations
		for blocknum, blockrep, blockID in blockinfo:

			fileconversion_block(experiment, blocknum, blockID, blockrep, coords = coords, cfs = cfs, v = v, electrodeinfo = electrodeinfo[electrodetype], stimulusinfo = stimulusinfo[prefix[blockID[0]]])

		# end block loop

	# end experiment loop

def fileconversion_block(experiment, blocknum, blockID, blockrep, coords = None, cfs = None, v = True, electrodeinfo = electrodeinfo['tungsten'], stimulusinfo = stimulusinfo['RF']):
	
	chanorder = electrodeinfo['chanorder']
	nchan = electrodeinfo['nchan']
	npen = electrodeinfo['npen']
	
	if cfs is None:
		cfs = load_cfs(experiment)
	if coords is None:
		coords = load_coords(experiment)
	
	tempfile = scipy.io.loadmat(os.path.join(basedir, experiment, 'data', blockID + '.mat')) # load the .mat file0
	Data0 = tempfile['Data'] # save the Data file0 to a temporary variable
	# Data0 = tempfile['data'] # when converting already separated penetrations

	nchan = Data0['trial'][0][0][0][0]['CH'][0].size

	print '\t%s' % blockID

	blockstr = '%s%3.3i' % (prefix[blockID[0]], blocknum)
	if blockrep is not None:
		blockstr = '%s_%1.1i' % (blockstr, blockrep)
	print '\t%s' % blockstr

	savepath = os.path.join(basedir, experiment, 'fileconversion', blockstr + '.h5')
	b_ = h5py.File(savepath, 'w')

	for cc in range(nchan):
	
		u_ = b_.create_group('site%2.2i' % cc)
		unit(u_, Data0, cc, blockID, nbins = stimulusinfo['nbins'])

		if chanorder is None:
			siteno = ((blocknum-1)*npen) + cc + 1
		elif len(chanorder.shape) == 2:
			siteno =((blocknum-1)*npen) + (chanorder == cc+1).nonzero()[1][0]+1
			print cc+1, siteno
			

		add_coords(u_, coords, siteno)
		add_cfs(u_, cfs, siteno)

	b_.close()

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
			# spktimes = np.append(spktimes, spktime)
			spktrials = np.append(spktrials, np.ones(spktime.size) * tt)
			spkwaveform = np.concatenate((spkwaveform, trial['CH'][0][cc]['spkwaveform'].T), 0)
			
		# add to Epoch_Value
		stimID[tt, :] = thisstimID

			# end if valid ID
	
	# do number of spikes in raster equal number of spike waveforms saved?
	assert rast.sum() == spkwaveform.shape[0]
	
	remID = np.array([0., 0.])
	spk_mask, trial_mask = RF.make_spk_and_trial_masks(spktrials, stimID, remID)
	rast = rast[~trial_mask, :]
	if spktime.size > 0:
		spkwaveform = spkwaveform[~spk_mask]
	stimID = stimID[~trial_mask, :]
	b_ = u_.parent
	if not 'stimID' in b_:
		b_.create_dataset('stimID', data = stimID)
	
	
	# save out to file
	u_.create_dataset('chan', data = cc)
	u_.create_dataset('blockID', data = blockID)
	# add stimulus ID datasets to this stimset on this unit
	u_.create_dataset('lfp', data = lfp, compression = 'gzip')
	u_.create_dataset('spkwaveform', data = spkwaveform, compression = 'gzip')
	u_.create_dataset('rast', data = rast, compression = 'gzip')
		
	if blockID.startswith('b'):
		rf = RF.calc_rf(rast, stimID)
		u_.create_dataset('rf', data = rf, compression = 'gzip')
	
def load_coords(experiment, v = True):
	
	try:
		coords = np.loadtxt(os.path.join(basedir, experiment, 'experimentfiles', experiment + '.txt'), ndmin = 1)
		if v:
			print 'Found coordinates at %s' % os.path.join(basedir, experiment, 'experimentfiles', experiment + '.txt')
	except:
		coords = np.nan
		if v:
			print 'Coordinates not found'
			
	return coords
	
def load_cfs(experiment, v = True):

	try:
		cfs = np.loadtxt(os.path.join(basedir, experiment, 'cfs.txt'), 'float32', ndmin = 1)
		if v:
			print 'Found CFs at %s' % os.path.join(basedir, experiment, 'cfs.txt')
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
	os.chdir(os.path.join(basedir, d, 'data'))
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