import os, glob, h5py, re
import numpy as np
import RF; reload(RF)
import RR; reload(RR)
import Spikes; reload(Spikes)
import misc; reload(misc)

verbose = True

ix2freq = 1000 * 2**(np.arange(0, 64)/10.)

basedir = '/Volumes/BOB_SAGET/Fmr1_RR/'
savepath = os.path.join(basedir, 'Analysis', 'RRTFs.h5')
fsave = h5py.File(savepath, 'a')

sessnames = glob.glob(os.path.join(basedir, 'Sessions/full_window/good', '*'))

for sessname in sessnames:
	fsess = fsave.create_group(os.path.split(sessname)[-1])
	rrpaths = glob.glob(os.path.join(basedir, 'Sessions', sessname, 'fileconversion', 'RR*.h5'))
	cfs = np.loadtxt(os.path.join(basedir, 'Sessions', sessname, 'cfs.txt'))

	for rrpath in rrpaths:
	
		absol, relat = os.path.split(rrpath)
		p = re.compile('RR')
		rfpath_ = p.sub('RF', relat)
		rfpath = os.path.join(absol, rfpath_)
		blockname = os.path.splitext(relat)[0]
		
		funit = fsess.create_group(blockname)
		
		p = re.compile('(\d+).h5')
		unitnum = np.int32(p.findall(relat))[0] # unit number
		if verbose:
			print rrpath
			print rfpath

		# open the RR file, get rast and stimparams, then close it
		frr = h5py.File(rrpath, 'r')
		rast = frr['rast'].value
		stimparams = frr['stimID'].value
		frr.close()

		# if not np.isnan(spktimes[0]):
		cf = cfs[cfs[:, 0]==unitnum, 1][0]
	
		freqs = stimparams[:, 0]
		rrs = stimparams[:, 1]
		ufreqs = np.unique(freqs)
		urrs = np.unique(rrs)
		nrrs = urrs.size
	
		# now we determine which of the frequencies we played is closest to this neuron's CF
		thisfreq, thisfreq_ix, thisfreq_err = misc.closest(ufreqs, ix2freq[20:][int(cf)], log = True)
		if np.abs(thisfreq_err) > 0.2:
			print 'No close frequency found!'
		thisfreq = ufreqs[thisfreq_ix]
	
		# isolate the parts of the raster for this frequency and build a psth for each RR
		ix = RF.get_trials(stimparams, np.array([thisfreq, np.nan]))
		thisrast = rast[ix, :]
		thisstims = stimparams[ix, :]
		psths, ustims = Spikes.calc_psth_by_stim(thisrast, thisstims)
	
		rrtf = RR.calc_rrtf_all(thisrast, thisstims, thisfreq, urrs)
	
	
		funit.create_dataset('rrtf', data = rrtf)
		funit.create_dataset('freq', data = thisfreq)
		funit.create_dataset('rrs', data = urrs)
		funit.create_dataset('cf', data = cf)

fsave.close()