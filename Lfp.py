import numpy as np
import matplotlib.pylab as plt

def calc_power_spectrum(lfp):
	
	ntrials = lfp.shape[1]
	
	
	spec = np.zeros(ntrials)
	Fs = 384.
	x = np.zeros((129, 7, ntrials))
	for i in range(ntrials):
		x[:, :, i] = specgram(lfp[:127, i], Fs = Fs)[0]