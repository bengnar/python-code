from matplotlib.widgets import Button, Slider
import os, glob
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from voc import plot_spec, load_spec, zscore_spec, lowpass_spec

plt.ion()
plt.rcParams['keymap.save'] = ''
plt.rcParams['keymap.yscale'] = ''
plt.rcParams['keymap.grid'] = ''

global hits, counter, savepaths, fnames

counter = -1

basedir = '/Volumes/BOB_SAGET/Vocalization/'
cages = ['KO5', 'WT5', 'WT6']

fig = plt.figure(figsize = [15, 8])
ax1 = fig.add_axes([0, 0.2, 1, 0.75])
ax2 = fig.add_subplot(616)

fig.subplots_adjust(bottom = 0.05)
fig.subplots_adjust(left = 0)
fig.subplots_adjust(right = 1)

fnames = []
savepaths = []
for cage in cages:
	sessdirs = glob.glob(os.path.join(basedir, 'Cages', cage, 'P*'))
	for sessdir in sessdirs:
		sess = os.path.split(sessdir)[-1]
		for i in np.arange(1, 10, 1):
			relat = '%s_%s_%2.2d' % (cage, sess, i)
			new_fname = os.path.join(basedir, 'Cages', cage, sess, relat+'.h5')
			new_savepath = os.path.join(basedir, 'analysis', 'gui_count_output', relat+'.txt')
			if not os.path.exists(new_savepath) and os.path.exists(new_fname):
				fnames.append(new_fname)
				savepaths.append(new_savepath)
				
fnames.sort()
savepaths.sort()
for (fname, savepath) in zip(fnames, savepaths):
	print fname, savepath

def keypress(event):

	global hits, counter, savepaths, fnames

	old_xlim = ax1.get_xlim()
	new_xlim = old_xlim
	xrange = old_xlim[1] - old_xlim[0]
	xcenter = np.mean(old_xlim)
	if event.key == 'd': #small forward
		ax1.set_xlim([old_xlim[0]+(xrange/4.), old_xlim[1]+(xrange/4.)])
	elif event.key == 's': #small back
		ax1.set_xlim([old_xlim[0]-(xrange/4.), old_xlim[1]-(xrange/4.)])
	elif event.key == 'f': #large forward
		ax1.set_xlim([old_xlim[1], old_xlim[1]+xrange])
	elif event.key == 'a': # large back
		ax1.set_xlim([old_xlim[0]-xrange, old_xlim[0]])

	elif event.key == 'e': #zoom in
		ax1.set_xlim([xcenter-0.25*xrange, xcenter+0.25*xrange])
	elif event.key == 'w': #zoom out
		ax1.set_xlim([xcenter-1.5*xrange, xcenter+1.5*xrange])

	elif event.key == 'g':
		xy = plt.ginput(n = 0, timeout = -1)
		hits = np.vstack((hits, xy))
		for xy_ in xy:
			ax1.axvline(xy_[0], color = 'r', ls = '--')
			ax2.axvline(xy_[0], color = 'r', ls = '--')
		plt.show()
		np.savetxt(savepaths[counter], hits)
		
	elif event.key == 'q':
		
		ax1.cla()
		ax2.cla()
		
		counter += 1
		
		print fnames[counter]
		print savepaths[counter]

		hits = np.empty((0, 2), dtype = 'float32')
		np.savetxt(savepaths[counter], hits)

		P, F, T = load_spec(fnames[counter])
		P, F = lowpass_spec(P, F)
		P = zscore_spec(P)
		P[P<0] = 0
		P = P**0.5

		wind = np.hamming(10)
		wind = wind / wind.sum()
		P_sum = (np.convolve(P.sum(0), wind, 'same'))

		nfreqs, ntimes = P.shape
		Tmax = T.max()

		plot_spec(P, F, T, ax = ax1)
		ax1.set_title(fnames[counter])
		ax1.set_xlim([0, 1.5])

		ax2.plot(T, P_sum)
		ax2.set_xlim([0, Tmax])
		
		# plt.show()

	plt.show()

fig.canvas.mpl_connect('key_press_event', keypress)

