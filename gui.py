import os, glob
import numpy as np
import matplotlib.pyplot as plt

plt.ion()
plt.rcParams['keymap.save'] = ''
plt.rcParams['keymap.yscale'] = ''
plt.rcParams['keymap.grid'] = ''


class gui_categorize():
	
	def __init__(self, data, plt_fcn = None):
		'''
		data : a list with nfactors elements per entry
		plt_fcn : the function that will be used to display the data, must take an argument ax that indicates which axis the data will be plotted to
		
		'''

		self.data = data
		self.counter = 0
		self.plt_fcn = plt_fcn
		self.fig = plt.figure()
		self.ax = []
		self.fig.canvas.mpl_connect('key_press_event', self.keypress)
		
		nfactors = 1#len(data[0])
		lenfactors = []
		for i in range(nfactors):
			lenfactors.append(len(self.data[i]))
			self.ax.append(self.fig.add_subplot(1, nfactors, i+1))
			
		assert len(np.unique(lenfactors)==1)
		self.categories = np.empty(lenfactors[0])
		self.plot_next()
		
	def plot_next(self):
	
		for i, ax_ in enumerate(self.ax):
			print self.counter, i
			ax_.cla()
			if self.plt_fcn is None:
				ax_.plot(self.data[self.counter])
			else:
				self.plt_fcn(self.data[self.counter], ax = ax_)
		plt.show()
				
	def keypress(self, event):
		
		
		if event.key == 'a': #small forward
			self.counter -= 1
		elif event.key == 'd': #small back
			self.counter += 1
		
		elif event.key == 'j':
			self.categories[self.counter] = True
			self.counter += 1
		elif event.key == 'l':
			self.categories[self.counter] = False
			self.counter += 1
			
		if event.key in 'adjl':
			self.plot_next()
		

