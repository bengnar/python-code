#!/usr/bin/env

import os
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.widgets import Button
from scipy.ndimage import imread
import wx


nr=[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,6.250000e-02,1.250000e-01,1.875000e-01,2.500000e-01,3.125000e-01,3.750000e-01,4.375000e-01,5.000000e-01,5.625000e-01,6.250000e-01,6.875000e-01,7.500000e-01,8.125000e-01,8.750000e-01,9.375000e-01,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,9.375000e-01,8.750000e-01,8.125000e-01,7.500000e-01,6.875000e-01,6.250000e-01,5.625000e-01,5.000000e-01,0.2]
ng=[0,0,0,0,0,0,0,0,6.250000e-02,1.250000e-01,1.875000e-01,2.500000e-01,3.125000e-01,3.750000e-01,4.375000e-01,5.000000e-01,5.625000e-01,6.250000e-01,6.875000e-01,7.500000e-01,8.125000e-01,8.750000e-01,9.375000e-01,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,9.375000e-01,8.750000e-01,8.125000e-01,7.500000e-01,6.875000e-01,6.250000e-01,5.625000e-01,5.000000e-01,4.375000e-01,3.750000e-01,3.125000e-01,2.500000e-01,1.875000e-01,1.250000e-01,6.250000e-02,0,0,0,0,0,0,0,0,0,0.2]
nb=[5.625000e-01,6.250000e-01,6.875000e-01,7.500000e-01,8.125000e-01,8.750000e-01,9.375000e-01,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,9.375000e-01,8.750000e-01,8.125000e-01,7.500000e-01,6.875000e-01,6.250000e-01,5.625000e-01,5.000000e-01,4.375000e-01,3.750000e-01,3.125000e-01,2.500000e-01,1.875000e-01,1.250000e-01,6.250000e-02,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.2]

ix2freq = np.int32(np.linspace(2, 64, len(nr)))

class Penetrations():
	
	def __init__(self, ax):
		
		self.pens = []
		self.cfs = []
		self.text = []
		self.npens = 0


		'''Get file path and display scopephoto'''
		# get scopephoto path
		dialog = wx.FileDialog(None, 'ScopePhoto:', style = wx.OK|wx.CANCEL)
		if dialog.ShowModal() == wx.ID_OK:
			impath = dialog.GetPath()
		dialog.Destroy()
		

		img = imread(impath) # load scopephoto
		ax.imshow(img) # show scopephoto
		
		# Set output directory for this experiment as the relative path of the scopephoto
		outputpath = os.path.split(impath)[0]

		'''Place the outputpath and note text on the axis'''
		xlim = np.array(ax.get_xlim())
		ylim = np.array(ax.get_ylim())
		
		colwidth = (xlim[1]-xlim[0]) / 20.
		rowwidth = (ylim[0]-ylim[1]) / 20.
		
		xloc = xlim[0] + np.abs(np.diff(xlim))/2.
		yloc = ylim[0] - rowwidth
		self.outputpath = ax.text(xloc, yloc, outputpath)
		
		yloc = ylim[0] - 2*rowwidth
		self.note = ax.text(xloc, yloc, 'Note')
		
		ax.figure.canvas.mpl_connect('button_press_event', self.button_press_callback)
		ax.figure.canvas.mpl_connect('motion_notify_event', self.motion_notify_callback)
			
		ax.callbacks.connect('xlim_changed', self.axis_update)
		ax.callbacks.connect('ylim_changed', self.axis_update)
		
		self.ax = ax
		self.selected = None
		self.textselected = None
		
	def button_press_callback(self, event):
		
		# if we are trying to move a penetration
		if event.key is None:
			if self.selected is None:
				selected = [p.contains(event)[0] for p in self.pens]				
				if sum(selected)==1:
					self.selected = np.int32(np.array(selected).nonzero()[0])
					return
				else: return
			else:
				self.selected = None
			
		# if we are trying to relabel a penetration
		elif event.key=='shift':
			selected = [p.contains(event)[0] for p in self.pens]
			if sum(selected)==1:
				self.textselected = np.int32(np.array(selected).nonzero()[0])
				self.edit_cf()
				return
			if self.note.contains(event)[0]:
				self.edit_note()
				return
			if self.outputpath.contains(event)[0]:
				self.edit_outputpath()
				return

		return
		
	def button_release_callback(self, event):
		self.selected = None
		return
		
	def motion_notify_callback(self, event):
		if event.inaxes!=self.ax: return
		if self.selected is None: return
		
		self.pens[self.selected].set_position((event.xdata, event.ydata))
		plt.draw()

		return

	def edit_cf(self):
		if self.textselected is None: return
		dialog = wx.TextEntryDialog(None, 'New CF: ', defaultValue = np.str(self.cfs[self.textselected]), style = wx.OK|wx.CANCEL)
		if dialog.ShowModal() == wx.ID_OK:
			try:
				freq = np.int32(dialog.GetValue())
			except ValueError: return
			clr_ix = np.abs(ix2freq-freq).argmin()
			clr = [nr[clr_ix], ng[clr_ix], nb[clr_ix]]
			self.cfs[self.textselected] = np.int32(dialog.GetValue())
			
		dialog.Destroy()
		
		self.pens[self.textselected].set_color(clr)
		plt.draw()
		self.textselected = None

	def edit_outputpath(self):
		
		dialog = wx.TextEntryDialog(None, 'Output dir: ', defaultValue = self.outputpath.get_text(), style = wx.OK|wx.CANCEL)
		if dialog.ShowModal() == wx.ID_OK:
			outputpath = dialog.GetValue()
			if os.path.exists(outputpath):
				self.outputpath.set_text(outputpath)
			else: return
		
		dialog.Destroy()
		
		plt.draw()

	def edit_note(self):
		
		dialog = wx.TextEntryDialog(None, 'Output dir: ', defaultValue = self.note.get_text(), style = wx.OK|wx.CANCEL)
		if dialog.ShowModal() == wx.ID_OK:
			self.note.set_text(dialog.GetValue())
			
		dialog.Destroy()
		
		plt.draw()


	def add_pen(self, event):

		# add number to the ScopePhoto
		# get current view
		xlim = self.ax.get_xlim()
		ylim = self.ax.get_ylim()
		colwidth = (xlim[1]-xlim[0]) / 20.
		rowwidth = (ylim[0]-ylim[1]) / 20.
		self.npens += 1
		x = xlim[0] + (rowwidth * ((self.npens%20)))
		y = ylim[1] + (colwidth * (np.floor(self.npens/20.)+1))
		self.pens.append(self.ax.text(x, y, np.str(self.npens)))
		self.cfs.append(0)
		
		# add text box to margin
		
		plt.draw()
		return


	def rem_pen(self, event):
		# Remove CF label box
		if self.npens==0: return
		# Remove movable text on ScopePhoto
		del self.pens[-1]
		del self.cfs[-1]
		self.npens -= 1
		plt.show()
		return
		
	def axis_update(self, event):

		xlim = np.array(self.ax.get_xlim())
		ylim = np.array(self.ax.get_ylim())
		
		colwidth = (xlim[1]-xlim[0]) / 20.
		rowwidth = (ylim[0]-ylim[1]) / 20.
		
		xloc = xlim[0] + np.abs(np.diff(xlim))/2.
		yloc = ylim[0] - rowwidth
		self.outputpath.set_position([xloc, yloc])
		
		yloc = ylim[0] - 2*rowwidth
		self.note.set_position([xloc, yloc])


	def save(self, event):
		
		if self.npens==0: return
		outputpath = self.outputpath.get_text()
		try:
			xy = np.empty((self.npens, 4), dtype = np.float32)
			for i in range(self.npens):
				xy_ = self.pens[i].get_position()
				xy[i, :] = np.array([i+1, xy_[0], xy_[1], self.cfs[i]])
			np.savetxt(os.path.join(outputpath, 'markpens.txt'), xy)
			ax.get_figure().savefig(os.path.join(outputpath, 'scopephotograph.png'))
			print xy

		except IOError:
			print '%s does not exist!' % outputpath

if __name__ == '__main__':
	
	fig = plt.figure(figsize = [7.175, 5.4])
	fig.subplots_adjust(bottom = 0)
	fig.subplots_adjust(left = 0)
	fig.subplots_adjust(right = 1)
	ax = fig.add_axes([0, 0.1, 1, 0.9])
	ax.set_xticks([])
	ax.set_yticks([])

	ax_plus = fig.add_axes([0.09, 0.02, 0.06, 0.04])
	ax_minus = fig.add_axes([0.02, 0.02, 0.06, 0.04])
	ax_save = fig.add_axes([0.9, 0.02, 0.06, 0.04])

	pens = Penetrations(ax)

	b_plus = Button(ax = ax_plus, label = '+')
	b_plus.on_clicked(pens.add_pen)
	b_minus = Button(ax = ax_minus, label = '-')
	b_minus.on_clicked(pens.rem_pen)
	b_save = Button(ax = ax_save, label = 'Save')
	b_save.on_clicked(pens.save)

	plt.show()