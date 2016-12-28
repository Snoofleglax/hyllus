import files as fi
import numpy as np
import scipy as sc
import matplotlib.pyplot as py
import matplotlib.mlab as mlab
import scipy.stats as stats
from pylab import *

def Deviation(rdev,bdev,vdev):
	py.subplots_adjust(hspace=0.5, wspace=0.5)
	ax1=py.subplot(221)
	ax2=py.subplot(222)
	ax3=py.subplot(223)

	ax1.hist(rdev,bins=10,histtype='step', range=(-2.5,2.5))
	ax2.hist(bdev,bins=10,histtype='step', range=(-2.5,2.5))
	ax3.hist(vdev,bins=14,histtype='step', range=(-70,70))
	
	ax1.axis([-2.5,2.5,0,5])
	ax2.axis([-2.5,2.5,0,5])
	ax3.axis([-70,70,0,5])
	#nbins = len(ax2.get_xticklabels())
	#setp(ax2.yaxis.get_ticklabels(), visible=False)
	#setp(ax3.xaxis.get_ticklabels(), visible=False)
	#setp(ax4.xaxis.get_ticklabels(), visible=False)
	#setp(ax4.yaxis.get_ticklabels(), visible=False)
	#ax1.xaxis.set_major_locator(MaxNLocator(nbins=nbins+1, prune='upper'))
	#ax3.xaxis.set_major_locator(MaxNLocator(nbins=nbins+1, prune='upper'))
	ax1.set_xlabel('r deviation', fontsize=16)
	ax2.set_xlabel('b deviation', fontsize=16)
	ax3.set_xlabel('vgsr deviation', fontsize=16)
	#ax3.set_ylabel('Heliocentric Distance (kpc)', fontsize=16)
	py.show()

if __name__ == "__main__": 
	
	data1=fi.read_file("rdev", ",")
	data2=fi.read_file("bdev", ",")
	data3=fi.read_file("vdev", ",")

	rdev=data1[:,0]
	bdev=data2[:,0]
	vdev=data3[:,0]
	
	Deviation(rdev,bdev,vdev)
