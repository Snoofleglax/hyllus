import math as ma
import matplotlib.pyplot as py
from pylab import *
import scipy as sp
import numpy as np
import files as fi
import os


#These first two definitions are from Matt's git repository used to convert RA,Dec into Galactic l,b  When comparing his method to my method they are in agreement down to then thousandth of a degree 
#which is negligable. The l and b values in Hermus.csv are therefore correct and left alone.
arr = sp.array([0.0])
def angle_bounds2 (theta, phi):
	''' Sets two spherical angles in bounds, -90<theta<90; 0<phi<360'''
	if type(theta) != type(arr): theta = sp.array([theta])
	if type(phi) != type(arr): phi = sp.array([phi])
	for i in range(len(theta)):
		theta[i] = angle_bounds(theta[i], -180.0, 180.0)
	for i in range(len(theta)):
		if sp.fabs(theta[i]) > 90.0:
			theta[i] = 180.0 - theta[i]
			phi[i] = phi[i] + 180.0
	for i in range(len(theta)):
		theta[i] = angle_bounds(theta[i], -180.0, 180.0)
		phi[i] = angle_bounds(phi[i], 0.0, 360.0)
	for i in range(len(theta)):
		if (sp.fabs(theta[i]) == 90.0): phi[i] = 0.0
	if len(theta)==1:
		theta, phi = theta[0], phi[0]
	return theta, phi


def angle_bounds (angle, min=0.0, max=360.0):
	''' Keeps an angle, in degrees, in a 360 degree region'''
	while angle < min: angle = angle + 360.0
	while angle > max: angle = angle - 360.0
	return angle

def Disk_velocity(l1,b1,l2,b2):

	p=np.polyfit(l2, b2, 3)
	l=np.linspace(27,67,1000)
	#b=p[0]*l**5+p[1]*l**4+p[2]*l**3+p[3]*l**2+p[4]*l+p[5]
	b = p[0] * l ** 3 + p[1] * l ** 2 + p[2] * l + p[3] 
	print p

	i=0
	
	newl1=[]
	newb1=[]
	newra=[]
	newdec=[]
	newg0=[]
	newumg=[]
	newgmr=[]
	newlogg=[]
	newfeh=[]
	newvgsr=[]
	newdist=[]
	#Nbody cut
	for i in range(len(b1)):
		#if (np.abs(b1[i]-(p[0]*l1[i]**5+p[1]*l1[i]**4+p[2]*l1[i]**3+p[3]*l1[i]**2+p[4]*l1[i]+p[5]))<2):
		if (np.abs(b1[i]-(p[0] * l1[i] ** 3 + p[1] * l1[i] ** 2 + p[2] * l1[i] + p[3]))<3):
			if(l1[i]>27):
				if(l1[i]<67):
					newl1.append(l1[i])
					newb1.append(b1[i])
					newra.append(ra[i])
					newdec.append(dec[i])
					newg0.append(g0[i])
					newumg.append(umg[i])
					newgmr.append(gmr[i])
					newlogg.append(logg[i])
					newfeh.append(feh[i])
					newvgsr.append(vgsr[i])
					newdist.append(dist[i])


	os.system("touch northern_BHBs_3deg_orb")
	test1=open('northern_BHBs_3deg_orb','r+') #name of output file
	j=0
	for j in range(len(newl1)):
	  print >>test1, newra[j], newdec[j], newl1[j], newb1[j], newdist[j], newvgsr[j], newg0[j], newumg[j], newgmr[j], newlogg[j], newfeh[j]
	test1.close()


if __name__ == "__main__": 	

	data1=fi.read_file("northern_BHBs", ",") #name of input nbody file
	
	data2=fi.read_file("hyllustest3.combined.100M.27to67.csv", ",") #names of input orbit files
	
	#for Pauls nbody data and orbit params
	l1=data1[:,5]
	b1=data1[:,6]
	ra=data1[:,3]
	dec=data1[:,4]
	g0=data1[:,0]
	umg=data1[:,1]
	gmr=data1[:,2]
	logg=data1[:,16]
	feh=data1[:,19]
	vgsr=data1[:,23]
	dist=data1[:,25]
	l2=data2[:,10]
	b2=data2[:,11]


	Disk_velocity(l1,b1,l2,b2)


