import math as ma
import matplotlib.pyplot as py
from pylab import *
import scipy as sp
import numpy as np
import files as fi


#These first two definitions are from Matt's git repository used to convert RA,Dec into Galactic l,b  When comparing his method to my method they are in agreement down to then thousandth of a degree which is negligable. The l and b values in Hermus.csv are therefore correct and left alone.
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

def Disk_velocity(l1,b1,x1,y1,z1,rgal1,rsol1,vgsr1,l2,b2,x2,y2,z2,rgal2,rsol2,vgsr2,l3,b3,x3,y3,z3,rgal3,rsol3,vgsr3):

	p=np.polyfit(l2, b2, 5)
	l=np.linspace(17,45,100)
	b=p[0]*l**5+p[1]*l**4+p[2]*l**3+p[3]*l**2+p[4]*l+p[5]

	q=np.polyfit(l3, b3, 5)
	ll=np.linspace(45,80,100)
	bb=q[0]*ll**5+q[1]*ll**4+q[2]*ll**3+q[3]*ll**2+q[4]*ll+q[5]	

	i=0
	
	newl1=[]
	newb1=[]
	newrsol1=[]
	newvgsr1=[]
	newrgal1=[]
	newx1=[]
	newy1=[]
	newz1=[]
	
	#Nbody cut
	for i in range(len(b1)):
		if (np.abs(b1[i]-(p[0]*l1[i]**5+p[1]*l1[i]**4+p[2]*l1[i]**3+p[3]*l1[i]**2+p[4]*l1[i]+p[5]))<2):
			if(l1[i]>17):
				if(l1[i]<45):
					newl1.append(l1[i])
					newb1.append(b1[i])
					newrsol1.append(rsol1[i])
					newvgsr1.append(vgsr1[i])
					newrgal1.append(rgal1[i])
					newx1.append(x1[i])
					newy1.append(y1[i])
					newz1.append(z1[i])


	i=0
	for i in range(len(b1)):
		if (np.abs(b1[i]-(q[0]*l1[i]**5+q[1]*l1[i]**4+q[2]*l1[i]**3+q[3]*l1[i]**2+q[4]*l1[i]+q[5]))<2):
			if(l1[i]>45):
				if(l1[i]<70):
					newl1.append(l1[i])
					newb1.append(b1[i])
					newrsol1.append(rsol1[i])
					newvgsr1.append(vgsr1[i])
					newrgal1.append(rgal1[i])
					newx1.append(x1[i])
					newy1.append(y1[i])
					newz1.append(z1[i])
	
	test1=open('sidd.5.05e6.globular.2.02B.polyfit.2deg','r+') #name of output file
	j=0
	for j in range(len(newl1)):
		print >>test1, newl1[j], newb1[j], newrgal1[j], newrsol1[j], newvgsr1[j], newx1[j], newy1[j], newz1[j]

	test1.close()
	py.scatter(newl1,newb1,c='r')
	py.plot(l3,b3)
	py.plot(l2,b2)
	py.show()

	
if __name__ == "__main__": 	

	data1=fi.read_file("sidd.5.05e6.globular.2.02B.params", " ") #name of input nbody file
	
	data2=fi.read_file("hermustest31.back.0.25B.params.small", ",") #names of input orbit files
	data3=fi.read_file("hermustest31.forward.0.25B.params.small", ",")
	
	#for Pauls nbody data and orbit params
	l1=data1[:,10]
	b1=data1[:,11]
	x1=data1[:,0]
	y1=data1[:,1]
	z1=data1[:,2]
	rgal1=data1[:,6]
	rsol1=data1[:,8]
	vgsr1=data1[:,9]

	l2=data2[:,10]
	b2=data2[:,11]
	x2=data2[:,0]
	y2=data2[:,1]
	z2=data2[:,2]
	rgal2=data2[:,6]
	rsol2=data2[:,8]
	vgsr2=data2[:,9]

	l3=data3[:,10]
	b3=data3[:,11]
	x3=data3[:,0]
	y3=data3[:,1]
	z3=data3[:,2]
	rgal3=data3[:,6]
	rsol3=data3[:,8]
	vgsr3=data3[:,9]
	
	Disk_velocity(l1,b1,x1,y1,z1,rgal1,rsol1,vgsr1,l2,b2,x2,y2,z2,rgal2,rsol2,vgsr2,l3,b3,x3,y3,z3,rgal3,rsol3,vgsr3)


