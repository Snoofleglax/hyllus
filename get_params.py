import math as ma
import matplotlib.pyplot as py
from pylab import *
import scipy as sp
import numpy as np
import files as fi
import os

def get_params(x,y,z,vx,vy,vz):
	xsol=x+8
	rad2deg = 180.0 / np.pi
	rgal=np.sqrt(np.power(x,2)+np.power(y,2)+np.power(z,2))
	rsol=np.sqrt(np.power(xsol,2)+np.power(y,2)+np.power(z,2))
	vgal=np.sqrt(np.power(vx,2)+np.power(vy,2)+np.power(vz,2)) #incorrect (value as calculated can't be negative); revisit
	vsol=(xsol*vx+y*vy+z*vz)/rsol
	l=(np.arctan2(y,xsol))*rad2deg
	i = 0
	for i in range(len(x)):
	  if l[i] < 0:
	    l[i] = l[i] + 360.0
	b=(np.arctan2(z,sqrt(np.power(xsol,2)+np.power(y,2))))*rad2deg
	os.system("touch " + str(filename) + ".params")
	params = open(str(filename) + ".params",'r+')
	params.write('#x,y,z,vx,vy,vz,rgal,vgal,r,vgsr,l,b\n')
	j = 0
	for j in range(len(x)):
	  print >> params, x[j], y[j],  z[j], vx[j], vy[j], vz[j], rgal[j], vgal[j], rsol[j], vsol[j], l[j], b[j]
	params.close()

filename = raw_input("Specify name of data file: ")
if __name__ == "__main__": 	
	data1=fi.read_file(str(filename), ",")
	x=data1[:,0] #pulls out a single column (note: i=0,1,2,3...)
	y=data1[:,1]
	z=data1[:,2]
	vx=data1[:,3]
	vy=data1[:,4]
	vz=data1[:,5]
	get_params(x,y,z,vx,vy,vz)