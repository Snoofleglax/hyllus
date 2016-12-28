import matplotlib.pyplot as py
from pylab import *
import scipy as sp
import numpy as np
import files as fi
import os

def get_params2(l,b,r):
	dr=np.pi/180
	sinl=np.sin(l*dr)
	cosl=np.cos(l*dr)
	sinb=np.sin(b*dr)
	cosb=np.cos(b*dr)
	x=r*cosl*cosb-8
	y=r*sinl*cosb
	z=r*sinb
	os.system("touch hyllus_pm3_5stars_xyz.csv")
	params=open('hyllus_pm3_5stars_xyz.csv','r+')
	j=0
	for j in range(len(l)):
		print >>params, x[j], y[j],  z[j]

	params.close()

if __name__ == "__main__": 	
	data1=fi.read_file("hyllus_pm3_5stars.csv", ",")
	l=data1[:,2] #pulls out a single column (note: i=0,1,2,3...)
	b=data1[:,3]
	r=data1[:,4]
	get_params2(l,b,r)

