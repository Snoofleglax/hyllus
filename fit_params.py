import math as ma
import matplotlib.pyplot as py
from pylab import *
import scipy as sp
import numpy as np
import files as fi
import itertools
from subprocess import call

def fit_params(l,b,vgsr,r):
	db=list(itertools.repeat(1,len(l)))
	dvgsr=list(itertools.repeat(20,len(l)))
	dr=list(itertools.repeat(1,len(l)))
	call([" touch orbit.input "], shell=True)
	params=open('orbit.input','r+') #terminal command to create file: touch [filename]
	j=0
	for j in range(len(lgedit)):
		print >>params, l[j], b[j], db[j], vgsr[j], dvgsr[j], r[j], dr[j]

	params.close()

if __name__ == "__main__": 	
	data1=fi.read_file("Hermus_orbit_stars.csv", ",")
	l=data1[:,2] #pulls out a single column (note: i=0,1,2,3...)
	b=data1[:,3]
	r=data1[:,4]
	vgsr=data1[:,5]
	fit_params(l,b,r,vgsr)

