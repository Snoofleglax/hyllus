import os
import scipy as sp
import numpy as np
import math as ma
import files as fi
import pylab as pl


filename = raw_input("Enter filename: ")
lmin = float(raw_input("Min l value: "))
lmax = float(raw_input("Max l value: "))
bmin = float(raw_input("Min b value: "))
bmax = float(raw_input("Max b value: "))
lrange = lmax - lmin
brange = bmax - bmin

def vdispplot(x,y,z):
  lenx = int(np.amax(x) - np.amin(x))
  leny = int(np.amax(y) - np.amin(y))
  vdispmatrix = np.zeros((leny,lenx))
  i = 0
  j = 0
  while i < lenx:
    while j < leny:
      counter = (180 * i) + j
      vdispmatrix[j][i] = z[counter]
      j = j + 1
    i = i + 1
    j = 0
  data = vdispmatrix
  dx = np.arange(lmin,lmax,lrange/10)
  dy = np.arange(bmin + 90,bmax + 90,brange/10)
  pl.xticks(dx)
  pl.yticks(dy)
  pl.xlim(lmin,lmax)
  pl.ylim(bmin + 90,bmax + 90)
  pl.xlabel('l')
  pl.ylabel('b + 90')
  pl.pcolor(data, vmin=0,vmax=100)
  pl.colorbar()
  pl.show()

if __name__ == "__main__": 	
	data1 = fi.read_file(str(filename))
	x = np.asarray(data1[:,1])
	y = np.asarray(data1[:,2])
	z = np.asarray(data1[:,0])
	vdispplot(x,y,z)