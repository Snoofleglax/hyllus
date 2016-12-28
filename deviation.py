import numpy as np
import scipy as sp
import os
import files as fi

def get_devs(b,v,r):
  N = len(b)
  bsq = np.power(b,2)
  bquot = np.divide(bsq,N)
  bsum = np.sum(bquot)
  bd = np.sqrt(bsum)
  print str("b error = " + str(bd))
  vsq = np.power(v,2)
  vquot = np.divide(vsq,N)
  vsum = np.sum(vquot)
  vd = np.sqrt(vsum)
  print str("v error = " + str(vd))
  rsq = np.power(r,2)
  rquot = np.divide(rsq,N)
  rsum = np.sum(rquot)
  rd = np.sqrt(rsum)
  print str("r error = " + str(rd))
  

if __name__ == "__main__": 	
	data1 = fi.read_file("bdev")
	data2 = fi.read_file("vdev")
	data3 = fi.read_file("rdev")
	b = data1[:]
	v = data2[:]
	r = data3[:]
	get_devs(b,v,r)