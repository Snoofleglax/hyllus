import os
import scipy as sp
import numpy as np
import math as ma
import files as fi

filename = raw_input("Enter filename: ")

def orbittime(rgal):
  for i, j in enumerate(rgal):
    if j > 15.9:
      print i

if __name__ == "__main__": 	
	data1 = fi.read_file(str(filename))
	rgal = np.asarray(data1[:,6])
	orbittime(rgal)