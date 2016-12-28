import os
import scipy as sp
import numpy as np
import math as ma
import matplotlib.pyplot as py
import files as fi

filename = raw_input("Enter filename: ")
lmin = float(raw_input("Min l value: "))
lmax = float(raw_input("Max l value: "))
bmin = float(raw_input("Min b value: "))
bmax = float(raw_input("Max b value: "))
lrange = lmax - lmin
brange = bmax - bmin
avgdisplist = []

def disprangeavg(l,b,disp):
  vdisplist = zip(disp,l,b)
  vdisparr = np.asarray(vdisplist)
  print vdisparr.shape
  vdisparr = vdisparr[np.all(vdisparr[:,0] >= 0, axis = 1)]
  print vdisparr
  print vdisparr.shape
  #i = 0
  #j = 0
  #ilim = lmax - lmin
  #jlim = bmax - bmin
  #k = lmin
  #l = bmin
  #while i < ilim - 1:
    #subsetl = vdisparr[vdisparr[:,1] == k]
    #while j < jlim - 1:
      #subsetb = vdisparr[vdisparr[:,2] == l]
      #templist.append[subsetb[:,0]]
      #j = j + 1
      #l = l + 1
if __name__ == "__main__": 	
	data1 = fi.read_file(str(filename))
	l = np.asarray(data1[:,1])
	b = np.asarray(data1[:,2])
	disp = np.asarray(data1[:,0])
	disprangeavg(l,b,disp)