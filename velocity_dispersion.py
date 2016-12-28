import os
import scipy as sp
import numpy as np
import math as ma
import matplotlib.pyplot as py
import files as fi

#filename = raw_input("Enter filename: ")
#binsizeb = float(raw_input("Size of b bin: "))
#binsizel = float(raw_input("Size of l bin: "))
#binsb = np.arange(-90.0,90.0,binsizeb)
#binsl = np.arange(0,360,binsizel)
#vdisplist = []
#np.set_printoptions(threshold='nan')

def mkbins(l,b,vlos):
  binnedb = np.digitize(b,binsb) #creates array of bin numbers based on b
  binnedl = np.digitize(l,binsl) #creates array of bin numbers based on l
  binarrlist = zip(vlos,binnedl,binnedb) #creates N x 3 array from vlos and b and l bin numbers, so each vlos value is with its corresponding bin
  binarr = np.asarray(binarrlist)
  ilim = len(binsl)
  jlim = len(binsb)
  size = ilim * jlim
  i = 1
  j = 1
  vdisplist = []
  while i < ilim + 1:
    subsetl = binarr[binarr[:,1] == i]
    while j < jlim + 1:
      templist = []
      subsetb = subsetl[subsetl[:,2] == j]
      #if len(subsetb) == 0:
	#out = ('No points', str(i), str(j))
	#vdisplist.extend(out)
      #else:
      templist.append(subsetb[:,0])
      x = np.asarray([templist])
      y = np.std(x[0])
      linit = 0 + (i - binsizel)
      binit = -90 + (j - binsizeb)
      outtup = (y,linit,binit)
      out = np.asarray(outtup)
      vdisplist.extend(out)
      del templist
      j = j + 1
    j = 1
    i = i + 1
  vdisp = np.asarray(vdisplist)
  vdisp.shape = (size,3)
  vdisp[np.isnan(vdisp)]=-1
  os.system("touch " + str(filename) + ".vdisp")
  vdispout=open(str(filename) + ".vdisp",'r+')
  s = 0
  while s < size:
    print >>vdispout, vdisp[s,0], vdisp[s,1], vdisp[s,2]
    s = s + 1
  vdispout.close()
  #Z = np.zeros((jlim,ilim))
  #xidx = (vdisp[:,1]-binsl[0]).astype(int)
  #yidx = (vdisp[:,2]-binsb[0]).astype(int)
  #Z[xidx, yidx] = vdisp[:,0]

if __name__ == "__main__": 	
	data1 = fi.read_file("hermus.1e6.175pc.null.1.nbody.params")
	l = np.array(data1[:,10])
	b = np.array(data1[:,11])
	vlos = np.array(data1[:,9])
	print np.std(vlos)
	#mkbins(l,b,vlos)