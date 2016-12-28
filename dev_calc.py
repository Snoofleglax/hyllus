import numpy as np
import scipy as sp
import mathplotlib.pyplot as py
import files as fi
import os

def fit_params(dev):
  i=0
  for i < len(dev):
    
	#dbraw=list(itertools.repeat(1,len(lraw)))
	#dvgsrraw=list(itertools.repeat(20,len(lraw)))
	#drraw=list(itertools.repeat(1,len(lraw)))
	#call([" touch orbit.input "], shell=True)
	#inputfile=open('orbit.input','r+')
	#j=0
	#for j in range(len(lraw)):
		#print >>inputfile, lraw[j], braw[j], dbraw[j], vgsrraw[j], dvgsrraw[j], rraw[j], drraw[j]

	#inputfile.close()

if __name__ == "__main__": 	
	rawdata=fi.read_file("bdev", ",")
	dev=rawdata[:,0] #pulls out a single column (note: i=0,1,2,3...)
	dev_calc(dev)