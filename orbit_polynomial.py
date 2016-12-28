import os
import scipy as sp
import numpy as np
import math as ma
import matplotlib.pyplot as py
import files as fi

filename = raw_input("Enter filename: ")
lcol = int(raw_input("Enter l column: "))
bcol = int(raw_input("Enter b column: "))

#lmin = float(raw_input("Min l value: "))
#lmax = float(raw_input("Max l value: "))
#bmin = float(raw_input("Min b value: "))
#bmax = float(raw_input("Max b value: "))

def orbit_polynomial(l,b):
  orbfit = np.polyfit(l, b, 3)
  print orbfit

if __name__ == "__main__": 	
	data1 = fi.read_file(str(filename), ",")
	l = np.array([data1[:,lcol]])
	b = np.array([data1[:,bcol]])
	orbit_polynomial(l,b)