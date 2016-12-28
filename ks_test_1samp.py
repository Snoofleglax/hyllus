import scipy as sp
from scipy import stats
from scipy.stats import norm
import numpy as np
import files as fi

def ks1sample(in1):
  D2,p2 = sp.stats.kstest(in1, 'norm', args=(0, 120))
  print D2, p2
  
filename = raw_input("Name of data file: ")
column = int(raw_input("Column number of data: "))

if __name__ == "__main__": 	
	data1 = fi.read_file(str(filename), ",")
	in1 = data1[:,column]
	ks1sample(in1)