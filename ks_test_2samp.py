import scipy as sp
from scipy import stats
from scipy.stats import norm
import numpy as np
import files as fi

def ks2sample(in1,in2):
 D1,p1 = sp.stats.ks_2samp(in1,in2)
 print D1
 print p1

filename1 = raw_input("Name of first data file: ")
filename2 = raw_input("Name of second data file: ")
column1 = int(raw_input("Column number in first file: "))
column2 = int(raw_input("Column number in second file: "))

if __name__ == "__main__": 	
	data1 = fi.read_file(str(filename1), ",")
	data2 = fi.read_file(str(filename2), ",")
	in1 = data1[:,column1]
	in2 = data2[:,column2]
	ks2sample(in1,in2)