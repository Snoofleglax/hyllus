import numpy as np
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
import files as fi

histdata = raw_input("Enter data file: ")
histcol = int(raw_input("Enter column: "))
binsize = float(raw_input("Bin size? "))

def hist(tohist):
  #data.shape = (len(data))
  print tohist.shape
  print tohist
  binl = np.divide(1,binsize)
  n, bins, patches = plt.hist(tohist, binl)
  plt.xlabel('b deviation')
  plt.ylabel('N')
  plt.axis([-5, 5, 0, 200])
  plt.show()

if __name__ == "__main__": 	
	datafile = fi.read_file(str(histdata))
	tohist = np.asarray(datafile)[:,histcol] #for Charles's polyfit files
	hist(tohist)