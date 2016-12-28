import math as ma
import numpy as np
import scipy as sc

def read_file(name, delimiter=None):  #Comment
	readfile = open(name, "r") 
	rows = []
	for line in readfile:
		if line.strip() == "": continue
		if line.strip()[0] == "#": continue
		rows.append(line)
	readfile.close()
	data = []
	for row in rows:
		holder = row.split(delimiter)
		for i in range(len(holder)):
			holder[i] = eval(holder[i].strip())
		data.append(holder)
	return sc.array(data)

if __name__ == "__main__":          #Asks if this is the main python script being called
    args = sys.argv[1:]             #Read Command line arguments, ignore first (this function)
    if len(args) == 1:              # Test to see if only the filename is given
        data = read_file(args[0])
    elif len(args) == 2:            # Test to see if both arguments are provided
        data = read_file(args[0], args[1])
    # If the wrong number of arguments is given, print out an error message and exit
    else:                           
        print "Invalid number of arguments, please give only two:\n  Filename (str), \
            Delimiter (str)"
        sys.exit(2)                 # Exit on internal error
    print data
