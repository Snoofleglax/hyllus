#! /usr/bin/python
import os
from subprocess import call
import matplotlib.pyplot as plt
from matplotlib import cm
import math
from pylab import *
import scipy as sp
import numpy as np
import files as fi
#import files as fi
#This creates a histrogram and matches to another histogram.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
                #/# # # # # # # # # # # # # # \#
                #          Control Panel       #
                #\# # # # # # # # # # # # # # /#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
y = True
n = False


get_data_and_make_parameters = y
plot = y
# # # # # # # # # # #
#    plot settings  #
# # # # # # # # # # #
x_range = [-10, 10]
y_range = [-10, 10]
z_range = [-15, 15]
plot_title = 'Hermus Orbit'
# # # # # # # #


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
                #/# # # # # # # # # # # # # # \#
                #          Engine Room         #
                #\# # # # # # # # # # # # # # /#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
def plot_my_shit():
    f = open('pauls_rad_gnuplot_script.gnuplot', 'w')
    f.write("reset\n")
    f.write("set xlabel 'X'\n")
    f.write("set ylabel 'Y'\n")
    f.write("set zlabel 'Z'\n")
    f.write("set xzeroaxis\n")
    f.write("set yzeroaxis\n")
    f.write("set zzeroaxis\n")
    f.write("set xrange [" + str(x_range[0]) + ":" + str(x_range[1]) + "]\n")
    f.write("set yrange [" + str(y_range[0]) + ":" + str(y_range[1]) + "]\n")
    f.write("set zrange [" + str(z_range[0]) + ":" + str(z_range[1]) + "]\n")
    f.write("set size square \n")
    f.write("set title '" + plot_title + "'\n")
    f.write("set term wxt persist size 900,900 \n\n")
    f.write("set output 'orbit.png'\n")
    f.write("splot 'orbit.1.tab' every 1:2002 using 1:2:3 with points pointtype 7 pointsize 0.1, \
	     'hermus.stars.csv' using 1:2:3 with points pointtype 7 pointsize 1 \n")
    f.close()
    os.system("gnuplot pauls_rad_gnuplot_script.gnuplot")
    
    
    
def get_params2(l, b, r):
     dr = np.pi / 180.0
     sinl = np.sin(l * dr)
     cosl = np.cos(l * dr)
     sinb = np.sin(b * dr)
     cosb = np.cos(b * dr)
     x = r * cosl * cosb - 8
     y = r * sinl * cosb
     z = r * sinb
     os.system("rm hermus.stars.csv")
     os.system("touch hermus.stars.csv")
     params = open('hermus.stars.csv','r+') #terminal command to create file: touch [filename]
     j = 0
     for j in range(len(l)):
          print >>params, x[j], y[j],  z[j]

     params.close()

def get_data_and_parameters():
    data1 = fi.read_file("Hermus_orbit_stars.csv", ",")
    l = data1[:,2] #pulls out a single column (note: i=0,1,2,3...)
    b = data1[:,3]
    r = data1[:,4]
    get_params2(l,b,r)
    
def main():
        
    if(get_data_and_make_parameters == True):  
        get_data_and_parameters()

    if(plot == True):
        plot_my_shit()
main()
