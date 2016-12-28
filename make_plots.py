import math as ma
import matplotlib.pyplot as py
from pylab import *
import scipy as sp
import numpy as np
import files as fi
import os

filename = raw_input("Specify file name of orbit to plot: ")
fortime = float(raw_input("How far forward did the orbit run (billions of years)? "))
backtime = float(raw_input("How far back did the orbit run? (billions of years)? "))

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