import math as ma
import matplotlib.pyplot as py
from pylab import *
import scipy as sp
import numpy as np
#import files as fi
import os
from subprocess import call
import itertools
import scipy.stats
from scipy import optimize

filename = raw_input("Specify file name for orbit: ")
nsteps = float(raw_input("How many timesteps? "))
fortime = float(raw_input("How far forward should the orbit run (billions of years)? "))
backtime = float(raw_input("How far back should the orbit run? (billions of years)? "))
steplengthfor = float(fortime / nsteps)
steplengthback = float(backtime / nsteps)
potname = raw_input("Enter potential file name: ")
potquery = str(raw_input("Do you want to use the standard potential parameters? y/n: "))
if potquery == 'y' or potquery == 'Y' or potquery:
  potvars = str("0,6.5,0.26,4.45865888E5,0.7,1.52954402E5,73,1.0,12.0")
  omega = 0.0
  a = 6.5
  b = 0.26
  miya_mass = 4.45865888E5
  plummer_rc = 0.7
  plummer_mass = 1.52954402E5
  vhalo = 73.0
  q = 1.0
  d = 12.0
if potquery == 'n' or potquery == 'N':
  omega = raw_input("Enter omega value: ")
  a = raw_input("Enter Miyamoto a value: ")
  b = raw_input("Enter Miyamoto b value: ")
  miya_mass = raw_input("Enter Miyamoto mass: ")
  plummer_rc = raw_input("Enter Plummer scale radius: ")
  plummer_mass = raw_input("Enter Plummer mass: ")
  vhalo = raw_input("Enter halo velocity: ")
  q = raw_input("Enter flattening parameter: ")
  d = raw_input("Enter scale thickness: ")

xinit = raw_input("Enter initial x: ")
yinit = raw_input("Enter initial y: ")
zinit = raw_input("Enter initial z: ")
vxinit = raw_input("Enter initial vx: ")
vyinit = raw_input("Enter initial vy: ")
vzinit = raw_input("Enter initial vz: ")

csv = raw_input("CSV on output file? t/f: ")

os.system("mkorbit out=" + str(filename) + " x=" + str(xinit) + " y=" + str(yinit) + " z=" + 
str(zinit) + " vx=" + str(vxinit) + " vy=" + str(vyinit) + " vz=" + str(vzinit) + " potname=" + 
str(potname) + " potpars=" + str(omega) + "," + str(a) + "," + str(b) + "," + str(miya_mass) + 
"," + str(plummer_rc) + "," + str(plummer_mass) + "," + str(vhalo) + "," + str(q) + "," + str(d))

os.system("orbint in=" + str(filename) + " out=" + str(filename) + ".forward." + str(fortime) + 
"B nsteps=" + str(nsteps) + " dt=" + str(steplengthfor) + " potname=" + str(potname) + " potpars=" 
+ str(omega) + "," + str(a) + "," + str(b) + "," + str(miya_mass) + 
"," + str(plummer_rc) + "," + str(plummer_mass) + "," + str(vhalo) + "," + str(q) + "," + str(d))

os.system("orbint in=" + str(filename) + " out=" + str(filename) + ".back." + str(backtime) +
"B nsteps=" + str(nsteps) + " dt=-" + str(steplengthback) + " potname=" + str(potname) + 
" potpars=" + str(omega) + "," + str(a) + "," + str(b) + "," + str(miya_mass) + 
"," + str(plummer_rc) + "," + str(plummer_mass) + "," + str(vhalo) + "," + str(q) + "," + str(d))

os.system("otos in=" + str(filename) + ".forward." + str(fortime) + "B out=" + str(filename) + 
".forward." + str(fortime) + "B.snap")

os.system("otos in=" + str(filename) + ".back." + str(backtime) + "B out=" + str(filename) + 
".back." + str(fortime) + "B.snap")

os.system("snapprint in=" + str(filename) + ".forward." + str(fortime) + "B.snap tab=" + 
str(filename) + ".forward." + str(fortime) + "B.tab csv=" + str(csv))

os.system("snapprint in=" + str(filename) + ".back." + str(fortime) + "B.snap tab=" + 
str(filename) + ".back." + str(fortime) + "B.tab csv=" + str(csv))


#os.system("orbint in=orbit1.temp out=orbit1.for.temp nsteps=" + str(orbsteps) + " dt=" + str(orbtimestep) + " potname=" + str(pot) + " potpars=" + str(potpars))
#os.system("orbint in=orbit1.temp out=orbit1.back.temp nsteps=" + str(orbsteps) + " dt=-" + str(orbtimestep) + " potname=" + str(pot) + " potpars=" + str(potpars))
#os.system("otos in=orbit1.for.temp out=orbit.1.for.snap.temp")
#os.system("otos in=orbit1.back.temp out=orbit.1.back.snap.temp")
#os.system("snapstack in1=orbit.1.for.snap.temp in2=orbit.1.back.snap.temp out=orbit.1.total.snap.temp")
#os.system("snapprint in=orbit.1.total.snap.temp tab=orbit.1.tab csv=t")
  
  