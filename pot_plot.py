import math as ma
import matplotlib.pyplot as py
from pylab import *
import scipy as sp
import numpy as np

#k2m = float(3.086E+19)
#SM2kg = float(1.98855E+30)
#G = float(6.6741E-11)
#nfwa = float(25)
#nfwamet = np.multiply(nfwa,k2m)
#rho0 = float(1.840E-24)
#pi = float(np.pi)

G = float(1)
nfwa = float(25)
NFWcoeff = float(960.746)
logcoeff = 
s2struct = 222288.47

#a = float(raw_input("Enter Miyamoto a scale (kpc): "))
#b = float(raw_input("Enter Miyamoto b scale (kpc): "))

#rc = float(raw_input("Enter bulge radius (kpc): "))

#d = float(raw_input("Enter halo scale parameter (kpc): "))
#dsq = np.power(d,2)

#Dmass = float(raw_input("Enter disk mass (solar masses): "))
#Dmassstruct = np.divide(Dmass,s2struct)

#Bmass = float(raw_input("Enter bulge mass (solar masses): "))
#Bmassstruct = np.divide(Bmass,s2struct)

#vhalo = float(raw_input("Enter halo velocity (km/s): "))

#q = float(raw_input("Enter flattening parameter: "))

r = np.arange(0,100,0.1)
NFWfrac = np.divide(r,nfwa)
NFWlogarg = np.add(1,NFWfrac)
NFWlog = np.log(NFWlogarg)
NFW = np.divide(NFWlog,NFWfrac)
log
py.plot(r,NFW)
py.show()
#x = float(raw_input("Enter x (kpc): "))
#y = float(raw_input("Enter y (kpc): "))
#z = float(raw_input("Enter z (kpc): "))
#xmet = x*k2m
#ymet = y*k2m
#zmet = z*k2m
#rsq = np.power(xmet,2) + np.power(ymet,2) + np.power(zmet,2)
#rcylsq = np.power(xmet,2) + np.power(ymet,2)

#Diskpotnum = G*Dmassmet*float(-1)
#DiskpotdenomR = np.power(xmet,2)+np.power(ymet,2)
#Diskpotdenomsum = np.power(zmet,2)+np.power(bmet,2)
#Diskpotsqrt = np.sqrt(Diskpotdenomsum)
#Diskpotdenomsum2 = np.add(amet,Diskpotsqrt)
#Diskpotdenomsq = np.power(Diskpotdenomsum2,2)
#Diskpotdenomsum3 = np.add(DiskpotdenomR,Diskpotdenomsq)
#Diskpotdenom = np.sqrt(Diskpotdenomsum3)
#Diskpot = np.divide(Diskpotnum,Diskpotdenom)

#Bulgepotnum = float(-1)*G*Bmassmet
#Bulgepotdenomsqrt = np.sqrt(rsq)
#Bulgepotdenom = np.add(rcmet,Bulgepotdenomsqrt)
#Bulgepot = np.divide(Bulgepotnum,Bulgepotdenom)

#logcoeff = np.power(vhalo,2)
#lnargz = np.divide(zmet,q)
#lnargzsq = np.power(lnargz,2)
#lnargsum = np.add(rcylsq,lnargzsq)
#lnarg = np.add(lnargsum,dsq)
#log = np.log(lnarg)
#Logpot = np.multiply(logcoeff,log)

#nfwcoeff = -4*pi*G*rho0*nfwa
#logquot = np.divide(rsq,nfwa)
#logarg = np.add(1,logquot)
#nfwlog = np.log(logarg)
#nfwfrac = np.divide(nfwlog,logquot)
#NFWpot = np.multiply(nfwcoeff,nfwfrac)

#Potsumlog = Diskpot + Bulgepot + Logpot

#Potsumnfw = Diskpot + Bulgepot + NFWpot

#print Potsumlog
#print Potsumnfw



