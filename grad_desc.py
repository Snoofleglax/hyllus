import math as ma
import matplotlib.pyplot as py
from pylab import *
import scipy as sp
import numpy as np
import files as fi
import os
from subprocess import call
import itertools
import scipy.stats
from scipy import optimize
from scipy.optimize import differential_evolution

np.set_printoptions(threshold=np.nan)

#NEMO parameters
pot = "mpn"
potpars = "0,5.0,0.26,329751.696,0.33,85024.6529,114.0,1.0,12.0"
orbsteps = 5000
orbtimestep = 0.000016

#Chi-squared parameters
berror = 0.33
vgsrerror = 13.27
rerror = 0.66
nparams = 5

#global test 
#test = float(1e05)
#testlist = []
#optbound = ([-90,90],[0,+inf],None,None,None)
#Read in raw data, rearrange to useful form. This is obviously subject to 
#modification depending on the format of the input file, and perhaps later
#I can rewrite this bit to allow users to specify the input file and 
#columns instead of assuming a specific form.

bdev = open("bdev","w")
vdev = open("vdev","w")
rdev = open("rdev","w")

def fit_params(lraw,braw,rraw,vgsrraw):
	dbraw=list(itertools.repeat(1,len(lraw)))
	dvgsrraw=list(itertools.repeat(20,len(lraw)))
	drraw=list(itertools.repeat(1,len(lraw)))
	call([" touch orbit.input "], shell=True)
	inputfile=open('orbit.input','r+')
	j=0
	for j in range(len(lraw)):
		print >>inputfile, lraw[j], braw[j], dbraw[j], vgsrraw[j], dvgsrraw[j], rraw[j], drraw[j]

	inputfile.close()

if __name__ == "__main__": 	
	rawdata=fi.read_file("hyllus_pm2_6stars_092816.csv", ",")
	lraw=rawdata[:,2] #pulls out a single column (note: i=0,1,2,3...)
	braw=rawdata[:,3]
	rraw=rawdata[:,4]
	vgsrraw=rawdata[:,5]
	fit_params(lraw,braw,rraw,vgsrraw)

#Read input file in, get the data points for fitting

data = fi.read_file("orbit.input", " ")
l_in = data[:,0]
b_in = data[:,1]
db_in = data[:,2]
vgsr_in = data[:,3]
dvgsr_in = data[:,4]
r_in = data[:,5]
dr_in = data[:,6]
npoints = len(l_in)
ndata = 3*npoints



def MakeOrbit(x):
  global test
  
  b, r, vx, vy, vz = x
  deg2rad = np.pi/180.0
  rad2deg = 180.0/np.pi
  l = 42
  sinlinit = np.sin(l*deg2rad)
  coslinit = np.cos(l*deg2rad)
  sinbinit = np.sin(b*deg2rad)
  cosbinit = np.cos(b*deg2rad)
  xinit = r*coslinit*cosbinit - 8.0
  yinit = r*sinlinit*cosbinit
  zinit = r*sinbinit
  vxinit = vx
  vyinit = vy
  vzinit = vz
  
  #Call NEMO functions to make first orbit, integrate it back and forwards in 
  #time, convert the orbit files to snapshots, and stack the snapshots together 
  #before printing them out in ASCII form.
  
  os.system("mkorbit out=orbit1.temp x=" + str(xinit) + " y=" + str(yinit) + " z=" + str(zinit) + " vx=" + str(vxinit) + " vy=" + str(vyinit) + " vz=" + str(vzinit) + " potname=" + str(pot) + " potpars=" + str(potpars))
  os.system("orbint in=orbit1.temp out=orbit1.for.temp nsteps=" + str(orbsteps) + " dt=" + str(orbtimestep) + " potname=" + str(pot) + " potpars=" + str(potpars))
  os.system("orbint in=orbit1.temp out=orbit1.back.temp nsteps=" + str(orbsteps) + " dt=-" + str(orbtimestep) + " potname=" + str(pot) + " potpars=" + str(potpars))
  os.system("otos in=orbit1.for.temp out=orbit.1.for.snap.temp")
  os.system("otos in=orbit1.back.temp out=orbit.1.back.snap.temp")
  os.system("snapstack in1=orbit.1.for.snap.temp in2=orbit.1.back.snap.temp out=orbit.1.total.snap.temp")
  os.system("snapprint in=orbit.1.total.snap.temp tab=orbit.1.tab csv=t")

  #Remove temporary orbit files, leaving just the ASCII file

  os.system("rm *.temp") 

  #Get the output data file, convert coordinates, and write to a new file

  data1 = fi.read_file("orbit.1.tab", ",")
  x1 = data1[:,0]
  y1 = data1[:,1]
  z1 = data1[:,2]
  vx1 = data1[:,3]
  vy1 = data1[:,4]
  vz1 = data1[:,5]

  e = 0
  a = 0
  xsol = np.add(x1,float(8))
  l1rad = np.arctan2(y1,xsol)
  l1 = np.multiply(l1rad,rad2deg)
  for e in range(len(l1)):
    if (l1[e] < 0):
      l1[e] = l1[e] + 360
  r1x = np.power(xsol,2)
  r1y = np.power(y1,2)
  r1z = np.power(z1,2)
  r1sum = np.add(r1x,r1y)
  r1sqrt = np.add(r1sum,r1z)
  r1 = np.sqrt(r1sqrt)
  b1arg = np.divide(z1,r1)
  b1rad = np.arcsin(b1arg)
  b1 = np.multiply(b1rad,rad2deg)
  #for a in range(len(z1)):
    #if (btest[a] > 0):
      #print xsol[a]
      #print y1[a]
      #print z1[a]
      #print r1[a]
      #sys.exit("Bullshit in z1!")
  vgsr1x = np.multiply(xsol,vx1)
  vgsr1y = np.multiply(y1,vy1)
  vgsr1z = np.multiply(z1,vz1)
  vgsr1numpart = np.add(vgsr1x,vgsr1y)
  vgsr1num = np.add(vgsr1numpart,vgsr1z)
  vgsr1 = np.divide(vgsr1num, r1)
  
  db = list(itertools.repeat(berror,len(l1)))
  dvgsr = list(itertools.repeat(vgsrerror,len(l1)))
  dr = list(itertools.repeat(rerror,len(l1)))

  os.system("touch orbit.1.fit")

  params1=open('orbit.1.fit','r+')
  j=0
  for j in range(len(x1)):
		  print >>params1, l1[j], b1[j], db[j], vgsr1[j], dvgsr[j], r1[j], dr[j]

  params1.close()

  os.system("rm orbit.1.tab")

  obsdata = np.genfromtxt("orbit.input")
  testdata1 = np.genfromtxt("orbit.1.fit")

  # To do the chi-square test, we first need to search the orbit data for two 
  #points on opposite sides of each of the observed points. We'll do a 
  #single-variable search in l. First, we make the lists and indices.

  lindex = 0 
  bindex = 0
  vgsrindex = 0
  rindex = 0
  lsearch = []
  bsearch = []
  vgsrsearch = []
  rsearch = []

  #This searches through the test orbit by l and gets the two values of l that are 
  #closest to and on the opposite side of the data point.

  while lindex < npoints:
    lsub = np.subtract(testdata1[:,0],obsdata[lindex,0])
    lsubplus =  [m for m in lsub if m > 0]
    lsubminus = [n for n in lsub if n < 0]
    lminplus = np.amin(lsubplus)
    lminminus = np.amax(lsubminus)

    lminus = np.add(lminminus,obsdata[lindex,0])
    lplus = np.add(lminplus,obsdata[lindex,0])
    lsearch.append(lminus)
    lsearch.append(lplus)
    lindex = lindex + 1

  #Write the list of l values to an array
  larray = np.asarray(lsearch)

  #This searches through the test orbit for the b values corresponding to the 
  #closest l values from above.

  while bindex < 2 * npoints:
    bminustup = testdata1[np.where((testdata1[:,0] - larray[bindex]) == 0)]
    bplustup = testdata1[np.where((testdata1[:,0] - larray[bindex + 1]) == 0)]
    bplus = bplustup[0,1]
    bminus = bminustup[0,1]
    bsearch.append(bminus)
    bsearch.append(bplus)
    bindex = bindex + 2
    
  #Ditto, except for vgsr.

  while vgsrindex < 2 * npoints:
    vgsrminustup = testdata1[np.where((testdata1[:,0] - larray[vgsrindex]) == 0)]
    vgsrplustup = testdata1[np.where((testdata1[:,0] - larray[vgsrindex + 1]) == 0)]
    vgsrplus = vgsrplustup[0,3]
    vgsrminus = vgsrminustup[0,3]
    vgsrsearch.append(vgsrminus)
    vgsrsearch.append(vgsrplus)
    vgsrindex = vgsrindex + 2
    
  ##One more time, for r.

  while rindex < 2 * npoints:
    rplustup = testdata1[np.where((testdata1[:,0] - larray[rindex]) == 0)]
    rminustup = testdata1[np.where((testdata1[:,0] - larray[rindex + 1]) == 0)]
    rplus = rplustup[0,5]
    rminus = rminustup[0,5]
    rsearch.append(rminus)
    rsearch.append(rplus)
    rindex = rindex + 2
    
  os.system("rm orbit.1.fit")


  #Now write the fit parameter lists to an array:

  barray = np.asarray(bsearch)
  vgsrarray = np.asarray(vgsrsearch)
  rarray = np.asarray(rsearch)
  lindex2 = 0
  lindex3 = 0
  lindex4 = 0
  bindex2 = 0
  vgsrindex2 = 0
  rindex2 = 0
  lcoll = []
  lcoll = obsdata[:,0]
  lcol = np.array(lcoll)
  bmodellist = []
  vgsrmodellist = []
  rmodellist = []
  
  #We begin the chi-squared part. We first interpolate the orbit using the fit 
  #parameters.

  while bindex2 < 2 * npoints:
    bintnum = np.subtract(barray[bindex2+1],barray[bindex2])
    bintdenom = np.subtract(larray[bindex2+1],larray[bindex2])
    bintsub = np.subtract(lcoll[lindex2],larray[bindex2])
    bintquot = np.divide(bintnum,bintdenom)
    bintprod = np.multiply(bintquot,bintsub)
    bmodel = np.add(bintprod,barray[bindex2])
    bmodellist.append(bmodel)
    bindex2 = bindex2 + 2
    lindex2 = lindex2 + 1

  while vgsrindex2 < 2 * npoints:
    vgsrintnum = np.subtract(vgsrarray[vgsrindex2+1],vgsrarray[vgsrindex2])
    vgsrintdenom = np.subtract(larray[vgsrindex2+1],larray[vgsrindex2])
    vgsrintsub = np.subtract(lcol[lindex3],larray[vgsrindex2])
    vgsrintquot = np.divide(vgsrintnum,vgsrintdenom)
    vgsrintprod = np.multiply(vgsrintquot,vgsrintsub)
    vgsrmodel = np.add(vgsrintprod,vgsrarray[vgsrindex2])
    vgsrmodellist.append(vgsrmodel)
    vgsrindex2 = vgsrindex2 + 2
    lindex3 = lindex3 + 1
    
  while rindex2 < 2 * npoints:
    rintnum = np.subtract(rarray[rindex2+1],rarray[rindex2])
    rintdenom = np.subtract(larray[rindex2+1],larray[rindex2])
    rintsub = np.subtract(lcol[lindex4],larray[rindex2])
    rintquot = np.divide(rintnum,rintdenom)
    rintprod = np.multiply(rintquot,rintsub)
    rmodel = np.add(rintprod,rarray[rindex2])
    rmodellist.append(rmodel)
    rindex2 = rindex2 + 2
    lindex4 = lindex4 + 1

  bmodelarray = np.array(bmodellist)
  vgsrmodelarray = np.array(vgsrmodellist)
  rmodelarray = np.array(rmodellist)
  
  #Now, let's calculate the chi-squareds for each fit variable. Each calculated 
  #chi-squared value is fed into an array, which is then summed to get a final 
  #chi-square for the variable.

  chibindex = 0
  chivgsrindex = 0
  chirindex = 0

  chibind = np.empty([npoints,1])
  chivgsrind = np.empty([npoints,1])
  chirind = np.empty([npoints,1])

  while chibindex < npoints:
    chibnum = np.subtract(bmodelarray[chibindex],obsdata[chibindex,1])
    chibdenom = berror
    chibfrac = np.divide(chibnum,chibdenom)
    chibind[chibindex] = np.power(chibfrac,2)
    bdev.write(str(chibnum) + "\n")
    #bdev.write(str(obsdata[chibindex,0]) + "\n")
    chibindex = chibindex + 1
    
  chisqb = np.sum(chibind)
  
  while chivgsrindex < npoints:
    chivgsrnum = np.subtract(vgsrmodelarray[chivgsrindex],obsdata[chivgsrindex,3])
    chivgsrdenom = vgsrerror
    chivgsrfrac = np.divide(chivgsrnum,chivgsrdenom)
    chivgsrind[chivgsrindex] = np.power(chivgsrfrac,2)
    vdev.write(str(chivgsrnum) + "\n")
    chivgsrindex = chivgsrindex + 1
  
  chisqvgsr = np.sum(chivgsrind)
  
  while chirindex < npoints:
    chirnum = np.subtract(rmodelarray[chirindex],obsdata[chirindex,5])
    chirdenom = rerror
    chirfrac = np.divide(chirnum,chirdenom)
    chirind[chirindex] = np.power(chirfrac,2)
    rdev.write(str(chirnum) + "\n")
    chirindex = chirindex + 1

  chisqr = np.sum(chirind)

  #Now that we have the chi-squareds for the three fit variables, we combine them
  #into a goodness of fit value for the stream as a whole:

  eta = float(ndata - nparams - 1)
  chisqstrcoeff = float(1 / eta)
  chisqsum = chisqb + chisqvgsr + chisqr
  chisqstr = chisqstrcoeff * chisqsum
  print chisqstr
  return chisqstr

x0 = np.asarray((41.99,16.95,-138.83,125.53,14.48))
res1 = optimize.fmin_cg(MakeOrbit, x0, gtol = 1e-05, epsilon = 1e-02)
#MakeOrbit(x0)

bdev.close()
vdev.close()
rdev.close()
print res1

os.system("rm orbit.input")

