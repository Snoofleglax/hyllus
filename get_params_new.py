import math as ma
import matplotlib.pyplot as py
from pylab import *
import scipy as sp
import numpy as np
import files as fi
import os

arr = sp.array([0.0])
deg = 180.0 / np.pi
i = 0
def get_params(x,y,z,vx,vy,vz):
	xsol=x+8
	rgal=np.sqrt(np.power(x,2)+np.power(y,2)+np.power(z,2))
	rsol=np.sqrt(np.power(xsol,2)+np.power(y,2)+np.power(z,2))
	vgal=np.sqrt(np.power(vx,2)+np.power(vy,2)+np.power(vz,2)) #incorrect (can't be negative); revisit
	vsol=(xsol*vx+y*vy+z*vz)/rsol
	l=(np.arctan2(y,xsol))*(180.0/np.pi)
	for i in range(len(x)):
	  if l[i] < 0:
	    l[i] = l[i] + 360.0
	b=(np.arctan2(z,sqrt(np.power(xsol,2)+np.power(y,2))))*(180.0/np.pi)
	x1,y1,z1 = xyz2plane(x,y,z, new_x=[23.156,2.270,-5.887], plane=[-0.064,0.970,0.233,0.232], origin=None)
	lamhc, betahc, rsaghc = xyz2longlat(x1,y1,z1)
	xsaggc = 0.998*x + 0.066*y
	ysaggc = -0.015*x + 0.232*y - 0.972*z
	os.system("touch " + str(filename) + ".params")
	params=open(str(filename) + ".params",'r+') #terminal command to create file:
	j=0
	for j in range(len(x)):
		print >>params, x[j], y[j],  z[j], vx[j], vy[j], vz[j], rgal[j], vgal[j], rsol[j], vsol[j], l[j], b[j], lamhc[j], betahc[j], rsaghc[j], xsaggc[j], ysaggc[j]

	params.close()
	
#def xyz2lambeta(x,y,z, new_x=[23.156,2.270,-5.887], plane=[-0.064,0.970,0.233,0.232], origin=None, verbose=1):
    #x1,y1,z1 = xyz2plane(x,y,z, new_x, plane, origin)
    #lam, beta, r = xyz2longlat(x1,y1,z1)
    #print len(lam),len(beta),len(r)
    #for i in range(len(lam)):
            #if lam[i] < 0.0:  lam[i] = lam[i] + 360.0
    ##else:
    ##    if lam < 0.0:  lam = lam + 360.0
    #print lam,beta,r
    #test1=open('Newby_output_2.csv','r+')
    #j=0
    #for j in range(len(lam)):
	    #print >>test1, x[j], y[j], z[j], lam[j], beta[j], r[j], g0[j]
    #test1.close()
    #return lam, beta, r

def xyz2plane(x,y,z, new_x=[], plane=[], origin=None):
    """ Converts galactic x,y,z into orbital plane x,y,z
        new_x is the x,y,z coordinates of new x-axis
        plane is a,b,c,d plane parameters: ax + by + cz + d = 0
        origin is x offset, use d_sun for sun-centered inputs"""
    # preliminary stuff
    if origin != None:  x = x - origin
    a,b,c,d = plane
    bottom = np.sqrt(a*a + b*b + c*c)  # normalize
    a,b,c,d = a/bottom, b/bottom, c/bottom, d/bottom
    px, py, pz = new_x
    bot = np.sqrt(px*px + py*py + pz*pz)  #normalize
    px, py, pz = px/bot, py/bot, pz/bot
    p0 = [px,py,pz]
    # do rotation
    z_hat = [a,b,c]
    y_hat = cross(z_hat, p0)
    x_hat = cross(y_hat, z_hat)
    if type(x)==type(arr) or type(x)==type([]):
        xp, yp, zp = [], [], []
        for i in range(len(x)):
            xp.append(dot([x[i],y[i],z[i]], x_hat))
            yp.append(dot([x[i],y[i],z[i]], y_hat))
            zp.append(dot([x[i],y[i],z[i]], z_hat))
    else:
        xp = dot([x,y,z], x_hat)
        yp = dot([x,y,z], y_hat)
        zp = dot([x,y,z], z_hat)
    return xp, yp, zp


def cross(a, b):
    """ returns the cross-product of two 3-vectors"""
    c1 = a[1]*b[2] - a[2]*b[1]
    c2 = a[2]*b[0] - a[0]*b[2]
    c3 = a[0]*b[1] - a[1]*b[0]
    return sp.array([c1,c2,c3])

def dot(a, b):
    """ returns the dot product of two 3-vectors"""
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]  

def xyz2longlat(x,y,z):
    """ converts cartesian x,y,z coordinates into spherical longitude and latitude """
    r=[]
    long=[]
    d=[]
    lat=[]
    longnew=[]
    latnew=[]
    i=0
    for i in range(len(x)):
            r.append(sp.sqrt(x[i]*x[i] + y[i]*y[i] + z[i]*z[i]))
     	    long.append(sp.arctan2(y[i],x[i]))
	    d.append(sp.sqrt(x[i]*x[i] + y[i]*y[i]))
	    lat.append(sp.arcsin(z[i]/r[i]))

    i=0
    for i in range(len(x)):
            longnew.append(long[i]*deg)
            latnew.append(lat[i]*deg)

    return longnew, latnew, r

filename = raw_input("Specify name of data file: ")
if __name__ == "__main__": 	
	data1=fi.read_file(str(filename) + ".csv", ",")
	x=data1[:,0] #pulls out a single column (note: i=0,1,2,3...)
	y=data1[:,1]
	z=data1[:,2]
	vx=data1[:,3]
	vy=data1[:,4]
	vz=data1[:,5]
	get_params(x,y,z,vx,vy,vz)