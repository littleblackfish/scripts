#! /usr/bin/env python2

#twist and write integration for circular DNA

from MDAnalysis import *
from numpy import *
from scipy import interpolate
from mpl_toolkits.mplot3d import Axes3D
from pylab import *
from math import sqrt
from sys import argv,exit
import pyximport; pyximport.install()
from writhe import * 
from twr_linear import twist


def cross(a,b) : return a[1]*b[2]-a[2]*b[1] , a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]

def dot(a,b) : return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]

def norm(a)   : return sqrt(a[0]*a[0] +a[1]*a[1]+a[2]*a[2])
  
def resample(r,bins=10) :
  interpolated=[]
  
  n=len(r)  #number of actual samples
  
  r_extended = concatenate([r[-5:],r,r[:5]]) #extended boundary condition for circular systems
  
  x =arange( 1 -5 ,  n+1 +5 )		# n+10 points   
  x_int= arange( 1, n+1, 1./bins  ) 
  
  for d in r_extended.T :
    #tck = interpolate.splrep(range(len(r)), d )
    #y= interpolate.splev(x,tck)
    tck=interpolate.UnivariateSpline(x,d,k=3, s=0)
    y= tck (x_int )
    interpolated.append(y)
  
  return array(interpolated).T

def plot3d (r) :
	fig=figure(1)
	ax=fig.gca(projection='3d')
	
#	ax.plot(r.T[0],r.T[1],r.T[2])
#	print r.T[0]

#	r=resample(r,bins)
	ax.plot(r.T[0],r.T[1],r.T[2])
	ax.plot(r.T[0],r.T[1],r.T[2], '.')
#	ax.plot(r.T[0][:5],r.T[1][:5],r.T[2][:5],'.', lw=2, label='beginning')
#	ax.plot(r.T[0][-5:],r.T[1][-5:],r.T[2][-5:],'.', lw=2, label='end')
	show()

def derivative (r) :#tangent calculation with special attention to first and last tangent
	n=len(r)
	T=empty((n,3))
	
	for x in range(n) :
		T[x]=  (r[(x+1)%n]-r[x-1]) /2 
	
	return T

def tangent (r) :	#tangent calculation with special attention to first and last tangent
	n=len(r)
	T=empty((n,3))
	
	for x in range (n):
		dr=r[(x+1)%n]-r[x-1]
		T[x]= dr/norm(dr)
	
	return T

def writhe_py(T,r)  :
  n=len(T)
  ds=empty(n)
  for x in range (n):		#calculating ds 
	  ds[x] =  (norm (r[(x+1)%n] - r[x] ) + norm(r[x]-r[x-1]) ) /2
    
  writhe=0.
  for x in range(n) :
    print  "%i percent complete  \r " %(x*100./n) ,     
    
    for y in range(x+1,n) :   
	
	b=r[x] - r[y]
	c=sqrt(dot(b,b))

	if c != 0:	
		writhe += ds[x]*ds[y]*dot( cross(T[x],T[y]) , b ) / (c*c*c)

  return writhe/(2*pi)
 

if __name__ == '__main__': 
	if len(argv) == 2 : 
		system=Universe(argv[1])
		nframes=1
	elif len(argv) == 3 :  
		system=Universe(argv[1],argv[2])
		nframes=system.trajectory.numframes
	elif len(argv) == 4 :  
		system=Universe(argv[1],argv[2])
		nframes=int(argv[3])
	else : 
		print 'gimme a system, damnit!'
		exit()

	P1 = system.selectAtoms('name B1 or resname 1B') 
	P2 = system.selectAtoms('name B2 or resname 2B')

	bins=1 # number of interpolated points (1 to turn off)
	npoints=(len(P1)-1)*bins+1
	
	ltw,lwr_zero,lwr_prev=[],[],[]
	
	T_zero=array( [zeros(npoints), zeros(npoints), ones(npoints)] ).T
	T=T_zero

	print 'Calculating twist and writhe for %i frames with %i bins...' %(nframes,bins)
	for t in range (nframes):
		print '%.1f percent complete   \r' %(t*100./nframes),
		T_prev=T

		if bins==1 : 	# no interpolation
			r1=P1.coordinates()
			r2=P2.coordinates()
		else :		# with interpolation
			r1=resample(P1.coordinates(),bins)
			r2=resample(P2.coordinates(),bins)

		r=(r1+r2)/2.	#calculating the centerline
		u=r1-r2		#vectors between comp. bases.

		T=tangent(r)
		
		if bins ==1 : 	# no interpolation	
			ltw.append (twist(T,u))
			lwr_zero.append (writhe_inc (T_zero, T))
			lwr_prev.append (writhe_inc (T_prev, T))
		else :		# with interpolation	
			ltw.append (collapse(twist(T,u),bins))
			lwr_zero.append (collapse (writhe_inc (T_zero, T),bins))
			lwr_prev.append (collapse (writhe_inc (T_prev, T),bins))
		
		if t < nframes-1:  system.trajectory.next()
	print ''
	del system

	tw, wr_zero,wr_prev=[],[0],[0]
	
	for x in ltw : 		
		tw.append(sum(x))
		x-=mean(x)
	
	for x in lwr_zero :	
		wr_zero.append(sum(x))
		x-=mean(x)

	for x in lwr_prev :	wr_prev.append(wr_prev[-1]+sum(x))

	print 'saving stuff..'
	savetxt('twist.dat', array([ range (len(tw)), tw]).T )
	savetxt('twist-traj.dat', array(ltw).T)
	savetxt('writhe.dat', array([range(len(wr_zero)),wr_prev,wr_zero]).T)
	savetxt('writhe-traj-zero.dat', array(lwr_zero).T )
	savetxt('writhe-traj-prev.dat', array(lwr_prev).T )

