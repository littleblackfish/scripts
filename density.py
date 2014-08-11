#! /usr/bin/env python2

from MDAnalysis import *
from numpy import *
from sys import argv
from pylab import *

b=200

def density (system) :
  all = system.selectAtoms('all')
#  coor=[]
  coor=  all.coordinates().T[2] 
  #coor/=system.dimensions[2]
  return histogram(coor, density=True, range=(0.,len(all)/4*3.2), bins=b)
  
def entropy (hist) :
	entropy=0
	for p in hist :	
		if p>0: entropy += -p*log(p)
	return entropy 

print "loading trajectory, this might take a while.."

system=Universe(argv[1], argv[2])
nframes=system.trajectory.numframes-1	
maxlength=system.dimensions[2]

stepsize =1 
density_timeline, ent=[],[]

print "calculating z-density for %i frames with %i stepsize and %i bins.." %(nframes, stepsize, b)

for x in range (0,nframes-stepsize, stepsize) :
  print  "%.1f percent complete  \r " %(x*100./nframes) ,
  for y in range (stepsize) : system.trajectory.next()
  h=density(system)
  h=array (h[0])
  density_timeline.append(h)
  ent.append(entropy(h))
  
savetxt( 'density.dat', array(density_timeline).T)
savetxt( 'entropy.dat', array([range(len(ent)),ent]).T)
