#! /usr/bin/env python2

#Script to calculate local and total twist from a trajectory.
#It does not do any fancy processing, clean and raw. 
#Third argument is stride

from numpy import *
from MDAnalysis import *
from twr_linear import tangent, twist
from math import sqrt
from sys import argv,exit

print 'Loading trajectory...'

system=Universe(argv[1],argv[2])
nframes=system.trajectory.numframes

if len(argv)>3 :
    skip = int(argv[3])
else :
    skip = 1 

frames = range(1, nframes, skip)

P1 = system.selectAtoms('name B1 or resname 1B') 
P2 = system.selectAtoms('name B2 or resname 2B')

print 'Calculating twist for %i frames and %i beads..' %(len(frames), len(P1))

ltw = empty([len(frames), len(P1)])     #local twist
tw = empty(len(frames))                 #total twist

for i in range(len(frames)):
    
    system.trajectory[frames[i]] 

    print '%.1f %% complete   \r' %(i*100./len(frames)),
    r1=P1.coordinates()
    r2=P2.coordinates()
    
    r=(r1+r2)/2.	#calculating the centerline
    u=r1-r2		#vectors between comp. bases.

    T=tangent(r)
    ltw[i] = twist(T,u)
    tw[i]  = sum(ltw[i])

    

print '\nSaving files..'
savetxt('twist.dat', array([ frames, tw]).T )
savetxt('twist-traj.dat', array(ltw).T)


    


