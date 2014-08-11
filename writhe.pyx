import numpy as np
from math import sqrt,pi

cdef inline double norm( a)   : return sqrt(a[0]*a[0] +a[1]*a[1]+a[2]*a[2])

def writhe_c(rin)  :
	cdef int x,y,n
	cdef double a[3],b[3],c,writhe
	n=len(rin)
	
	cdef double ds[100000]		#be careful with hard-coded limits. 
	cdef double  T[100000][3]	#implement dynamic allocation later on
	cdef double  r[100000][3]

	for x in range (n):		#calculating T and ds 
		dr = ( rin[(x+1)%n] - rin[x-1] )	#calculating normal tangent vectors 
		dr/= norm(dr)		#normalizing tangent vectors
		T[x][0]=dr[0]
		T[x][1]=dr[1]
		T[x][2]=dr[2]

		ds[x]= ( (norm (rin[(x+1)%n] - rin[x] ) + norm(rin[x]-rin[x-1]) ) /2)  	

	for x in range(n) : 		#copying numpy arrays to c arrays
		for y in range(3): 
			r[x][y]=rin[x][y]
	
	writhe=0
	for x in range(n) :		#main loop (expensive)
		for y in range(x+1,n) :   
			
			a[0]= T[x][1]*T[y][2]-T[x][2]*T[y][1]	#cross product
			a[1]= T[x][2]*T[y][0]-T[x][0]*T[y][2]	#a=T[x]xT[y]
			a[2]= T[x][0]*T[y][1]-T[x][1]*T[y][0]

			b[0]=r[x][0] - r[y][0]	#b=r[x]-r[y]
			b[1]=r[x][1] - r[y][1]
			b[2]=r[x][2] - r[y][2]

			c=sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]) #norm(r[x]-r[y])
			if c != 0 : writhe += ds[x]*ds[y]* (a[0]*b[0]+a[1]*b[1]+a[2]*b[2]) /(c*c*c)

	return writhe/(2*pi) #2pi instead of 4pi because of symmetric integral

