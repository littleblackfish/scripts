#! /usr/bin/env julia

# Code for calculating r^2 vs. time curves for 
# calculation of diffusion constant
# 2014, murozturk@ku.edu.tr

# function that smooths data by window averaging

function smooth(data::Vector, winSize=3)
	n = size(data,1)
	sdata = similar(data, n-winSize)
	for i=1:n-winSize
		sdata[i]=mean(data[i:i+winSize])
	end
	return sdata
end

# function that returns the index of center of the plectoneme
# this is defined as point of maximum curvature

function center(data::Matrix)
	nframes = size(data,2)
	c= Array(Int, nframes)

	for i=1:nframes
		c[i]=sortperm(smooth(data[:,i], 10) )[end]
	end
	return c
end

#function that returns r^2 vs t

function r2t(center::Vector)
	
	nframes=length(center)
	maxtau = int(nframes/4)
	maxtau=100
	rsqt = Array(Float64,maxtau)

	for tau=1:maxtau

		tmp = Array(Int, nframes-tau )
		
		for tzero=1:nframes-tau
			tmp[tzero] = (center[tzero]-center[tzero+tau])^2
		end
		rsqt[tau] = mean(tmp)
	end
	return rsqt
end

data=readdlm("curvature.dat");
rsqt=r2t(center(data))
writedlm("r2t.dat", rsqt)





