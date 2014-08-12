#! /usr/bin/env julia

# function to calculate autocorrelation of twist over time
# rows are beads and colums are time
# produce bound region as
#bound = data[xboundlow:xboundhigh,tboundlow:tboundhigh];

function correlation (data)
	
	#normalize wrt mean
	data -= mean(data);
	nbeads, nframes = size(data)

	maxtau=500

	c= Array(Float64, maxtau)

	for tau=1:maxtau
		tmp = Array(Float64, nbeads, nframes-tau )
	
		for tzero=1:nframes-tau
			tmp[:,tzero] = data[:,tzero].*data[:,tzero+tau] ; 
		end
	
		c[tau] = mean(tmp); 
	end
	return c
end

