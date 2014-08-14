#! /usr/bin/env julia

using PyCall
pyinitialize("python2") 

# function that locates plectonemes
# Nc is minimum plectoneme size in bases
# raise < actual raise (3.4) 
# plectoneme edges have i,j where i-j>Nc and dist<Nc*raise

function plocate(r::Matrix, Nc=40, raise=1.8)
	
	nbeads = size(r,1)
	d0=Nc*raise
	
	began,center, width=false, -1,0

	for i=1:nbeads-Nc		#start at the beginning
		for j=i+Nc:nbeads	#look at the rest of the chain
			dist = norm(r[i,:]-r[j,:])

			if !began && dist < d0
				#flag begining of plectoneme
				began=true

			elseif began && dist>d0
				#flag end of plectoneme
				center = (i+j)/2
				width  = j-i
				break
			end
		end
		if began
			break 
		end
	end

	return center,width
end

function plocate(mol::String, traj::String )

	@pyimport MDAnalysis

	system = MDAnalysis.Universe(mol, traj)

	P1 = system[:selectAtoms]("name B1 or resname 1B") 
	P2 = system[:selectAtoms]("name B2 or resname 2B")


	nframes=int(system[:trajectory][:numframes])
	nbeads = size(P1[:coordinates](),1)

	#plectoneme center
	pcenter = ones(nframes)*-1
	#plectoneme width
	pwidth     = zeros(nframes)

	#rewind just in case
	system[:trajectory][:rewind]()

	for t=1:nframes-1

		println (t, "\r")
		# skip the first frame (dummy.gro)
		system[:trajectory][:next]()
		
		# use centerline coordinates
		r= (P1[:coordinates]()+P2[:coordinates]())/2
		#print (r)

		pcenter[t], pwidth[t] = plocate(r)
	end
	return pcenter, pwidth
end


# function that returns deltar^2 vs. t assuming a continuous trajectory 

function r2t(center::Vector)
	
	nframes=length(center)
	maxtau = int(nframes/10)
	rsqt = Array(Float64,maxtau)

	for tau=1:maxtau

		tmp = Array(Float64, nframes-tau )
		
		for tzero=1:nframes-tau
			tmp[tzero] = (center[tzero]-center[tzero+tau])^2
		end
		rsqt[tau] = mean(tmp)
	end
	return rsqt
end

center, width = plocate("minimized.gro", "minimized.gro")
center, width = plocate("minimized.gro", "traj-sparse.xtc")

writedlm("center.dat", center)
writedlm("width.dat" , width)
writedlm("r2t.dat", r2t(center))
