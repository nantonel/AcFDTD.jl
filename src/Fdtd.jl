function Fdtd(Nx,Ny,Nz,Nt,A,Qp,Qm,s,pos,posm)

	#=
	Finite Difference Time Domain simulator for 
	frequency indepenent boundaries 
	
	Inputs: Nx,Ny,Nz Dimensions (in samples) of 
	                 the axis where the geometry
			 is defined
		Nt       Number of Time Samples	 
		A        Geometry Matrix
		Qp,Qm    Loss matrices
		s        Source vector
		pos      Index of position of source 
		posm     Indices of position of 
		         K Microphones

	Outputs: p_out   array K×Nt containing the 
	                 simulated sound pressure
	=#

	p0 = zeros(Nx*Ny*Nz) #initialize sound pressure p^n-1
	p1 = zeros(Nx*Ny*Nz) #                          p^n
	p2 = zeros(Nx*Ny*Nz) #                          p^n+1


	if(size(pos,1) ==  3)
		Ks = size(pos,2)   #Ks number of sources
		indpos = zeros(Int32,Ks)
		for k = 1:Ks
		indpos[k] = sub2ind((Nx,Ny,Nz),pos[1,k],pos[2,k],pos[3,k])
                end
        end

	if(size(posm,1) ==  3)
		K = size(posm,2)   #K number of mics
		indposm = zeros(Int32,K)
		for k = 1:K
		indposm[k] = sub2ind((Nx,Ny,Nz),posm[1,k],posm[2,k],posm[3,k])
                end
        end

	p_out = zeros(Nt,K) # initialize output

	for n = 2:Nt+1
	
		p2 = (A*p1+Qm.*p0).*Qp  #compute next time step sound pressure
		
		p2[indpos,:] = p2[indpos,:] + s[n-1,:]'
		# inject source

		p0 = copy(p1)  # update sound pressures
		p1 = copy(p2)
		p_out[n-1,:] = copy(p2[indposm])'  #record microphones
	end

	return p_out

end

