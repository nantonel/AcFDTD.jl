function GetLossMatrix(Gx,Gy,Gz, υ, λ)

	#=
	Creates the Loss Matrix Qm, Qp that the Fdtd simulations 
	need Loss Matrices contains boundary condition.
	This functions assignes costant frequency independent υ 
	admittance to all surfaces

	Inputs: Gx, Gy, Gz containing boundary points information 
	                   (from GetGeometryMatrix)
	        υ          scalar value of the impedance
		λ          Courant Number

	Outputs: Qm, Qp    Loss Matrices 
	=#

	Nx, Ny, Nz, = size(Gx)
	Gbc = sum(Gx,4)+sum(Gy,4)+sum(Gz,4)
	# Gbc will contain 3 in corners
	#                  2 in edges
	#                  1 in walls
	# that is the number of ghost 
	# points substituted at that 
	# position

	Qp = (1+(λ/2.*υ).*Gbc) #create Loss matrix
	                        #note: the factor 2
				#is due to forward diff.
	Qp = Qp.^(-1)
	Qm = ((λ/2.*υ).*Gbc-1)

	Qm = reshape( Qm, Nx*Ny*Nz ) 
	Qp = reshape( Qp, Nx*Ny*Nz )
	
	return Qp, Qm
end


function GetLossMatrixFaces(Gx,Gy,Gz, υ, λ)

	#=
	Creates the Loss Matrix Qm, Qp that the Fdtd simulations 
	need Loss Matrices contains boundary condition.
	This functions assignes different frequency independent υ 
	admittances to all surfaces depending on direction.
	E.g. in a cubic room all different surfaces will have a
	different impedance
	
	Inputs: Gx, Gy, Gz containing boundary points information 
	                   (from GetGeometryMatrix)
	        υ          vector of length 6 containing the 
		           admittances of each surface
		λ          Courant Number

	Outputs: Qm, Qp    Loss Matrices 
	=#

	Nx, Ny, Nz, = size(Gx)
	
	Qp = (1+(λ/2*(υ[1])).*Gx[:,:,:,1]+
	        (λ/2*(υ[2])).*Gx[:,:,:,2]+
	        (λ/2*(υ[3])).*Gy[:,:,:,1]+
	        (λ/2*(υ[4])).*Gy[:,:,:,2]+
	        (λ/2*(υ[5])).*Gz[:,:,:,1]+
	        (λ/2*(υ[6])).*Gz[:,:,:,2])
	Qm = Qp-2	
	Qp = Qp.^(-1)

	Qm = reshape( Qm, Nx*Ny*Nz )
	Qp = reshape( Qp, Nx*Ny*Nz )
	
	return Qp, Qm
end


