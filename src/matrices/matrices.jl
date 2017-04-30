export get_L, 
       get_Q, 
       get_B,
       get_Sel,
       get_Exp

"""
## Create FDTD Laplatian matrix
* usage: L = get_L(g::CuboidRoom)
"""
function get_L(g::CuboidRoom)
	
	Dx = spdiagm((-2*ones(g.Nx), ones(g.Nx-1), ones(g.Nx-1)),(0,1,-1),g.Nx,g.Nx )
	Dy = spdiagm((-2*ones(g.Ny), ones(g.Ny-1), ones(g.Ny-1)),(0,1,-1),g.Ny,g.Ny )
	Dz = spdiagm((-2*ones(g.Nz), ones(g.Nz-1), ones(g.Nz-1)),(0,1,-1),g.Nz,g.Nz )

	#perfect reflection in phase (forward diff)
	Dx[1,2],Dy[1,2],Dz[1,2] =  1, 1, 1                
	Dx[1,1],Dy[1,1],Dz[1,1] =  -1, -1, -1                
	Dx[end, end-1], Dy[end, end-1], Dz[end, end-1] =   1, 1, 1
	Dx[end, end], Dy[end, end], Dz[end, end] =   -1, -1, -1

	Dxx = kron(speye(g.Nz),kron(speye(g.Ny), Dx))
	Dyy = kron(speye(g.Nz),kron(Dy, speye(g.Nx)))
	Dzz = kron(Dz, kron(speye(g.Ny),speye(g.Nx)))

	L = Dxx+Dyy+Dzz+g.env.scheme.a*(Dxx*Dyy+Dxx*Dzz+Dyy*Dzz)+g.env.scheme.b*(Dxx*Dyy*Dzz)

	return L
end

"""
## Create FDTD loss matrix Q
* usage: Q = get_Q(g::CuboidRoom)
"""
function get_Q(g::CuboidRoom)
	
	qx = spdiagm([1/(g.ξ[1]*2);zeros(g.Nx-2);1/(g.ξ[2]*2)]) #factor 2 due to fwd diff
	qy = spdiagm([1/(g.ξ[3]*2);zeros(g.Ny-2);1/(g.ξ[4]*2)])
	qz = spdiagm([1/(g.ξ[5]*2);zeros(g.Nz-2);1/(g.ξ[6]*2)])

	Qx = kron(speye(g.Nz),kron(speye(g.Ny),qx))
	Qy = kron(speye(g.Nz),kron(qy,speye(g.Nx)))
	Qz = kron(qz,kron(speye(g.Ny),speye(g.Nx)))

	return Qx+Qy+Qz
end

"""
## Create FDTD matrix B
* usage: B = get_B(g::CuboidRoom, Nt::Int64)
"""
function get_B(g::CuboidRoom, Nt::Int64)

	L = get_L(g)
	Q = get_Q(g)

	λ2 = (g.env.scheme.λ)^2 
	Nxyz = size(L,1)

	Qm,A,Qp = get_QmAQp(g)
        
	B =(kron(speye(Nt+2),spdiagm(Qp)) 
	+kron(spdiagm(ones(Nt+2-1),-1,Nt+2,Nt+2),-A) 
	+kron(spdiagm(ones(Nt+2-2),-2,Nt+2,Nt+2),-spdiagm(Qm)))

	B[1:2*Nxyz,1:2*Nxyz] = speye(2*Nxyz) #initial conditions

	return B
	
end

"""
## Create selection matrix
* usage: F = get_Sel(g::CuboidRoom,Nt::Int,xr::Vector{Tuple{Int,Int,Int}})
"""
function get_Sel(g::CuboidRoom,Nt::Int, xr::Vector{Tuple{Int,Int,Int}})

	Ks = length(xr)   #Ks number of mics
		
	indxr = [sub2ind((g.Nx,g.Ny,g.Nz),x...) for x in xr]

	Nxyz = g.Nx*g.Ny*g.Nz
	 
	indposvec = zeros(Int64,Ks*Nt)
	for k =1:Ks
		for i = 1:Nt  
			indposvec[i+(k-1)*Nt] = indxr[k]
			indxr[k] += Nxyz
		end
	end
	indposvec = indposvec+2*Nxyz
	F =  sparse((1:Ks*Nt),indposvec,ones(Int64,Ks*Nt),Ks*Nt,Nxyz*(Nt+2))
	
	return F

end

"""
## Create selection matrix
* usage: E = get_Exp(g::CuboidRoom,Nt::Int64,pos::Vector{Tuple{Int,Int,Int}})
"""
function get_Exp(g::CuboidRoom,Nt::Int,pos::Vector{Tuple{Int,Int,Int}})

		
	Ks = length(pos)   #Ks number of sources
	
	indposm = [sub2ind((g.Nx,g.Ny,g.Nz) ,x...) for x in pos]

	Nxyz = g.Nx*g.Ny*g.Nz
	 
	indposvec = zeros(Int64,Ks*Nt)
	for k =1:Ks
		for i = 1:Nt  
			indposvec[i+(k-1)*(Nt)] = indposm[k]
			indposm[k] = indposm[k]+Nxyz
		end
	end
	indposvec = indposvec +2*Nxyz #remove ICs
	E =  sparse(indposvec,(1:Ks*Nt),ones(Int64,Ks*Nt),Nxyz*(Nt+2),Ks*Nt)
	
	return E


end

function get_QmAQp(geo::CuboidRoom)

	λ2 = (geo.env.scheme.λ)^2 
	Q = get_Q(geo)
	L = get_L(geo)

	Qm = geo.env.scheme.λ*diag(Q)-1 
	A  = λ2*L+2*speye(geo.Nx*geo.Ny*geo.Nz)
	Qp = geo.env.scheme.λ*diag(Q)+1 

	return Qm,A,Qp
end
