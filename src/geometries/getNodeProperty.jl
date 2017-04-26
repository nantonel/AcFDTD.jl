function getNodeProperty(G::AbstractArray{Bool,3})
	Nx,Ny,Nz = size(G,1),size(G,2),size(G,3)
	I = Array{UInt8}(Nx-2,Ny-2,Nz-2)
	for l = 2:Nx-1, m = 2:Ny-1, i = 2:Nz-1
		if G[l-1,m,i] == 1 && G[l,m-1,i] == 1 && G[l,m,i-1] == 1 && 
		   G[l+1,m,i] == 1 && G[l,m+1,i] == 1 && G[l,m,i+1] == 1
		   I[l-1,m-1,i-1] = 1


			#walls
		elseif G[l-1,m,i] == 0 && G[l,m-1,i] == 1 && G[l,m,i-1] == 1 && 
		       G[l+1,m,i] == 1 && G[l,m+1,i] == 1 && G[l,m,i+1] == 1
		       I[l-1,m-1,i-1] = 2 #wall l = 1
		elseif G[l-1,m,i] == 1 && G[l,m-1,i] == 1 && G[l,m,i-1] == 1 && 
		       G[l+1,m,i] == 0 && G[l,m+1,i] == 1 && G[l,m,i+1] == 1
		       I[l-1,m-1,i-1] = 3 #wall l = Nx
		elseif G[l-1,m,i] == 1 && G[l,m-1,i] == 0 && G[l,m,i-1] == 1 && 
		       G[l+1,m,i] == 1 && G[l,m+1,i] == 1 && G[l,m,i+1] == 1
		       I[l-1,m-1,i-1] = 4 #wall m = 1
		elseif G[l-1,m,i] == 1 && G[l,m-1,i] == 1 && G[l,m,i-1] == 1 && 
		       G[l+1,m,i] == 1 && G[l,m+1,i] == 0 && G[l,m,i+1] == 1
		       I[l-1,m-1,i-1] = 5 #wall m = Ny
		elseif G[l-1,m,i] == 1 && G[l,m-1,i] == 1 && G[l,m,i-1] == 0 && 
		       G[l+1,m,i] == 1 && G[l,m+1,i] == 1 && G[l,m,i+1] == 1
		       I[l-1,m-1,i-1] = 6 #wall i = 1
		elseif G[l-1,m,i] == 1 && G[l,m-1,i] == 1 && G[l,m,i-1] == 1 && 
		       G[l+1,m,i] == 1 && G[l,m+1,i] == 1 && G[l,m,i+1] == 0
		       I[l-1,m-1,i-1] = 7 #wall i = Nz
			

			#edges
		elseif G[l-1,m,i] == 0 && G[l,m-1,i] == 0 && G[l,m,i-1] == 1 && 
		       G[l+1,m,i] == 1 && G[l,m+1,i] == 1 && G[l,m,i+1] == 1
		       I[l-1,m-1,i-1] = 8 #edge l = 1, m = 1
		elseif G[l-1,m,i] == 0 && G[l,m-1,i] == 1 && G[l,m,i-1] == 1 && 
		       G[l+1,m,i] == 1 && G[l,m+1,i] == 0 && G[l,m,i+1] == 1
		       I[l-1,m-1,i-1] = 9 #edge l = 1, m = Ny
		elseif G[l-1,m,i] == 0 && G[l,m-1,i] == 1 && G[l,m,i-1] == 0 && 
		       G[l+1,m,i] == 1 && G[l,m+1,i] == 1 && G[l,m,i+1] == 1
		       I[l-1,m-1,i-1] = 10 #edge l = 1, i = 1
		elseif G[l-1,m,i] == 0 && G[l,m-1,i] == 1 && G[l,m,i-1] == 1 && 
		       G[l+1,m,i] == 1 && G[l,m+1,i] == 1 && G[l,m,i+1] == 0
		       I[l-1,m-1,i-1] = 11 #edge l = 1, i = Nz

		elseif G[l-1,m,i] == 1 && G[l,m-1,i] == 0 && G[l,m,i-1] == 1 && 
		       G[l+1,m,i] == 0 && G[l,m+1,i] == 1 && G[l,m,i+1] == 1
		       I[l-1,m-1,i-1] = 12 #edge l = Nx, m = 1
		elseif G[l-1,m,i] == 1 && G[l,m-1,i] == 1 && G[l,m,i-1] == 1 && 
		       G[l+1,m,i] == 0 && G[l,m+1,i] == 0 && G[l,m,i+1] == 1
		       I[l-1,m-1,i-1] = 13 #edge l = Nx, m = Ny
		elseif G[l-1,m,i] == 1 && G[l,m-1,i] == 1 && G[l,m,i-1] == 0 && 
		       G[l+1,m,i] == 0 && G[l,m+1,i] == 1 && G[l,m,i+1] == 1
		       I[l-1,m-1,i-1] = 14 #edge l = Nx, i = 1
		elseif G[l-1,m,i] == 1 && G[l,m-1,i] == 1 && G[l,m,i-1] == 1 && 
		       G[l+1,m,i] == 0 && G[l,m+1,i] == 1 && G[l,m,i+1] == 0
		       I[l-1,m-1,i-1] = 15 #edge l = Nx, i = Nz

		elseif G[l-1,m,i] == 1 && G[l,m-1,i] == 0 && G[l,m,i-1] == 0 && 
		       G[l+1,m,i] == 1 && G[l,m+1,i] == 1 && G[l,m,i+1] == 1
		       I[l-1,m-1,i-1] = 16 #edge m = 1, i = 1
		elseif G[l-1,m,i] == 1 && G[l,m-1,i] == 0 && G[l,m,i-1] == 1 && 
		       G[l+1,m,i] == 1 && G[l,m+1,i] == 1 && G[l,m,i+1] == 0
		       I[l-1,m-1,i-1] = 17 #edge m = 1, i = Nz
		elseif G[l-1,m,i] == 1 && G[l,m-1,i] == 1 && G[l,m,i-1] == 0 && 
		       G[l+1,m,i] == 1 && G[l,m+1,i] == 0 && G[l,m,i+1] == 1
		       I[l-1,m-1,i-1] = 18 #edge m = Ny, i = 1
		elseif G[l-1,m,i] == 1 && G[l,m-1,i] == 1 && G[l,m,i-1] == 1 && 
		       G[l+1,m,i] == 1 && G[l,m+1,i] == 0 && G[l,m,i+1] == 0
		       I[l-1,m-1,i-1] = 19 #edge m = Ny, i = Nz


			#corners
		elseif G[l-1,m,i] == 0 && G[l,m-1,i] == 0 && G[l,m,i-1] == 0 && 
		       G[l+1,m,i] == 1 && G[l,m+1,i] == 1 && G[l,m,i+1] == 1
		       I[l-1,m-1,i-1] = 20 #edge l = 1, m = 1, i = 1
		elseif G[l-1,m,i] == 0 && G[l,m-1,i] == 0 && G[l,m,i-1] == 1 && 
		       G[l+1,m,i] == 1 && G[l,m+1,i] == 1 && G[l,m,i+1] == 0
		       I[l-1,m-1,i-1] = 21 #edge l = 1, m = 1, i = Nz
		elseif G[l-1,m,i] == 0 && G[l,m-1,i] == 1 && G[l,m,i-1] == 0 && 
		       G[l+1,m,i] == 1 && G[l,m+1,i] == 0 && G[l,m,i+1] == 1
		       I[l-1,m-1,i-1] = 22 #edge l = 1, m = Ny, i = 1
		elseif G[l-1,m,i] == 0 && G[l,m-1,i] == 1 && G[l,m,i-1] == 1 && 
		       G[l+1,m,i] == 1 && G[l,m+1,i] == 0 && G[l,m,i+1] == 0
		       I[l-1,m-1,i-1] = 23 #edge l = 1, m = Ny, i = Nz

		elseif G[l-1,m,i] == 1 && G[l,m-1,i] == 0 && G[l,m,i-1] == 0 && 
		       G[l+1,m,i] == 0 && G[l,m+1,i] == 1 && G[l,m,i+1] == 1
			I[l-1,m-1,i-1] = 24 #edge l = Nx, m = 1, i = 1
		elseif G[l-1,m,i] == 1 && G[l,m-1,i] == 0 && G[l,m,i-1] == 1 && 
		       G[l+1,m,i] == 0 && G[l,m+1,i] == 1 && G[l,m,i+1] == 0
			I[l-1,m-1,i-1] = 25 #edge l = Nx, m = 1, i = Nz
		elseif G[l-1,m,i] == 1 && G[l,m-1,i] == 1 && G[l,m,i-1] == 0 && 
		       G[l+1,m,i] == 0 && G[l,m+1,i] == 0 && G[l,m,i+1] == 1
			I[l-1,m-1,i-1] = 26 #edge l = Nx, m = Ny, i = 1
		elseif G[l-1,m,i] == 1 && G[l,m-1,i] == 1 && G[l,m,i-1] == 1 && 
		       G[l+1,m,i] == 0 && G[l,m+1,i] == 0 && G[l,m,i+1] == 0
			I[l-1,m-1,i-1] = 27 #edge l = Nx, m = Ny, i = Nz
		end
	end
	return I
end
