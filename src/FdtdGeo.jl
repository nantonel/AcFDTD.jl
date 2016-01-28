function CreateGeometry(Nx,Ny,Nz)
	#=
        Create cuboid geometry
        inputs: Nx, number of samples x direction
        Ny, number of samples y direction
        Nz, number of samples z direction
        output: G, 3D tensor containing 1 inside domain and 0 outside
	=#

	G = ones(Int32,Nx,Ny,Nz)                 # creates the tensor with 1 where domain is present 

	G2 = zeros(Int32,Nx+2,Ny+2,Nz+2)# puts a shell of zeros around the tensor
	G2[2:Nx+1, 2:Ny+1, 2:Nz+1] = G
	
	return G2
end

function CreateGeometryLshaped(Nx,Ny,Nz)
	#=
        Create cuboid geometry
        inputs: Nx, number of samples x direction
        Ny, number of samples y direction
        Nz, number of samples z direction
        output: G, 3D tensor containing 1 inside domain and 0 outside
	=#

	G = ones(Int32,Nx,Ny,Nz)                 # creates the tensor with 1 where domain is present 
	G[end-round(Int32,Nx/2):end,end-round(Int32,Ny/3):end,:] = 0     # this change geometry from pure cuboid 

	G2 = zeros(Int32,Nx+2,Ny+2,Nz+2)# puts a shell of zeros around the tensor
	G2[2:Nx+1, 2:Ny+1, 2:Nz+1] = G
	
	return G2
end

function GetGeometryMatrix(G, a, b, λ2)
	#=
	Create matrix A which will be used in the FDTD
	Inputs: Nx,Ny,Nz, number of samples of the room
		a, b parameters of the FDTD scheme
		λ2 Courant number squared
	Output: A matrix, 
	        Gxx,Gyy, Gzz 4D tensors containing the ghost points
		             positions
			     the 4th dimension indicates positive/
			     negative direction over x,y,z
			     e.g. to see the ceiling surfaces
			     type Gzz[:,:,:,2]
	=#
	Nx,Ny,Nz = size(G)
	Gx = zeros(Int32,Nx,Ny,Nz,19)
	Gy = zeros(Int32,Nx,Ny,Nz,19)
	Gz = zeros(Int32,Nx,Ny,Nz,19)

	Gxy = zeros(Int32,Nx,Ny,Nz,19)
	Gxz = zeros(Int32,Nx,Ny,Nz,19) 
	Gyz = zeros(Int32,Nx,Ny,Nz,19)
	
	Gxyz = zeros(Int32,Nx,Ny,Nz,19)

	Gbcx = zeros(Int32,Nx,Ny,Nz)
	Gbcy = zeros(Int32,Nx,Ny,Nz)
	Gbcz = zeros(Int32,Nx,Ny,Nz)
	
	for l = 2:Nx-1, m = 2:Ny-1,i = 2:Nz-1
	if (G[l,m,i] == 1)
		# compute δ2x: Gx[:,:,:,1] diag 1, Gx[:,:,:,2] diag 0, Gx[:,:,:,3] diag +1
	 
		ax = [G[l-1,m,i],  G[l,m,i],  G[l+1,m,i]]
                Gx[l,m,i,[ 7,10,13]], Gbcx[l,m,i] = CheckAx(ax)

		# compute δ2y: Gy[:,:,:,1] diag -Nx, Gy[:,:,:,2] diag 0, Gy[:,:,:,3] diag +Nx
		ax = [G[l,m-1,i],  G[l,m,i],  G[l,m+1,i]]
                Gy[l,m,i,[ 9,10,11]], Gbcy[l,m,i] = CheckAx(ax)
		
		# compute δ2z: Gz[:,:,:,1] diag -Nx*Ny, Gz[:,:,:,2] diag 0, Gz[:,:,:,3] diag +Nx*Ny
		ax = [G[l,m,i-1],  G[l,m,i],  G[l,m,i+1]]
                Gz[l,m,i,[3,10,17]], Gbcz[l,m,i] = CheckAx(ax)
		
		diagxy =  [G[l+1,m+1,i],G[l+1,m-1,i],G[l-1,m+1,i],G[l-1,m-1,i]]
		Gxy[l,m,i,[14,13,12,11,10, 9, 8, 7, 6]] = CheckDiag(diagxy, Gbcx[l,m,i],Gbcy[l,m,i],Gbcz[l,m,i],l,m,i)

		diagxz =  [G[l+1,m,i+1],G[l+1,m,i-1],G[l-1,m,i+1],G[l-1,m,i-1]]
		Gxz[l,m,i,[19,13, 5,17,10, 3,15, 7, 1]] = CheckDiag(diagxz, Gbcx[l,m,i],Gbcz[l,m,i],Gbcz[l,m,i],l,m,i)

		diagyz =  [G[l,m+1,i+1],G[l,m+1,i-1],G[l,m-1,i+1],G[l,m-1,i-1]]
		Gyz[l,m,i,[18,11, 4,17,10, 3,16, 9, 2]] = CheckDiag(diagyz, Gbcy[l,m,i],Gbcz[l,m,i],Gbcz[l,m,i],l,m,i)


	end
	end

	#removing the shell of outer zeros
	G  =  G[2:end-1, 2:end-1, 2:end-1]	

	Gx = Gx[2:Nx-1,2:Ny-1,2:Nz-1,:]
	Gy = Gy[2:Nx-1,2:Ny-1,2:Nz-1,:]
	Gz = Gz[2:Nx-1,2:Ny-1,2:Nz-1,:]

	Gxy = Gxy[2:Nx-1,2:Ny-1,2:Nz-1,:]
	Gxz = Gxz[2:Nx-1,2:Ny-1,2:Nz-1,:]
	Gyz = Gyz[2:Nx-1,2:Ny-1,2:Nz-1,:]

	Gbcx = Gbcx[2:Nx-1,2:Ny-1,2:Nz-1,:]
	Gbcy = Gbcy[2:Nx-1,2:Ny-1,2:Nz-1,:]
	Gbcz = Gbcz[2:Nx-1,2:Ny-1,2:Nz-1,:]

	

	Nx,Ny,Nz = Nx-2,Ny-2,Nz-2 #getting back to original dimension

	L = Gx+Gy+Gz +a*(Gxy+Gxz+Gyz)        #Laplatian Tensor


	#here we create the laplatian matrix from the laplatian tensor
	Lm =      spdiagm(reshape(L[:,:,:, 1],Nx*Ny*Nz)[Nx*Ny+2:end]  ,        -Nx*Ny-1, Nx*Ny*Nz,Nx*Ny*Nz)
	#l-1, m  , i-1
	Lm = Lm + spdiagm(reshape(L[:,:,:, 2],Nx*Ny*Nz)[Nx*(Ny+1)+1:end]  ,  -Nx*(Ny+1), Nx*Ny*Nz,Nx*Ny*Nz)
	#l  , m-1, i-1
	Lm = Lm + spdiagm(reshape(L[:,:,:, 3],Nx*Ny*Nz)[Nx*Ny+1:end]      ,      -Nx*Ny, Nx*Ny*Nz,Nx*Ny*Nz)
	#l  , m  , i-1
	Lm = Lm + spdiagm(reshape(L[:,:,:, 4],Nx*Ny*Nz)[Nx*(Ny-1)+1:end]  ,  -Nx*(Ny-1), Nx*Ny*Nz,Nx*Ny*Nz)
	#l  , m+1, i-1
	Lm = Lm + spdiagm(reshape(L[:,:,:, 5],Nx*Ny*Nz)[Nx*Ny:end]        ,    -Nx*Ny+1, Nx*Ny*Nz,Nx*Ny*Nz)
	#l+1, m  , i-1

	Lm = Lm + spdiagm(reshape(L[:,:,:, 6],Nx*Ny*Nz)[Nx+2:end]         ,       -Nx-1, Nx*Ny*Nz,Nx*Ny*Nz)
	#l-1, m-1, i
	Lm = Lm + spdiagm(reshape(L[:,:,:, 7],Nx*Ny*Nz)[2:end]            ,          -1, Nx*Ny*Nz,Nx*Ny*Nz)
	#l-1, m  , i
	Lm = Lm + spdiagm(reshape(L[:,:,:, 8],Nx*Ny*Nz)[1:end-Nx+1]       ,        Nx-1, Nx*Ny*Nz,Nx*Ny*Nz)
	#l-1, m+1, i
	Lm = Lm + spdiagm(reshape(L[:,:,:, 9],Nx*Ny*Nz)[Nx+1:end]         ,         -Nx, Nx*Ny*Nz,Nx*Ny*Nz)
	#l  , m-1, i
	Lm = Lm + spdiagm(reshape(L[:,:,:,10],Nx*Ny*Nz)                   ,           0, Nx*Ny*Nz,Nx*Ny*Nz)
	#l  , m  , i
	Lm = Lm + spdiagm(reshape(L[:,:,:,11],Nx*Ny*Nz)[1:end-Nx]         ,          Nx, Nx*Ny*Nz,Nx*Ny*Nz)
	#l  , m+1, i
	Lm = Lm + spdiagm(reshape(L[:,:,:,12],Nx*Ny*Nz)[Nx:end]           ,       -Nx+1, Nx*Ny*Nz,Nx*Ny*Nz)
	#l+1, m-1, i
	Lm = Lm + spdiagm(reshape(L[:,:,:,13],Nx*Ny*Nz)[1:end-1]          ,           1, Nx*Ny*Nz,Nx*Ny*Nz)
	#l+1, m  , i
	Lm = Lm + spdiagm(reshape(L[:,:,:,14],Nx*Ny*Nz)[1:end-Nx-1]       ,        Nx+1, Nx*Ny*Nz,Nx*Ny*Nz)
	#l+1, m+1, i

	Lm = Lm + spdiagm(reshape(L[:,:,:,15],Nx*Ny*Nz)[1:end-Nx*Ny+1]    ,     Nx*Ny-1, Nx*Ny*Nz,Nx*Ny*Nz)
	#l-1, m  , i+1
	Lm = Lm + spdiagm(reshape(L[:,:,:,16],Nx*Ny*Nz)[1:end-Nx*(Ny-1)]  ,   Nx*(Ny-1), Nx*Ny*Nz,Nx*Ny*Nz)
	#l  , m-1, i+1
	Lm = Lm + spdiagm(reshape(L[:,:,:,17],Nx*Ny*Nz)[1:end-Nx*Ny]      ,      +Nx*Ny, Nx*Ny*Nz,Nx*Ny*Nz)
	#l  , m  , i+1
	Lm = Lm + spdiagm(reshape(L[:,:,:,18],Nx*Ny*Nz)[1:end-(Ny+1)*Nx]  ,   (Ny+1)*Nx, Nx*Ny*Nz,Nx*Ny*Nz)
	#l  , m+1, i+1
	Lm = Lm + spdiagm(reshape(L[:,:,:,19],Nx*Ny*Nz)[1:end-Nx*Ny-1]   ,      Nx*Ny+1, Nx*Ny*Nz,Nx*Ny*Nz)
	#l+1, m  , i+1


	L = λ2*L
	L[:,:,:,10] = L[:,:,:,10]+2.*G

	Gxx = zeros(Int32,Nx,Ny,Nz,2)  
	Gyy = zeros(Int32,Nx,Ny,Nz,2)
	Gzz = zeros(Int32,Nx,Ny,Nz,2)

	Gxx[:,:,:,1], Gxx[:,:,:,2] = -(Gx[:,:,:,7]-1).*G,-(Gx[:,:,:,13]-1).*G  #these contains where ghost points
	Gyy[:,:,:,1], Gyy[:,:,:,2] = -(Gy[:,:,:,9]-1).*G,-(Gy[:,:,:,11]-1).*G  #were substituted for each axis
	Gzz[:,:,:,1], Gzz[:,:,:,2] = -(Gz[:,:,:,3]-1).*G,-(Gz[:,:,:,17]-1).*G  #positive [:,1] and negative [:,2]

	return  λ2*Lm+2*spdiagm(reshape(G,Nx*Ny*Nz)), L,Gxx,Gyy,Gzz 


end


function CheckAx(ax)

	#check if any axial point is a ghost point
		if    (ax[1] == 1 && ax[2]  == 1 && ax[3] == 1 ) #air
			return [1, -2, 1], 0   # central finite difference, no boundary
		else
		if    (ax[1] == 0 && ax[3] == 1)                 # wall
			return [0, -1, 1], 1   # forward finite difference, boundary
		
		elseif(ax[1] == 1 && ax[3] == 0)                 # wall
			return [1, -1, 0], 1
		end
		end
end



function CheckDiag(diag,Gbcx,Gbcy, Gbcz, l,m,i)

	#check if any diagonal point is a ghost point
	sumdiag = sum(diag)
	if    (sumdiag == 4 || sumdiag == 0)  #most likely case
                return  [1, -2,  1,  
		        -2,  4, -2,
		  	 1, -2,  1]
	elseif(sumdiag != 0)
		# here  boundary conditions are enforced
		if    (sumdiag == 2) # wall

			if    (diag == [0,0,1,1])
                
				return  [0,  0,  0,  
		                        -1,  2, -1,
		  	                 1, -2,  1]
		
			elseif(diag == [0,1,0,1])

				return  [0, -1,  1,  
		                         0,  2, -2,
		  	                 0, -1,  1]

			elseif(diag == [1,1,0,0])

				return  [1, -2,  1,  
		                        -1,  2, -1,
		  	                 0,  0,  0]

			elseif(diag == [1,0,1,0])

				return  [1, -1,  0,  
		                        -2,  2,  0,
		  	                 1, -1,  0]
			end

		elseif(sumdiag == 1) #inner edge/corner

			if    (diag ==  [1,0,0,0])

				return  [1, -1,  0,  
		                        -1,  1,  0,
		  	                 0,  0,  0]

			elseif(diag ==  [0,1,0,0])

				return  [0, -1,  1,  
		                         0,  1, -1,
		  	                 0,  0,  0]

			elseif(diag ==  [0,0,1,0])

				return  [0,  0,  0,  
		                        -1,  1,  0,
		  	                 1, -1,  0]

			elseif(diag ==  [0,0,0,1])

				return  [0,  0,  0,  
		                         0,  1, -1,
		  	                 0, -1,  1]
			end

		elseif(sumdiag == 3) #outer edge/corner

			if (Gbcx==0 && Gbcy ==0)
		# here coherence conditions are enforced
	        		if     (diag== [0,1,1,1]) 
                
					return  [0, -1,  1,  
		                                -1,  3, -2,
		  	                         1, -2,  1]


				elseif (diag== [1,1,1,0])

					return  [1, -2,  1,  
		                                -2,  3, -1,
		  	                         1, -1,  0]
		 
	        		elseif (diag== [1,0,1,1]) 

					return  [1, -1,  0,  
		                                -2,  3, -1,
		  	                         1, -2,  1]

			        elseif (diag== [1,1,0,1])

					return  [1, -2,  1,  
		                                -1,  3, -2,
		  	                         0, -1,  1]
				end
		        else
				if    (Gbcx == 1)  #it's a wall!
					if    (diag == [1,1,0,1]) 
					
					       return  [1, -2,  1,  
		                                       -1,  2, -1,
		  	                                0,  0,  0]

					elseif (diag == [1,1,1,0])
				          
					       return  [1, -2,  1,  
		                                       -1,  2, -1,
		  	                                0,  0,  0]

				        elseif(diag == [0,1,1,1]) 

					       return  [0,  0,  0,  
		                                       -1,  2, -1,
		  	                                1, -2,  1]

					elseif(diag == [1,0,1,1])

					       return  [0,  0,  0,  
		                                       -1,  2, -1,
		  	                                1, -2,  1]

				        end
				elseif(Gbcy == 1) 
					if    (diag == [1,0,1,1]) 

					       return  [1, -2,  0,  
		                                       -2,  2,  0,
		  	                                1, -1,  0]

					elseif(diag == [1,1,1,0])
					       
					       return  [1, -1,  0,  
		                                       -2,  2,  0,
		  	                                1, -1,  0]

				        elseif(diag == [0,1,1,1]) 

					       return  [0, -1,  1,  
		                                        0,  2, -2,
		  	                                0, -1,  1]

					elseif(diag == [1,1,0,1])
				          
					       return  [0, -1,  1,  
		                                        0,  2, -2,
		  	                                0, -1,  1]
			        	end
				end
			end
		end
	end

end


