export LShapedRoom

"""
`L-shaped cuboid room object`

Returns a `LShapedRoom` object containing the geometrical information of the room.

# Usage 

* LShapedRoom(Lx::Float64,Ly::Float64,Lz::Float64,xedge::Float64,yedge::Float64,ξ::Array{Float64,1}, env::FDTDEnv)
* LShapedRoom(Nx::Int64,Ny::Int64,Nz::Int64,xedge::Int64,yedge::Int64,ξ::Array{Float64,1}, env::FDTDEnv)
  * `ξ::Array{Float64}`: must be 6 element Array containg acoustic impedance of the 6 walls
  * `xedge` and `yedge` indicates the position of the edges, can also be expressed in samples

"""
immutable LShapedRoom <: AbstractGeometry
	Lx::Float64   #x dimension
	Ly::Float64   #y dimension
	Lz::Float64   #z dimension

	Nx::Int64   #x-axis samples
	Ny::Int64   #y-axis samples
	Nz::Int64   #z-axis samples

	xedge::Float64  #edge x position
	yedge::Float64  #edge y position

	G::Array{Bool,3}      #Geometry
	ξ::Array{Float64,1}   #Impedance
	env::FDTDEnv          #FDTD acoustic env.
end
	
function LShapedRoom(Lx::Float64,Ly::Float64,Lz::Float64,
		     xedge::Float64,yedge::Float64,
		    ξ::Array{Float64,1}, env::FDTDEnv)

	if(any([Lx;Ly;Lz].< 0)) error("room dimensions L should be positive") end
	if(length(ξ)!= 6) error("length(ξ) must be 6") end

	Nx = round(Int64,Lx/env.X)
	Ny = round(Int64,Ly/env.X)
	Nz = round(Int64,Lz/env.X)

	Ex = round(Int64,xedge/env.X)
	Ey = round(Int64,yedge/env.X)

	G = CreateGeometryLshaped(Nx,Ny,Nz,Ex,Ey)

	LShapedRoom(Lx,Ly,Lz,Nx,Ny,Nz,xedge,yedge,G,ξ,env)
end

function LShapedRoom(Nx::Int64,Ny::Int64,Nz::Int64,
		     Ex::Int64,Ey::Int64,
		    ξ::Array{Float64,1}, env::FDTDEnv)

	if(any([Nx;Ny;Nz].< 0)) error("room dimensions L should be positive") end
	if(length(ξ)!= 6) error("length(ξ) must be 6") end

	Lx =Nx*env.X
	Ly =Ny*env.X
	Lz =Nz*env.X

	xedge =Ex*env.X
	yedge =Ey*env.X
	G = CreateGeometryLshaped(Nx,Ny,Nz,Ex,Ey)

	LShapedRoom(Lx,Ly,Lz,Nx,Ny,Nz,xedge,yedge,G,ξ,env)
end

function CreateGeometryLshaped(Nx::Int64,Ny::Int64,Nz::Int64,
			       Ex::Int64,Ey::Int64)
	#=
        Create cuboid geometry
        inputs: Nx, number of samples x direction
        Ny, number of samples y direction
        Nz, number of samples z direction
        output: G, 3D tensor containing 1 inside domain and 0 outside
	=#

	if(Ex>Nx || Ey>Ny) error("edge outside room!") end

	G = ones(Bool,Nx,Ny,Nz)                 # creates the tensor with 1 where domain is present 
	G[end-Ex:end,end-Ey:end,:] = 0     # this change geometry from pure cuboid 

	G2 = zeros(Bool,Nx+2,Ny+2,Nz+2)# puts a shell of zeros around the tensor
	G2[2:Nx+1, 2:Ny+1, 2:Nz+1] = G
	
	return G2
end


fun_name(f::LShapedRoom) = "L-Shaped room"
fun_dim( f::LShapedRoom) = @sprintf("Lx x Ly x Lz = %.2f x %.2f x %.2f m^3", f.Lx, f.Ly, f.Lz)
fun_edge(f::LShapedRoom) = @sprintf("(x,y) = (%.2f, %.2f) m", f.xedge, f.yedge)
fun_Nxyz(f::LShapedRoom) = @sprintf("Nx x Ny x Nz = %d x %d x %d", f.Nx, f.Ny, f.Nz)
fun_ξ(   f::LShapedRoom) =
@sprintf("frequency independent   \n\n |  ξx1 |  ξx2 |  ξy1 |  ξy2 |  ξz1 | ξz2  | \n | % 02.1f| % 02.1f| % 02.1f| % 02.1f| % 02.1f| % 02.1f|\n" , f.ξ[1],f.ξ[2],f.ξ[3],f.ξ[4],f.ξ[5],f.ξ[6])
fun_X(  f::LShapedRoom) = @sprintf(" %.2f cm", f.env.X*100)
fun_Fs( f::LShapedRoom) = @sprintf(" %.2f kHz", f.env.Fs/1000)

function Base.show(io::IO, f::LShapedRoom)

	println(io, "geometry     : ", fun_name(f))
	println(io, "dimensions   : ", fun_dim(f))
	println(io, "edge at      : ", fun_edge(f))
	println(io, "samples      : ", fun_Nxyz(f))
	println(io, "ξ            : ", fun_ξ(f))
	println(io, "spatial step : ", fun_X(f))
	println(io, "sampling frq.: ", fun_Fs(f))
	
end
