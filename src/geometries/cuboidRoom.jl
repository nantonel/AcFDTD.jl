export CuboidRoom

"""
`cuboid room object`

Returns a `CuboidRoom` object containing the geometrical information of the room.

# Usage 

* CuboidRoom(Lx::Float64,Ly::Float64,Lz::Float64,ξ::Array{Float64,1}, env::FDTDEnv)
* CuboidRoom(Nx::Int64,Ny::Int64,Nz::Int64,ξ::Array{Float64,1}, env::FDTDEnv)
  * `ξ::Array{Float64}`: must be 6 element Array containg acoustic impedance of the 6 walls

"""
immutable CuboidRoom <: AbstractGeometry
	Lx::Float64   #x dimension
	Ly::Float64   #y dimension
	Lz::Float64   #z dimension

	Nx::Int64   #x-axis samples
	Ny::Int64   #y-axis samples
	Nz::Int64   #z-axis samples

	ξ::Array{Float64,1}   #Impedance
	env::FDTDEnv          #FDTD acoustic env.
end
	
function CuboidRoom(Lx::Float64,Ly::Float64,Lz::Float64,
		    ξ::Array{Float64,1}, env::FDTDEnv)

	if(any([Lx;Ly;Lz].< 0)) error("room dimensions L should be positive") end
	if(length(ξ)!= 6) error("length(ξ) must be 6") end

	Nx = round(Int64,Lx/env.X)
	Ny = round(Int64,Ly/env.X)
	Nz = round(Int64,Lz/env.X)

	CuboidRoom(env.X*Nx,env.X*Ny,env.X*Nz,Nx,Ny,Nz,ξ,env)
end

function CuboidRoom(Nx::Int64,Ny::Int64,Nz::Int64,
		    ξ::Array{Float64,1}, env::FDTDEnv)

	if(any([Nx;Ny;Nz].< 0)) error("room dimensions L should be positive") end
	if(length(ξ)!= 6) error("length(ξ) must be 6") end

	Lx = Nx*env.X
	Ly = Ny*env.X
	Lz = Nz*env.X

	CuboidRoom(Lx,Ly,Lz,Nx,Ny,Nz,ξ,env)
end

fun_name(f::CuboidRoom) = "cuboid room"
fun_dim(f::CuboidRoom) = @sprintf("Lx x Ly x Lz = %.2f x %.2f x %.2f m^3", f.Lx, f.Ly, f.Lz)
fun_Nxyz(f::CuboidRoom) = @sprintf("Nx x Ny x Nz = %d x %d x %d", f.Nx, f.Ny, f.Nz)
fun_ξ(f::CuboidRoom) =
@sprintf("frequency independent   \n\n |  ξx1 |  ξx2 |  ξy1 |  ξy2 |  ξz1 | ξz2  | \n | % 02.1f| % 02.1f| % 02.1f| % 02.1f| % 02.1f| % 02.1f|\n" , f.ξ[1],f.ξ[2],f.ξ[3],f.ξ[4],f.ξ[5],f.ξ[6])
fun_X(f::CuboidRoom) = @sprintf(" %.2f cm", f.env.X*100)
fun_Fs(f::CuboidRoom) = @sprintf(" %.2f kHz", f.env.Fs/1000)
