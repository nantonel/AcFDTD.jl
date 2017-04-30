abstract AbstractGeometry   

import Base: size

include("getNodeProperty.jl")
include("cuboidRoom.jl")
#include("LShapedRoom.jl")

size(f::AbstractGeometry) = (f.Nx,f.Ny,f.Nz)

function Base.show(io::IO, f::AbstractGeometry)

	println(io, "geometry     : ", fun_name(f))
	println(io, "dimensions   : ", fun_dim(f))
	println(io, "samples      : ", fun_Nxyz(f))
	println(io, "ξ            : ", fun_ξ(f))
	println(io, "spatial step : ", fun_X(f))
	println(io, "sampling frq.: ", fun_Fs(f))
	
end
