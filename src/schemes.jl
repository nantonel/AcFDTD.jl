export IISO,
       SLF,
       IWB
       
abstract FDTDScheme 

immutable IISO <: FDTDScheme
	λ::Float64              # Courant number
	a::Float64              # scheme parameter a
	b::Float64              # scheme parameter b
	c::Float64              # scheme parameter c
	name::AbstractString    # name of scheme
	IISO() = new(sqrt(3/4), 1/6, 0., 0., "Interp. Isotropic (IISO)")
end

immutable SLF <: FDTDScheme
	λ::Float64              # Courant number
	a::Float64              # scheme parameter a
	b::Float64              # scheme parameter b
	c::Float64              # scheme parameter c
	name::AbstractString    # name of scheme
	SLF() = new(sqrt(1/3), 0., 0., 0., "Standard Leapfrog Scheme (SLF)")
end

immutable IWB <: FDTDScheme
	λ::Float64              # Courant number
	a::Float64              # scheme parameter a
	b::Float64              # scheme parameter b
	c::Float64              # scheme parameter c
	name::AbstractString    # name of scheme
	IWB() = new(1., 1/4, 1/16, 0., "Interp. Wideband (IWB)")
end

function Base.show(io::IO, scm::FDTDScheme)
	println(io, "FDTD Scheme    : $(scm.name)")
end
